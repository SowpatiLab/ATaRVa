#!/usr/bin/env python3
"""
Merge many single-sample tandem-repeat VCFs into one multi-sample VCF,
splitting work into chunks with ~equal numbers of loci.

Strategy
  * inputs: bgzipped + tabix-indexed single-sample VCFs from the SAME locus catalog
  * catalog: catalog file
  * a reference fasta file is required with the catalog is provided 
  * a reference VCF is scanned once to place chunk boundaries every K loci;
    chunks may span contigs (list of (contig, lo, hi) segments, lo <= POS < hi)
  * each worker opens all inputs once, k-way merges its segments with a heap,
    and writes a BGZF-compressed part file
  * parts are concatenated byte-wise (valid BGZF) in chunk order and tabix-indexed
  * locus key = (POS, END, REF, TRID/ID); ALT alleles unioned by sequence, GT remapped
  * missing loci -> ./.

Requirements: pysam.   For thousands of inputs:  ulimit -n 65536
Memory note: each open TabixFile holds its index in RAM -> ~ threads x samples x index size.

Usage:
  python merge_tr_vcfs_chunked.py -l vcf_list.txt -o merged.vcf.gz -t 16
  python merge_tr_vcfs_chunked.py -l list.txt -o out.vcf.gz -t 16 --chunks-per-thread 1 --ref-vcf best.vcf.gz
"""
import argparse
import heapq
import os
import re
import shutil
import sys
import tempfile
from array import array
from multiprocessing import Pool
from collections import Counter

import pysam

END_RE       = re.compile(r"(?:^|;)END=(\d+)")
START_RE     = re.compile(r"(?:^|;)START=(\d+)")
MOTIF_RE     = re.compile(r"(?:^|;)MOTIF=([^;]+)")
TRID_RE      = re.compile(r"(?:^|;)ID=([^;]+)")
GT_SPLIT     = re.compile(r"([/|])")
DEFAULT_KEYS = ['GT', 'AL', 'CN', 'AR', 'SD', 'DP', 'SN', 'SQ', 'MA', 'MR', 'DS', 'MV', 'OL', 'FV', 'PS']
FLUSH_EVERY  = 100

_H = []  # per-process tabix handles, set by init_worker


def merge_parser(subparsers):
    """
    Merge ATaRVa VCF files over specified regions.

    :param subparsers: ArgumentParser subparsers object
    :return: None
    """

    parser = subparsers.add_parser("merge", help="merging multiple ATaRVa VCF files over specified regions", description="Merge ATaRVa VCF files",
                                   formatter_class=argparse.RawDescriptionHelpFormatter)
    parser._action_groups.pop()

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--vcfs', nargs='+', required=True, metavar="<FILE>", help="ATaRVa output VCF files to be merged.")
    required.add_argument('-r', '--regions',         required=True, metavar='<FILE>', help="bgzipped + tabix-indexed regions file. Prepare: sort with bedtools"
                                                                                           "→ bgzip → tabix index. Include with reference fasta file.")
    required.add_argument('-f', '--fasta',           required=True, metavar='<FILE>', help="input reference fasta file. requed if --regions is provided. "
                                                                                           "Prepare: samtools faidx.")
    required.add_argument('-o', '--out',  type=str,  required=True, metavar='<FILE>', help="name of the output file, output is in vcf format.")

    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--contigs',       nargs='+', help="Contigs to get merged [chr1 chr12 chr22 ..]. If not mentioned every contigs in the region file will be merged.")
    optional.add_argument('--format-keys',   nargs='+', metavar='<STR>', default=','.join(DEFAULT_KEYS), help=f"FORMAT keys to keep in output. GT is by default part of the"
                                                                                                               " output. (default: {','.join(DEFAULT_KEYS)})")
    optional.add_argument('-t', '--threads', type=int,  metavar='<INT>', default=1,    help="number of threads. [default: 1]")

    if (len(sys.argv) == 2) and (sys.argv[1] == 'merge'):
        parser.print_help()
        sys.exit()

    parser.set_defaults(func=run_merge)


def init_worker(paths):
    global _H
    _H = [pysam.TabixFile(p) for p in paths]


def read_header(path):
    tbx     = pysam.TabixFile(path)
    hdr     = list(tbx.header)
    contigs = list(tbx.contigs)
    tbx.close()
    return hdr, contigs


def sample_name(header_lines):
    chrom_line = [h for h in header_lines if h.startswith("#CHROM")][0]
    return chrom_line.rstrip("\n").split("\t")[9]


def vcf_iter(tbx, contig, lo, hi):
    """Yield (key, fields) for records with lo <= POS < hi (hi=None -> open end)."""
    try:
        it = tbx.fetch(contig, lo, hi + 1)
    except ValueError:  # contig absent in this file
        return
    for line in it:
        f = line.split("\t", 10)
        pos = int(f[1])
        start  = START_RE.search(f[7])
        start  = int(start.group(1)) if start else pos
        if start < lo: continue
        if start > hi: break
        end    = END_RE.search(f[7])
        end    = int(end.group(1)) if end else pos + len(f[3])
        motif  = MOTIF_RE.search(f[7])
        motif  = motif.group(1) if motif else f[3]
        _id    = TRID_RE.search(f[7])
        _id    = _id.group(1) if _id else f"{f[0]}-{start}-{end}-{motif}"
        yield (start, end, motif, _id), f


def catalog_iter(tbx, contig, lo, hi, fasta):
    """Yield (key, fields) for records with lo <= POS < hi (hi=None -> open end)."""
    try:
        it = tbx.fetch(contig, lo, hi + 1)
    except ValueError:  # contig absent in this file
        return
    for line in it:
        f = line.split("\t")
        start = int(f[1])
        end   = int(f[2])
        motif = f[3]
        locus_id = f'{f[0]}-{start}-{end}-{motif}'

        if start < lo: continue
        if start > hi: break

        ref_seq = fasta.fetch(f[0], start, end)
        f.append("CATALOG") # add a keyword to distinguish catalog records from sample records
        yield (start, end, motif, locus_id), f 


def remap_gt(gt, amap):
    parts = GT_SPLIT.split(gt)
    for i in range(0, len(parts), 2):
        a = parts[i]
        if a not in (".", "0", ""):
            parts[i] = str(amap[int(a)])
    return "".join(parts)

def count_alleles(gt, counter):
    parts = GT_SPLIT.split(gt)
    for i in range(0, len(parts), 2):
        a = parts[i]
        if a not in (".", "0", ""):
            counter[int(a)] += 1

def allele_len(allele):
    if allele == "<DEL>": return 0
    else: return len(allele)

def alts_onlen(alts):
    """
    Calculate the alternates based on their lengths.
    """

    alts_len = {}
    for a in alts:
        alen = allele_len(a)
        if alen not in alts_len:
            alts_len[alen] = [a]
        else:
            alts_len[alen].append(a)

    alens = sorted(alts_len.keys())

    return alens, alts_len


def tamatr(args):

    cid, catalog, ref_file, segments, keys, tmpdir = args
    tbx   = pysam.TabixFile(catalog)
    fasta = pysam.FastaFile(ref_file)

    n = len(_H)
    missing_str = ":".join("./." if k == "GT" else "." for k in keys)
    fmt_cache = {}
    out_path = os.path.join(tmpdir, f"chunk_{cid:06d}.part.gz")
    buf = []

    with pysam.libcbgzf.BGZFile(out_path, "wb") as out:
        for segment in segments:
            contig, lo, hi = segment
            lo = lo[0]
            hi = hi[0]
            iters = [catalog_iter(tbx, contig, lo, hi, fasta)] + [vcf_iter(h, contig, lo, hi) for h in _H]
            heap = []
            for i, it in enumerate(iters):
                nxt = next(it, None)
                if nxt:
                    heap.append((nxt[0], i, nxt[1]))
            heapq.heapify(heap)

            while heap:
                key = heap[0][0]
                group = []
                locus = None
                while heap and heap[0][0] == key:
                    _, i, f = heapq.heappop(heap)
                    if f[-1] == "CATALOG": locus = f
                    else: group.append((i, f))
                    nxt = next(iters[i], None)
                    if nxt:
                        heapq.heappush(heap, (nxt[0], i, nxt[1]))

                first = None
                alts  = []
                cols  = [missing_str] * n
                if len(group) > 0:
                    first = group[0][1]
                    alts  = set()
                    for _, f in group:
                        if f[4] != ".":
                            alts.update(f[4].split(","))
                    alts       = sorted(alts, key=lambda s: (allele_len(s), s))
                    alt_counts = {}
                    alt_idx = {a: j + 1 for j, a in enumerate(alts)}
                    # if args.length_based:
                    #     alens, alts_len = alts_onlen(alts)
                    #     ref_len = len(first[3])
                    #     if ref_len in alens: alens.remove(ref_len)
                    #     alt_idx = {a: j + 1 for j, a in enumerate(alens)}  # remap GT to the shortest ALT alleles
                    #     alt_idx[ref_len] = 0
                    #     for x in alt_idx: alt_counts[x] = Counter()

                    cols = [missing_str] * n
                    for i, f in group:
                        format  = f[8]
                        sel_keys  = fmt_cache.get(format)
                        if sel_keys is None:
                            fmt_keys = format.split(":")
                            sel_keys = [fmt_keys.index(k) if k in fmt_keys else -1 for k in keys]
                            fmt_cache[format] = sel_keys
                        sample   = f[9].rstrip("\n").split(":")
                        loc_alts = f[4].split(",") if f[4] != "." else []
                        # if args.length_based:
                        #     amap = {j + 1: alt_idx[len(a)] for j, a in enumerate(loc_alts)}
                        # else:
                        amap = {j + 1: alt_idx[a] for j, a in enumerate(loc_alts)}
                        row = []
                        for k, s in zip(keys, sel_keys):
                            v = sample[s] if 0 <= s < len(sample) else "."
                            if k == "GT" and v != ".":
                                v = remap_gt(v, amap)
                            row.append(v)
                        cols[i-1] = ":".join(row)

                elif len(group) == 0:
                    ref_seq = key[-1]
                    info = f'AC=0;AN=0;MOTIF={key[-1]};START={key[0]};END={key[1]};ID=.;REFCN={(key[1]-key[0])//len(ref_seq)}'
                    first = [locus[0], str(int(locus[1])+1), '.', ref_seq, '0', 'FAIL', '.', info]

                entry = first[:4] + [",".join(alts) or "."] + [".", "PASS", first[7], ":".join(keys)] + cols
                buf.append( "\t".join(entry) + "\n")
                if len(buf) >= FLUSH_EVERY:
                    out.write("".join(buf).encode())
                    buf.clear()
        if buf:
            out.write("".join(buf).encode())

    tbx.close()
    fasta.close()
    return cid, out_path


def splitfile_threads(tbx, contigs, threads):
    """
    split the input region file into per-thread chunks for parallel processing.

    :param tbx:        pysam Tabixfile object
    :param contigs:    list of contigs to process
    :param threads:    number of parallel threads
    :return:           list of region chunk tuples, one per thread
    """

    tbx = pysam.TabixFile(tbx)
    total_loci = 0
    for contig in contigs:
        for _ in tbx.fetch(contig):
            total_loci += 1

    split_size  = total_loci // threads
    if split_size == 0:
        split_size = total_loci
    fetcher     = []
    line_count  = 0
    cur_chunk   = []
    chrom       = None
    start_coord = None

    _start, _end = None, None

    for contig in contigs:
        first_in_contig = True

        fields = None
        for row in tbx.fetch(contig):
            fields = row.split('\t')
            line_count += 1
            if first_in_contig and _start == int(fields[1]) and _end == int(fields[2]):
                continue  # skip duplicate rows

            if first_in_contig:
                chrom       = fields[0]
                start_coord = (int(fields[1]), int(fields[2]))
                first_in_contig = False

            if len(fetcher) < threads - 1 and line_count % split_size == 0:
                end_coord = (int(fields[1]), int(fields[2]))
                cur_chunk.append([chrom, start_coord, end_coord])
                fetcher.append(tuple(cur_chunk))
                cur_chunk       = []
                line_count      = 0
                first_in_contig = True

            _start, _end = int(fields[1]), int(fields[2])

        if not first_in_contig:
            end_coord = (int(fields[1]), int(fields[2]))
            cur_chunk.append([chrom, start_coord, end_coord])

    fetcher.append(tuple(cur_chunk)) # pad with empty chunks if needed
    return fetcher, total_loci


def run_merge(args):
    
    paths = list(args.vcfs or [])
    if not paths: sys.exit("Input VCFs are not provided!")

    keys = args.format_keys.split(",")
    if keys[0] != "GT":
        keys = ['GT'] + keys
        sys.exit("--format-keys must start with GT")

    catalog = args.regions
    fasta   = args.fasta

    hdr0 = None
    contigs, seen, samples = [], set(), []
    for p in paths:
        h, cs = read_header(p)
        if hdr0 is None: hdr0 = h
        samples.append(sample_name(h))
        for c in cs:
            if c not in seen:
                seen.add(c)
                contigs.append(c)
    if args.contigs:
        contigs = args.contigs.split(",")
    if len(set(samples)) != len(samples):
        sys.exit("Error: Duplicate sample names among inputs")

    avail  = {l.split("ID=")[1].split(",")[0]: l for l in hdr0 if l.startswith("##FORMAT=<ID=")}
    header = [l.rstrip("\n") for l in hdr0
              if l.startswith("##") and not l.startswith("##FORMAT")
              and not l.startswith("##bcftools") and not l.startswith("##source=")]
    header += [avail[k].rstrip("\n") for k in keys if k in avail]
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))

    fetcher, total_loci = splitfile_threads(catalog, contigs, args.threads)

    print(f"Reference loci: {total_loci}; threads: {args.threads} (~{-(-total_loci // len(fetcher))} loci each)", file=sys.stderr)

    tmpdir = tempfile.mkdtemp(prefix="trmerge_", dir=os.path.dirname(os.path.abspath(args.out)))
    try:
        jobs = [(i, catalog, fasta, segs, keys, tmpdir) for i, segs in enumerate(fetcher)]
        with Pool(min(args.threads, len(jobs)), initializer=init_worker, initargs=(paths,)) as pool:
            results = dict(pool.imap_unordered(tamatr, jobs))

        if not args.out.endswith("gz"):
            args.out += ".gz"

        with pysam.libcbgzf.BGZFile(args.out, "wb") as out:
            out.write(("\n".join(header) + "\n").encode())
        with open(args.out, "ab") as out:
            for i in range(len(fetcher)):  # chunk order == contig/position order
                with open(results[i], "rb") as part:
                    shutil.copyfileobj(part, out)
        pysam.tabix_index(args.out, preset="vcf", force=True)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
    print(f"Merged {len(paths)} samples -> {args.out}", file=sys.stderr)
