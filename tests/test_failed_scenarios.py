
import sys
from pathlib import Path

if __name__ == "__main__":                      
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ATARVA import decompose as mt
from ATARVA.decompose import motif_decomposition, get_repeat_components

SCENARIOS = [
    ('A', 'CTCCATCTCCATCTCCATCTCCATCTCCATCTCCATTTGGGGGCTA', 6, '(CTCCAT)6-TT-(G)5-CTA'),
    ('A', 'AAAAAAATTTTTTTC', 2, '(A)7-(T)7-C'),
    ('A', 'ACAACAAGACCCCCT', 3, '(ACA)2-AGA-(C)5-T'),
    ('A', 'ACCACCACCACCACCGTTTTTTTTTTGC', 7, '(ACC)5-G-(T)10-GC'),
    ('A', 'ATCTGAATCTGAATCTGAATCTGAATCTGAATCTGAATCTGACGAGAGAGAGAGAGACCTTACTTA', 5, '(ATCTGA)7-C-(GA)7-C-(CTTA)2'),
    ('B', 'AAAA', 3, '(A)4'),
    ('B', 'AAAAAA', 2, '(A)6'),
    ('B', 'AAAAAAAAAAAA', 4, '(A)12'),
    ('B', 'AAAAAAAAAAAA', 6, '(A)12'),
    ('C', 'AAAAAAAAAACTA', 3, '(A)9-ACTA'),
    ('C', 'AAAAACCAAAAGTG', 2, '(A)5-CC-(A)4-GTG'),
    ('C', 'AAAAAACTT', 2, '(A)6-CTT'),
    ('C', 'CCCCCCCCCACCCCCCCCCCCCCCAAAA', 2, '(C)9-A-(C)14-(A)4'),
    ('D', 'GAGAGAGAGAGAGAACCTGTACCTGTACCTGTACCTGTTG', 4, '(GA)7-(ACCTGT)4-TG'),
    ('D', 'TCTCTCTCTCTCTCTCCGCGACCACGACCACGACCACGACCAC', 2, '(TC)8-CG-(CGACCA)4-C'),
    ('D', 'ATTTTCATTTTCATTTTCTGTGTGTGTGTGTGTGAGCTTAGCTTA', 8, '(ATTTTC)3-(TG)8-(AGCTT)2-A'),
    ('E', 'GAGTTCGAGTTCGAGTTCCGCGCGCGCGCGTC', 7, '(GAGTTC)3-(CG)6-TC'),
    ('E', 'AGTAAGAGTAAGAGTAAGGCACACACACACACAATA', 2, '(AGTAAG)3-G-(CA)7-ATA'),
    ('E', 'ATACCCATACCCATACCCGTCTCTCTCTCTCTCCGC', 5, '(ATACCC)3-G-(TC)7-CGC'),
    ('E', 'CGAGAACGAGAACGAGAATCTCTCTCTCTCTCGGGTTA', 7, '(CGAGAA)3-(TC)7-GGGTTA'),
    ('OK', 'CAGCAGCAGTTCAGCAG', 3, '(CAG)3-TT-(CAG)2'),
    ('OK', 'ACGTTGCA', 2, 'ACG-(T)2-GCA'),
    ('OK', 'ATATATATAT', 2, '(AT)5'),
    ('OK', 'TATCTGATATCTGATATCTGATATCTGATATCTGATATCTGAATCATCATCATC', 7, '(TATCTGA)6-(ATC)4'),
]


def expand(decomp):
    out = []
    for part in decomp.split("-"):
        if part.startswith("("):
            motif, n = get_repeat_components(part)
            out.append(motif * n)
        else:
            out.append(part)
    return "".join(out)


def is_atomic(motif):
    return motif not in (motif + motif)[1:-1]


def analyse(seq, size):
    mt.cyclic_motif_registry.clear()
    decomp, non_rep = motif_decomposition(seq, size)
    blocks = [get_repeat_components(p)[0] for p in decomp.split("-") if p.startswith("(")]
    return dict(decomp=decomp, non_rep=non_rep,
                roundtrip=expand(decomp) == seq,
                atomic=all(is_atomic(m) for m in blocks))


def describe(cat, seq, size, expected=None):
    r = analyse(seq, size)
    lines = [f"[{cat}] motif_size={size}  length={len(seq)}",
             f"    sequence      : {seq}",
             f"    decomposition : {r['decomp']}   (non_repeat={r['non_rep']})",
             f"    expands back to the sequence: {'yes' if r['roundtrip'] else 'NO  -> ' + expand(r['decomp'])}",
             f"    all blocks atomic           : {'yes' if r['atomic'] else 'NO'}"]
    if expected is not None:
        lines.append(f"    matches expected            : "
                     f"{'yes' if r['decomp'] == expected else 'NO  -> expected ' + expected}")
    return "\n".join(lines), r



try:
    import pytest
except ImportError:                              
    pytest = None

if pytest is not None:
    IDS = [f"{c}{i}-size{z}" for i, (c, s, z, e) in enumerate(SCENARIOS)]

    @pytest.fixture(autouse=True)
    def clean_registry():
        mt.cyclic_motif_registry.clear()

    @pytest.mark.parametrize("cat,seq,size,expected", SCENARIOS, ids=IDS)
    def test_scenario(cat, seq, size, expected):
        text, r = describe(cat, seq, size, expected)
        print("\n" + text)
        assert r["roundtrip"], f"expands to {expand(r['decomp'])}"
        assert r["atomic"], f"non-atomic block in {r['decomp']}"
        assert r["decomp"] == expected



if __name__ == "__main__":
    if len(sys.argv) == 3:                       # test a new sequence 
        print(describe("custom", sys.argv[1].upper(), int(sys.argv[2]))[0])
        sys.exit(0)
    bad = 0
    for cat, seq, size, expected in SCENARIOS:
        text, r = describe(cat, seq, size, expected)
        print(text + "\n")
        bad += not (r["roundtrip"] and r["atomic"] and r["decomp"] == expected)
    print(f"{len(SCENARIOS) - bad} of {len(SCENARIOS)} scenarios OK")