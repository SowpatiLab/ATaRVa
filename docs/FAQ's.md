<details>
<summary><b>Whats the difference between --read-wise and --loci-wise mode?</b></summary>

By default, ATaRVa processes alignment files (BAM/CRAM) in `--read-wise` mode. In this mode, reads are processed sequentially in coordinate-sorted order, and all tandem repeat (TR) loci spanned by a read are analyzed together. Since long reads often span multiple nearby TR regions, this approach efficiently reuses read information across loci. The `--read-wise` mode is optimized for genome-wide analyses involving millions of TR loci and significantly reduces runtime compared to locus-by-locus processing.

ATaRVa also provides a `--loci-wise`, which follows the conventional approach of processing TR loci individually. This mode is recommended for targeted analyses involving a relatively small number of loci (e.g., fewer than 500 regions), particularly when the loci are distributed across different genomic locations. For such use cases, `--loci-wise` can provide faster execution by avoiding the overhead of scanning reads genome-wide.

</details>


