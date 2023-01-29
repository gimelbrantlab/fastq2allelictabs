# Preprocessing steps for alelle-specific expression on RNA-seq
Includes:
1. pseudoreference preparation: `pseudoreferences_creation/prepare_pseudoreference.py`
2. alignment: `STAR`
3. alelic reads resolving: `pseudoreferences_creation/prepare_pseudoreference.py`
4. counting reads per gene: `featureCounts`

For an example wrapper function for steps (1-3) see `fastq2allelicbams.sh`; for stats collection (like # of raw reads, # of aligned reads, spike-in reads proportion) see `fastq2allelicbams_stats.sh`; for step (4) see `allelicbams2genecounts.sh`. See `example` directory for sample butch table example, and Wiki page for more details and usecases.

![pic](https://github.com/gimelbrantlab/fastq2allelictabs/blob/main/schemes/ase-preprocessing-outline.png)
