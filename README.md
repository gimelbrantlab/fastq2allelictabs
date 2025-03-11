# Preprocessing protocol for alelle-specific expression on RNA-seq
Includes:
1. pseudoreference preparation: `pseudoreferences_creation/prepare_pseudoreference.py`
2. alignment: `STAR`
3. alelic reads resolving: `fastq_to_allelic_counts_tabs/alleleseparation.py`
4. counting reads per gene: `featureCounts`

See also: [controlFreq](https://github.com/gimelbrantlab/controlFreq), R-package for calculating overdispersion in RNA-seq samples, in presence of technical replication or spike-ins.
### Installation
`docker pull sasha/ase`

### Description
For an example wrapper function for steps (1-3) see `fastq2allelicbams.sh`; for stats collection (like # of raw reads, # of aligned reads, spike-in reads proportion) see `fastq2allelicbams_stats.sh`; for step (4) see `allelicbams2genecounts.sh`. See `example` directory for sample butch table example, and Wiki page for more details and usecases, motivation of pipeline choice, and QC.

![pic](https://github.com/gimelbrantlab/fastq2allelictabs/blob/main/schemes/ase-preprocessing-outline.png)

Note: step (1) is the same as in [ASEReadCounter*](https://github.com/gimelbrantlab/ASEReadCounter_star) (see [Wiki](https://github.com/gimelbrantlab/ASEReadCounter_star/wiki/1.-Reference-preparation#preparing-files-for-using-as-reference)), when (2-4) have been evolved.

Scheme is made in BioRender.
