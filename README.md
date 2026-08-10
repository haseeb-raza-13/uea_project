# uea_project

Tools, pipelines, and outputs for a chromosomal rearrangement study on *Acinetobacter
baumannii*, using short-read and long-read whole-genome sequencing.

This repository holds the processing pipeline, quality-control results, read mapping
outputs, and analysis figures. Raw sequencing data and internal working documents are
not included here and are kept private; the tools in this repository are written to
work against that data using the directory layout described below.

## Repository layout

```
tools/       Pipeline scripts (bash entry points + Python drivers)
scripts/     Shared utilities for logging pipeline steps
qc_output/   FastQC and MultiQC reports, raw and adapter-trimmed reads
mapping/     Reference index, deduplicated alignments, and marker-frequency tables
analysis/    Downstream analysis scripts and publication-format figures
```

## Pipeline overview

The pipeline runs in a WSL/Linux environment using two conda environments:

- `seqqc`: FastQC, fastp, Trimmomatic, cutadapt, BBMap (BBDuk), MultiQC, Bowtie2,
  samtools, bedtools, and supporting Python packages.
- `rplots`: R with ggplot2, for the final figures.

Each stage below is a script in `tools/`, run from the project root. Most scripts have
a `.sh` entry point that activates the right environment and a `.py` driver
underneath.

1. **Quality control** (`tools/run_short_read_qc.sh`): runs FastQC on all short-read
   FASTQ files and aggregates the results with MultiQC. Accepts `--input-dir` and
   `--out-name` so the same script covers both raw and post-trimming data.

2. **Adapter trimming** (`tools/run_bbduk_trim.sh` / `run_bbduk_trim.py`): pairs R1/R2
   FASTQ files by sample and removes Illumina/Nextera adapter contamination with
   BBDuk, using BBMap's standard adapter-trimming preset.

3. **Origin/terminus localization** (`tools/find_ori_ter.sh` / `find_ori_ter.py`):
   locates the chromosomal origin of replication (oriC) and terminus (ter) on a
   reference genome using cumulative GC-skew analysis, independent of any externally
   supplied coordinates.

4. **Reference indexing** (`tools/build_reference_index.sh`): builds a Bowtie2 index
   for a reference genome.

5. **Read mapping** (`tools/run_mapping.sh` / `run_mapping.py`): aligns paired-end
   reads to the reference with Bowtie2, then removes PCR duplicates with samtools
   (sort, fixmate, markdup).

6. **Marker-frequency computation** (`tools/compute_marker_frequency.sh` /
   `compute_marker_frequency.py`): computes windowed read-start counts from the
   deduplicated alignments, normalizes them to the genome-wide median, and folds
   genomic position onto an oriC-centered coordinate (oriC = 0, terminus = plus or
   minus 1), following the marker-frequency-analysis method of Skovgaard et al.
   (2011, Genome Research).

7. **Sample metadata** (`tools/build_sample_metadata.sh` / `build_sample_metadata.py`):
   parses the sample-naming spreadsheet into a tidy CSV (strain, biological
   replicate, treatment, growth phase) used to group samples for plotting.

8. **Visualization** (`tools/plot_replication_profile.sh`, running
   `analysis/replication_profile/<strain>/plot_replication_profile.R`): generates
   replication-profile plots comparing biological replicates and treatment
   conditions, with the origin of replication centered and the terminus at the plot
   edges.

Every step above is also recorded automatically, with the exact command, tool
version, and outcome, in a project log kept outside this repository.

## Outputs included in this repository

- `qc_output/fastqc`, `qc_output/fastqc_trimmed`: per-file FastQC reports, before and
  after adapter trimming.
- `qc_output/multiqc`, `qc_output/multiqc_trimmed`: aggregated QC summaries.
- `mapping/reference/<strain>`: Bowtie2 index, reference FASTA, and the GC-skew
  diagnostic plot and derived oriC/ter coordinates.
- `mapping/bam/<strain>`: deduplicated, indexed BAM alignments (tracked with Git LFS).
- `mapping/coverage/<strain>`: windowed marker-frequency tables (TSV) and sample
  metadata.
- `analysis/replication_profile/<strain>`: the R plotting script and the resulting
  figures (PNG and PDF).

## Notes on scope

This repository currently covers strain **mff**. The pipeline is written to be
reused for additional strains by supplying their reference genome and running the
same sequence of tools.

Git LFS is required to fetch the BAM files in `mapping/bam/`. Run `git lfs install`
once, then clone or pull as usual.
