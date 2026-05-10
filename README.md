# ONT Long-Read SV/CNV Pipeline

A reproducible pipeline for calling structural variants (SVs) and copy number variants (CNVs) from Oxford Nanopore Technology (ONT) long-read sequencing data.

## Workflow

```
Raw FASTQ files (ONT basecalled reads)
        ↓
Concatenation (cat)
        ↓
Read Length Filtering — min 1000 bp (filtlong)
        ↓
Alignment to Reference Genome (minimap2)
        ↓
BAM Conversion, Sorting & Indexing (samtools)
        ↓
Alignment QC — gates on minimum alignment rate (samtools flagstat)
        ↓
Structural Variant Calling (Sniffles2)
        ↓
VCF Summary Statistics (bcftools)
        ↓
output_structural_variants.vcf + output_vcf_stats.txt
```

## Repository Contents

| File | Description |
|------|-------------|
| `ont_sv_pipeline.sh` | End-to-end pipeline script |

## Dependencies

All tools are installable via conda. Run inside a dedicated conda environment.

| Tool | Purpose | Install |
|------|---------|---------|
| [minimap2](https://github.com/lh3/minimap2) | Long-read alignment | `conda install bioconda::minimap2` |
| [samtools](http://www.htslib.org/) | SAM/BAM conversion, sorting, QC | `conda install bioconda::samtools` |
| [filtlong](https://github.com/rrwick/Filtlong) | Read length filtering | `conda install bioconda::filtlong` |
| [Sniffles2](https://github.com/fritzsedlazeck/Sniffles) | SV/CNV calling from long reads | `conda install bioconda::sniffles` |
| [bcftools](https://samtools.github.io/bcftools/) | VCF statistics | `conda install bioconda::bcftools` |
| [moreutils](https://joeyh.name/code/moreutils/) | Timestamped logging (`ts`) | `conda install conda-forge::moreutils` |

## Environment Setup

```bash
# Add channels (run once)
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Create a dedicated environment
conda create -n ont_sv minimap2 samtools filtlong sniffles bcftools moreutils
conda activate ont_sv
```

## Reference Genome

Download and unzip the hg38 reference genome before running:

```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
```

## Usage

**1. Edit the configuration block** at the top of the script:

```bash
FASTQ_DIR="."          # directory containing your .fastq or .fastq.gz files
REFERENCE="hg38.fa"    # path to reference genome FASTA
THREADS=16             # adjust to your available CPU cores
MIN_READ_LENGTH=1000   # reads shorter than this are discarded before alignment
MIN_ALIGNMENT_RATE=70  # pipeline exits if alignment rate falls below this (%)
MIN_DISK_GB=50         # pipeline exits if available disk space is below this (GB)
OUT_PREFIX="output"    # prefix for all output files
```

**2. Run the pipeline:**

```bash
chmod +x ont_sv_pipeline.sh
./ont_sv_pipeline.sh
```

**3. To run in the background** (recommended for long jobs — persists after terminal closes):

```bash
nohup bash ont_sv_pipeline.sh &

# Monitor progress via the timestamped log file
tail -f pipeline_YYYYMMDD_HHMMSS.log
```

## Output Files

| File | Description |
|------|-------------|
| `*_concatenated.fastq` | All input reads merged into a single file |
| `*_filtered.fastq` | Reads passing minimum length filter |
| `*_genomic_sorted.bam` | Coordinate-sorted BAM alignment |
| `*_genomic_sorted.bam.bai` | BAM index |
| `*_flagstat.txt` | Alignment QC statistics (mapped rate, read counts) |
| `*_structural_variants.vcf` | Called SVs/CNVs in VCF format |
| `*_vcf_stats.txt` | Per-SV-type breakdown (DEL, INS, DUP, INV, BND) |
| `pipeline_TIMESTAMP.log` | Timestamped log of the full run |

## Key Design Decisions

**Read length filtering (`filtlong --min_length 1000`)**
Very short ONT reads (< 1 kb) align ambiguously and inflate false positive SV calls. Filtering before alignment improves both alignment rate and SV precision.

**`--secondary=no` in minimap2**
Secondary alignments are not used by Sniffles and can inflate BAM file size by 2–3×. Suppressing them speeds up sorting and SV calling with no loss of information.

**Alignment rate gate**
The pipeline will exit before SV calling if the overall alignment rate falls below `MIN_ALIGNMENT_RATE`. This prevents Sniffles from running on a low-quality or mismatched alignment and producing unreliable variant calls.

**Sniffles2 syntax**
This pipeline uses the Sniffles2 (v2.x) command-line interface (`--input`, `--vcf`). If you are running Sniffles v1.x, replace the Sniffles call with:
```bash
sniffles -m "$SORTED_BAM" -v "$VCF"
```

**Intermediate file cleanup**
The unsorted SAM and BAM are removed after the sorted BAM is produced. Long-read SAM files are large; keeping them would roughly triple disk usage with no downstream benefit.
