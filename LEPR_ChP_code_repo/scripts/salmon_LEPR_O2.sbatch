#!/bin/bash
#SBATCH --job-name=lepr_salmon
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=lepr_salmon_%j.out
#SBATCH --error=lepr_salmon_%j.err
# =====================================================================
# TOTAL LEPR EXPRESSION — Salmon, full human transcriptome (GENCODE v44)
# Harvard O2 (x86_64; salmon runs natively here).
# GSE228458 human choroid plexus, 13 donors (5 control / 4 AD / 4 DM1).
#
# Submit:   sbatch salmon_LEPR_O2.sbatch
# Needs internet on compute node? O2 'short' has it for ftp/https. If your
# node lacks internet, pre-download FASTQ+refs on a login node first (see README).
# =====================================================================
set -euo pipefail

# ---- salmon: O2 has NO salmon module. Use precompiled x86_64 binary in scratch.
#      Download once on a login node (see install block in README_O2.md):
#        wget .../salmon-1.10.0_linux_x86_64.tar.gz && tar xzf ...
export PATH="/n/scratch/users/s/sak2343/lepr_salmon/salmon-latest_linux_x86_64/bin:$PATH"
module load gcc/9.2.0 2>/dev/null || true   # libstdc++ etc; harmless if absent
command -v salmon >/dev/null || { echo "ERROR: salmon not on PATH — download the binary to scratch first"; exit 1; }
salmon --version

# ---- workspace = submission dir on scratch (scripts + sample_map.csv live here) ----
WORK=${WORK:-$PWD}; THREADS=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$WORK"/{ref,idx,quant,fastq}; cd "$WORK"

# =====================================================================
# 1) REFERENCE — GENCODE v44 (GRCh38.p14); decoy-aware
# =====================================================================
cd ref
[ -f gencode.v44.transcripts.fa.gz ]      || wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
[ -f GRCh38.primary_assembly.genome.fa.gz ] || wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
[ -f gencode.v44.annotation.gtf.gz ]      || wget -q https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
zcat GRCh38.primary_assembly.genome.fa.gz | grep "^>" | cut -d" " -f1 | sed 's/>//' > decoys.txt
cat gencode.v44.transcripts.fa.gz GRCh38.primary_assembly.genome.fa.gz > gentrome.fa.gz
zcat gencode.v44.annotation.gtf.gz | awk -F'\t' '$3=="transcript"{
  match($9,/transcript_id "([^"]+)"/,t); match($9,/gene_id "([^"]+)"/,g);
  match($9,/gene_name "([^"]+)"/,n); print t[1]"\t"g[1]"\t"n[1]}' > ../tx2gene.tsv
cd ..

# 2) INDEX  (reuse if already built: info.json is salmon's completion marker)
[ -f idx/info.json ] || salmon index -t ref/gentrome.fa.gz -d ref/decoys.txt -i idx -k 31 -p $THREADS

# =====================================================================
# 3) SAMPLES — download R1+R2 from ENA, quant each
# =====================================================================
declare -A RUN=( [ctrl1]=SRR24003390 [ctrl2]=SRR24003388 [ctrl3]=SRR24003385 [ctrl4]=SRR24003383 [ctrl5]=SRR24003378
  [AD1]=SRR24003376 [AD2]=SRR24003374 [AD3]=SRR24003372 [AD4]=SRR24003370
  [DM1_1]=SRR24003392 [DM1_2]=SRR24003381 [DM1_3]=SRR24003377 [DM1_4]=SRR24003393 )
ena(){ local s=$1; echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${s:0:6}/0${s: -2}/${s}"; }

# ---- Phase 1: PARALLEL download all 26 mates (6 concurrent; wget -c resumes
#      partials and skips already-complete files). Sequential ENA pulls were the
#      bottleneck (~2 MB/s each); 6 streams multiply throughput. ----
: > urls.txt
for nm in ctrl1 ctrl2 ctrl3 ctrl4 ctrl5 AD1 AD2 AD3 AD4 DM1_1 DM1_2 DM1_3 DM1_4; do
  s=${RUN[$nm]}; b=$(ena $s)
  echo "$b/${s}_1.fastq.gz"; echo "$b/${s}_2.fastq.gz"
done > urls.txt
echo ">>> downloading $(wc -l < urls.txt) files, 6 parallel"
xargs -P 6 -n 1 wget -c -q -P fastq/ < urls.txt || true
echo ">>> download phase done"; ls -la fastq/

# ---- Phase 2: quant each sample (sequential; reuses local fastq). Idempotent:
#      skips samples already quantified, single-end (R2) fallback if a mate is missing. ----
for nm in ctrl1 ctrl2 ctrl3 ctrl4 ctrl5 AD1 AD2 AD3 AD4 DM1_1 DM1_2 DM1_3 DM1_4; do
  s=${RUN[$nm]}
  if [ -f quant/${nm}/quant.sf ]; then echo ">>> $nm already quantified, skip"; continue; fi
  echo ">>> quant $nm ($s)"
  if [ -s fastq/${s}_1.fastq.gz ] && [ -s fastq/${s}_2.fastq.gz ]; then
    salmon quant -i idx -l A -1 fastq/${s}_1.fastq.gz -2 fastq/${s}_2.fastq.gz \
      -p $THREADS --validateMappings --gcBias -o quant/${nm}
  elif [ -s fastq/${s}_2.fastq.gz ]; then  # single-end fallback (R2 only) — less ideal
    salmon quant -i idx -l A -r fastq/${s}_2.fastq.gz \
      -p $THREADS --validateMappings -o quant/${nm}
  else
    echo "!!! $nm ($s): no fastq present, skipped"
  fi
done

# 4) extract LEPR gene + transcript TPM/counts -> CSV
python3 extract_lepr.py
echo "DONE. -> LEPR_expression_FINAL.csv  (+ quant/*/quant.sf)"
