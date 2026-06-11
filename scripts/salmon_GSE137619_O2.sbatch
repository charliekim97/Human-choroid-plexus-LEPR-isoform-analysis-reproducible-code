#!/bin/bash
#SBATCH --job-name=lepr_gse137619
#SBATCH --partition=medium
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=gse137619_%j.out
#SBATCH --error=gse137619_%j.err
# =====================================================================
# INDEPENDENT REPLICATION — GSE137619 human choroid plexus (6 control + 6 MS)
# Same pipeline as GSE228458: REUSES that decoy-aware GENCODE v44 salmon index,
# same salmon 1.10.0 binary, same extract_lepr.py (same transcript IDs).
# Paired-end (verified). Two read-outs, kept separate:
#   - expression : salmon quant -> LEPR_expression_GSE137619.csv
#   - usage      : junction reads (seqkit 30M seed42 + grep) -> lepr_junction_GSE137619.csv
# Submit from a fresh dir, e.g. /n/scratch/users/s/sak2343/lepr_GSE137619 :
#   sbatch salmon_GSE137619_O2.sbatch
# =====================================================================
set -euo pipefail

OLD=/n/scratch/users/s/sak2343/lepr_salmon            # GSE228458 working dir (reuse assets)
IDX=$OLD/idx                                          # reuse decoy-aware index
export PATH="$OLD/salmon-latest_linux_x86_64/bin:$PATH"
module load gcc/9.2.0 2>/dev/null || true
command -v salmon >/dev/null && salmon --version
[ -f "$IDX/info.json" ] || { echo "ERROR: reuse index not found at $IDX (was GSE228458 run purged?)"; exit 1; }

WORK=${WORK:-$PWD}; THREADS=${SLURM_CPUS_PER_TASK:-12}
mkdir -p "$WORK"/{quant,fastq}; cd "$WORK"
cp "$OLD/extract_lepr.py" .                           # identical isoform annotation

# ---- seqkit binary for junction subsampling (O2 has no seqkit module) ----
if ! command -v seqkit >/dev/null 2>&1; then
  [ -x ./seqkit ] || { wget -q https://github.com/shenwei356/seqkit/releases/download/v2.8.2/seqkit_linux_amd64.tar.gz && tar xzf seqkit_linux_amd64.tar.gz; }
  export PATH="$PWD:$PATH"
fi
seqkit version 2>/dev/null || true

# ---- samples (donor -> SRR ; R2 read_count for 30M subsample proportion ; group) ----
declare -A RUN=( [ctrl1]=SRR10134643 [ctrl2]=SRR10134644 [ctrl3]=SRR10134645 [ctrl4]=SRR10134646 [ctrl5]=SRR10134647 [ctrl6]=SRR10134648
  [MS1]=SRR10134649 [MS2]=SRR10134650 [MS3]=SRR10134651 [MS4]=SRR10134652 [MS5]=SRR10134653 [MS6]=SRR10134654 )
declare -A NRD=( [ctrl1]=32827137 [ctrl2]=30467885 [ctrl3]=25742013 [ctrl4]=30490491 [ctrl5]=34429360 [ctrl6]=32507852
  [MS1]=28737771 [MS2]=28588457 [MS3]=30792166 [MS4]=32868076 [MS5]=29789252 [MS6]=32036129 )
declare -A GRP=( [ctrl1]=control [ctrl2]=control [ctrl3]=control [ctrl4]=control [ctrl5]=control [ctrl6]=control
  [MS1]=MS [MS2]=MS [MS3]=MS [MS4]=MS [MS5]=MS [MS6]=MS )
ORDER="ctrl1 ctrl2 ctrl3 ctrl4 ctrl5 ctrl6 MS1 MS2 MS3 MS4 MS5 MS6"
ena(){ local s=$1; echo "https://ftp.sra.ebi.ac.uk/vol1/fastq/${s:0:6}/0${s: -2}/${s}"; }

# sample_map.csv for extract_lepr.py (generated; no separate upload needed)
echo "donor,group,SRR" > sample_map.csv
for nm in $ORDER; do echo "$nm,${GRP[$nm]},${RUN[$nm]}" >> sample_map.csv; done

# ---- Phase 1: parallel download 24 mates (6 concurrent, resume) ----
: > urls.txt
for nm in $ORDER; do s=${RUN[$nm]}; b=$(ena $s); echo "$b/${s}_1.fastq.gz"; echo "$b/${s}_2.fastq.gz"; done > urls.txt
echo ">>> downloading $(wc -l < urls.txt) files, 6 parallel"
xargs -P 6 -n 1 wget -c -q -P fastq/ < urls.txt || true
echo ">>> download phase done"; ls -la fastq/

# ---- Phase 2: salmon quant (paired; REUSE index) ----
for nm in $ORDER; do s=${RUN[$nm]}
  if [ -f quant/${nm}/quant.sf ]; then echo ">>> $nm already quantified"; continue; fi
  echo ">>> quant $nm ($s)"
  if [ -s fastq/${s}_1.fastq.gz ] && [ -s fastq/${s}_2.fastq.gz ]; then
    salmon quant -i "$IDX" -l A -1 fastq/${s}_1.fastq.gz -2 fastq/${s}_2.fastq.gz \
      -p $THREADS --validateMappings --gcBias -o quant/${nm}
  elif [ -s fastq/${s}_2.fastq.gz ]; then
    salmon quant -i "$IDX" -l A -r fastq/${s}_2.fastq.gz -p $THREADS --validateMappings -o quant/${nm}
  else echo "!!! $nm: no fastq"; fi
done

# ---- Phase 3: junction usage (R2; seqkit 30M seed42; single-pass awk count) ----
# probes (forward + reverse-complement); same as GSE228458 (validated vs Ensembl cDNA)
echo "donor,group,SRR,total_reads,seed,prop_sampled,n_reads_sampled,LepRa_reads,LepRb_reads,LepRa_CPM,LepRb_CPM,LepRb_fraction" > lepr_junction_GSE137619.csv
for nm in $ORDER; do
  s=${RUN[$nm]}; tot=${NRD[$nm]}; g=${GRP[$nm]}
  [ -s fastq/${s}_2.fastq.gz ] || { echo "$nm,$g,$s,$tot,42,NA,0,0,0,0,0,NA" >> lepr_junction_GSE137619.csv; continue; }
  p=$(awk -v t=$tot 'BEGIN{p=30000000/t; if(p>1)p=1; printf "%.6f",p}')
  read n a b < <(seqkit sample -s 42 -p $p fastq/${s}_2.fastq.gz 2>/dev/null | awk 'NR%4==2{
      n++;
      if (index($0,"ATTTTCAGAAGAGAACGGAC")||index($0,"GTCCGTTCTCTTCTGAAAAT")) a++;
      if (index($0,"ATTTTCAGAAGCCAGAAACG")||index($0,"CGTTTCTGGCTTCTGAAAAT")) bb++;
    } END{print n+0, a+0, bb+0}')
  acpm=$(awk -v x=$a -v n=$n 'BEGIN{printf "%.3f", (n>0)? x/n*1e6:0}')
  bcpm=$(awk -v x=$b -v n=$n 'BEGIN{printf "%.3f", (n>0)? x/n*1e6:0}')
  frac=$(awk -v x=$a -v y=$b 'BEGIN{printf "%.4f", (x+y>0)? y/(x+y):0}')
  echo "$nm,$g,$s,$tot,42,$p,$n,$a,$b,$acpm,$bcpm,$frac" >> lepr_junction_GSE137619.csv
  echo ">>> junction $nm: LepRa=$a LepRb=$b (sampled $n)"
done

# ---- Phase 4: expression extract (same script/IDs) ----
python3 extract_lepr.py || true
[ -f LEPR_expression_FINAL.csv ] && mv -f LEPR_expression_FINAL.csv LEPR_expression_GSE137619.csv
echo "DONE. -> LEPR_expression_GSE137619.csv + lepr_junction_GSE137619.csv + quant/*/quant.sf"
echo "(extract_lepr.py console summary lists control/AD/DM1; this dataset is control/MS — use the CSV; stats computed locally.)"
