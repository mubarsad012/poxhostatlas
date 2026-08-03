#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Yang & Moss 2015 (SRP056975) ingestion: download FASTQ, trim, align to a
# combined human + Vaccinia-WR reference, and produce a HOST gene count matrix
# in the same format as data/processed/counts.csv so it can be added as a study
# to the cross-study meta-analysis (scripts/13_run_meta_analysis.py).
#
# Two read types are handled separately:
#   * RNA-seq (3 runs)            -> host transcript counts (feeds DGE/meta)
#   * Ribosome profiling / RPF    -> footprint counts (feeds translation
#     (12 runs, library=OTHER)       efficiency = RPF/RNA)
#
# Run on a machine WITH network + the conda env (see environment.yml).
# This is intentionally NOT part of run_all.sh (needs network, big refs, hours).
#
# REQUIRED before a mock-vs-infected contrast: fill sample_conditions.tsv (see
# README.md) — the SRA metadata has no sample titles, so conditions are unknown.
# ---------------------------------------------------------------------------
set -euo pipefail

# ---- knobs ----------------------------------------------------------------
THREADS="${THREADS:-8}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${WORK:-$HERE/work}"          # all intermediates + outputs land here
REF="${REF:-$WORK/ref}"
FASTQ="$WORK/fastq"; TRIM="$WORK/trim"; BAM="$WORK/bam"; CNT="$WORK/counts"
mkdir -p "$REF" "$FASTQ" "$TRIM" "$BAM" "$CNT"

# Ribosome-profiling 3' adapter — VERIFY against the SRA run metadata / paper.
# (classic Ingolia/Illumina small-RNA linker; override if the report differs)
RPF_ADAPTER="${RPF_ADAPTER:-CTGTAGGCACCATCAAT}"
RPF_MIN=26; RPF_MAX=34            # keep canonical footprint lengths

# featureCounts strandedness (0=unstranded, 1=forward, 2=reverse). DO NOT TRUST
# the defaults — the pipeline prints an empirical strand test (detect_strand) on
# the first RNA-seq BAM; set these to whichever maximizes 'Assigned'. 2015 HeLa
# RNA-seq is most often unstranded (0); Ingolia-style RPF is typically forward (1).
RNASEQ_STRAND="${RNASEQ_STRAND:-0}"
RPF_STRAND="${RPF_STRAND:-1}"
# One shared STAR index; 100 is a robust general overhang for ~50 bp RNA reads and
# short (26-34 nt) RPF reads alike (override per read length if desired).
SJDB_OVERHANG="${SJDB_OVERHANG:-100}"

# Optional contaminant (rRNA/tRNA) bowtie2 index prefix for RPF cleanup.
# Build once from human rRNA/tRNA FASTA: bowtie2-build contam.fa $REF/contam
CONTAM_INDEX="${CONTAM_INDEX:-$REF/contam}"

RNASEQ_RUNS=(SRR1959018 SRR1959019 SRR1959020)
RPF_RUNS=(SRR1959021 SRR1959022 SRR1959023 SRR1959024 SRR1959025 SRR1959026 \
          SRR1959027 SRR1959028 SRR1959029 SRR1959030 SRR1959031 SRR1959032)

# ---- references -----------------------------------------------------------
GENCODE=https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44
HUMAN_FA_URL=$GENCODE/GRCh38.primary_assembly.genome.fa.gz
GTF_URL=$GENCODE/gencode.v44.annotation.gtf.gz
# Vaccinia virus WR RefSeq genome (VERIFY accession): NC_006998.1
VACV_EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_006998.1&rettype=fasta&retmode=text"

build_reference() {
  echo "[ref] downloading + building combined human+VACV STAR index"
  [ -f "$REF/human.fa" ]  || { wget -nc "$HUMAN_FA_URL" -O "$REF/human.fa.gz"; gunzip -f "$REF/human.fa.gz"; }
  [ -f "$REF/anno.gtf" ]  || { wget -nc "$GTF_URL" -O "$REF/anno.gtf.gz"; gunzip -f "$REF/anno.gtf.gz"; }
  [ -f "$REF/vacv.fa" ]   || curl -s "$VACV_EFETCH" -o "$REF/vacv.fa"
  # guard: VACV-WR (NC_006998.1) is ~194,711 bp — fail loudly on an empty/wrong fetch
  local vbp; vbp=$(grep -v '^>' "$REF/vacv.fa" | tr -d '\n' | wc -c | tr -d ' ')
  if [ "$vbp" -lt 150000 ] || [ "$vbp" -gt 230000 ]; then
    echo "[ref] ERROR: VACV genome is $vbp bp, outside expected ~194,711 bp — check the NC_006998.1 fetch"; exit 1
  fi
  [ -f "$REF/genome.fa" ] || cat "$REF/human.fa" "$REF/vacv.fa" > "$REF/genome.fa"
  if [ ! -d "$REF/star" ]; then
    mkdir -p "$REF/star"
    STAR --runMode genomeGenerate --genomeDir "$REF/star" --genomeFastaFiles "$REF/genome.fa" \
         --sjdbGTFfile "$REF/anno.gtf" --sjdbOverhang "$SJDB_OVERHANG" --runThreadN "$THREADS" \
         --genomeSAindexNbases 14
  fi
}

get_fastq() { # $1 = SRR
  local r="$1"
  if [ ! -f "$FASTQ/$r.fastq.gz" ]; then
    echo "[dl] $r"
    # fasterq-dump is primary; ENA https fallback (00X subdir logic holds for these
    # 7-digit SRR accessions — verified against the vendored download script).
    fasterq-dump "$r" -e "$THREADS" -O "$FASTQ" --skip-technical 2>/dev/null && pigz -f "$FASTQ/$r.fastq" \
      || wget -nc "https://ftp.sra.ebi.ac.uk/vol1/fastq/${r:0:6}/00${r: -1}/$r/$r.fastq.gz" -O "$FASTQ/$r.fastq.gz"
  fi
}

align_count() { # $1=SRR $2=fastq $3=feature(exon|CDS) $4=strand(0|1|2) $5..=extra STAR args
  local r="$1" fq="$2" feat="$3" strand="$4"; shift 4
  STAR --runMode alignReads --genomeDir "$REF/star" --readFilesIn "$fq" --readFilesCommand zcat \
       --runThreadN "$THREADS" --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix "$BAM/${r}." "$@"
  featureCounts -T "$THREADS" -a "$REF/anno.gtf" -t "$feat" -g gene_id -s "$strand" \
       -o "$CNT/${r}.featureCounts.txt" "$BAM/${r}.Aligned.sortedByCoord.out.bam"
}

detect_strand() { # $1 = a BAM; prints 'Assigned' for -s 0/1/2 so you can set *_STRAND
  local bam="$1"
  echo "[strand] strandedness test on $(basename "$bam") — set RNASEQ_STRAND to the -s with most Assigned:"
  for s in 0 1 2; do
    featureCounts -T "$THREADS" -a "$REF/anno.gtf" -t exon -g gene_id -s "$s" \
      -o "$CNT/_strandtest_s${s}.txt" "$bam" >/dev/null 2>&1 || true
    awk -v s="$s" '/Assigned/{print "    -s "s"  Assigned="$2}' "$CNT/_strandtest_s${s}.txt.summary" 2>/dev/null || true
  done
}

process_rnaseq() {
  for r in "${RNASEQ_RUNS[@]}"; do
    get_fastq "$r"
    fastp -i "$FASTQ/$r.fastq.gz" -o "$TRIM/$r.trim.fastq.gz" -w "$THREADS" \
          -j "$TRIM/$r.fastp.json" -h "$TRIM/$r.fastp.html"
    align_count "$r" "$TRIM/$r.trim.fastq.gz" exon "$RNASEQ_STRAND"
  done
  detect_strand "$BAM/${RNASEQ_RUNS[0]}.Aligned.sortedByCoord.out.bam" || true
}

process_rpf() {
  for r in "${RPF_RUNS[@]}"; do
    get_fastq "$r"
    # 3' adapter trim + footprint length window
    cutadapt -a "$RPF_ADAPTER" -m "$RPF_MIN" -M "$RPF_MAX" --discard-untrimmed -j "$THREADS" \
             -o "$TRIM/$r.cut.fastq.gz" "$FASTQ/$r.fastq.gz"
    # remove rRNA/tRNA if a contaminant index was built
    if [ -f "${CONTAM_INDEX}.1.bt2" ]; then
      bowtie2 -x "$CONTAM_INDEX" -U "$TRIM/$r.cut.fastq.gz" -p "$THREADS" \
              --un-gz "$TRIM/$r.clean.fastq.gz" -S /dev/null 2>"$TRIM/$r.contam.log"
    else
      echo "[rpf] WARNING: no contaminant index at ${CONTAM_INDEX}; skipping rRNA/tRNA depletion"
      cp "$TRIM/$r.cut.fastq.gz" "$TRIM/$r.clean.fastq.gz"
    fi
    # footprints: end-to-end, no soft-clip; count over CDS
    align_count "$r" "$TRIM/$r.clean.fastq.gz" CDS "$RPF_STRAND" --alignEndsType EndToEnd \
                --outFilterMismatchNmax 2 --seedSearchStartLmax 15 --outFilterMultimapNmax 1
  done
}

merge_counts() { # $1 = label (rnaseq|rpf), shift = runs
  local label="$1"; shift
  # pass the GTF so we can attach gene_symbol (the meta-analysis keys on symbol)
  python3 - "$label" "$CNT" "$REF/anno.gtf" "$@" <<'PY'
import sys, re, pandas as pd
label, cntdir, gtf = sys.argv[1], sys.argv[2], sys.argv[3]; runs = sys.argv[4:]
mat = None
for r in runs:
    d = pd.read_csv(f"{cntdir}/{r}.featureCounts.txt", sep="\t", comment="#")
    s = d.set_index("Geneid")[d.columns[-1]].rename(r)
    mat = s.to_frame() if mat is None else mat.join(s, how="outer")
# Geneid -> gene_symbol from the GENCODE GTF gene lines
sym = {}
with open(gtf) as fh:
    for line in fh:
        if line.startswith("#") or "\tgene\t" not in line:
            continue
        gid = re.search(r'gene_id "([^"]+)"', line)
        gnm = re.search(r'gene_name "([^"]+)"', line)
        if gid and gnm:
            sym[gid.group(1)] = gnm.group(1)
out = mat.reset_index().rename(columns={"Geneid": "gene_id", "index": "gene_id"})
out.insert(1, "gene_symbol", out["gene_id"].map(sym))
# match data/processed/counts.csv column order: gene_id, gene_symbol, <run columns>
out.to_csv(f"{cntdir}/counts_{label}.csv", index=False)
print(f"[merge] wrote counts_{label}.csv  shape={out.shape}  cols=gene_id,gene_symbol,{len(runs)} runs")
PY
}

main() {
  build_reference
  process_rnaseq
  process_rpf
  merge_counts rnaseq "${RNASEQ_RUNS[@]}"
  merge_counts rpf "${RPF_RUNS[@]}"
  echo "[done] host count matrices in $CNT/counts_rnaseq.csv and counts_rpf.csv"
  echo "       next: map conditions (sample_conditions.tsv), then run a PyDESeq2"
  echo "       contrast like scripts/03_run_pydeseq2.py and add to the meta-analysis."
}
main "$@"
