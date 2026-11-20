#!/bin/bash
#SBATCH --job-name=BRCA_ST
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vnikoone@asu.edu
#SBATCH --output=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/05_Downstream_Analysis/02_80_patients/05_logs/04_Neoantigen_Expression/04_StringTie_%A_%a.log

set -euo pipefail

BASE=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/05_Downstream_Analysis/02_80_patients
RESULTS="$BASE/03_results/04_Neoantigen_Expression"
MANIFEST="$RESULTS/00_StringTie_80_Genomic_RNA_BAMs.tsv"
GTF=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/06_NextNEOpi/nextNEOpi/resources/references/hg38/annotation/gencode.v47.primary_assembly.annotation.gtf

export PATH="/scratch/vnikoone/999_envs/stringtie/bin:$PATH"

module load samtools-1.21-gcc-12.1.0

echo "===== JOB INFORMATION ====="
date
hostname
id
echo "Job: ${SLURM_JOB_ID:-NA}"
echo "Array task: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "StringTie: $(stringtie --version)"
echo "Samtools: $(samtools --version | head -1)"
echo "==========================="

TASK=${SLURM_ARRAY_TASK_ID:?Missing array task ID}

[[ -s "$MANIFEST" ]] || {
    echo "ERROR: Missing manifest: $MANIFEST"
    exit 1
}

[[ -s "$GTF" ]] || {
    echo "ERROR: Missing GTF: $GTF"
    exit 1
}

LINE=$(sed -n "$((TASK + 1))p" "$MANIFEST")
[[ -n "$LINE" ]] || {
    echo "ERROR: No manifest entry for task $TASK"
    exit 1
}

IFS=$'\t' read -r PID BAM <<< "$LINE"
BAM=${BAM%$'\r'}

echo "Patient: $PID"
echo "BAM: $BAM"
echo "GTF: $GTF"

[[ "$PID" =~ ^TCGA-[A-Za-z0-9]{2}-[A-Za-z0-9]{4}$ ]] || {
    echo "ERROR: Invalid patient ID: $PID"
    exit 1
}

echo "===== BAM ACCESS CHECK ====="

ls -ld "$(dirname "$BAM")" || true
ls -lh "$BAM" || true

for attempt in 1 2 3; do
    if [[ -s "$BAM" && -r "$BAM" ]]; then
        echo "BAM accessible on attempt $attempt"
        break
    fi

    echo "WARNING: BAM not accessible on attempt $attempt"
    sleep 15
done

if [[ ! -s "$BAM" || ! -r "$BAM" ]]; then
    echo "ERROR: BAM remains inaccessible"
    echo "Compute node: $(hostname)"
    stat "$BAM" || true
    exit 1
fi

echo "===== BAM INTEGRITY ====="
samtools quickcheck -v "$BAM"

SORT_ORDER=$(samtools view -H "$BAM" | awk -F'\t' '$1=="@HD" {for (i=2;i<=NF;i++) if ($i ~ /^SO:/) {sub(/^SO:/,"",$i); print $i}}' | head -1)

echo "BAM sort order: $SORT_ORDER"

[[ "$SORT_ORDER" == coordinate ]] || {
    echo "ERROR: BAM is not coordinate sorted"
    exit 1
}

DIR="$RESULTS/00_StringTie_Quantification/$PID"
mkdir -p "$DIR"

TRANSCRIPTS="$DIR/${PID}.transcripts.gtf"
GENES="$DIR/${PID}.gene_abundance.tsv"
DONE="$DIR/${PID}.complete"

if [[ -s "$DONE" && -s "$TRANSCRIPTS" && -s "$GENES" ]]; then
    echo "Already completed: $PID"
    exit 0
fi

TMPDIR_RUN=$(mktemp -d "$DIR/.tmp.XXXXXXXX")
trap 'rm -rf "$TMPDIR_RUN"' EXIT

echo "===== STRINGTIE START ====="
date

stringtie "$BAM" \
    -p "${SLURM_CPUS_PER_TASK:-4}" \
    -e \
    -G "$GTF" \
    -o "$TMPDIR_RUN/transcripts.gtf" \
    -A "$TMPDIR_RUN/genes.tsv"

echo "===== OUTPUT VALIDATION ====="

[[ -s "$TMPDIR_RUN/transcripts.gtf" ]] || {
    echo "ERROR: Missing transcript GTF"
    exit 1
}

[[ -s "$TMPDIR_RUN/genes.tsv" ]] || {
    echo "ERROR: Missing gene abundance TSV"
    exit 1
}

grep -q $'\ttranscript\t' "$TMPDIR_RUN/transcripts.gtf" || {
    echo "ERROR: No transcript features found"
    exit 1
}

head -1 "$TMPDIR_RUN/genes.tsv" | grep -q 'TPM' || {
    echo "ERROR: Gene abundance output lacks TPM column"
    exit 1
}

mv "$TMPDIR_RUN/transcripts.gtf" "$TRANSCRIPTS"
mv "$TMPDIR_RUN/genes.tsv" "$GENES"

printf 'patient_id\t%s\ninput_bam\t%s\nstringtie_version\t%s\n' \
    "$PID" "$BAM" "$(stringtie --version)" > "$DONE"

echo "===== COMPLETED SUCCESSFULLY ====="
echo "Patient: $PID"
echo "Transcript output: $TRANSCRIPTS"
echo "Gene output: $GENES"
date
