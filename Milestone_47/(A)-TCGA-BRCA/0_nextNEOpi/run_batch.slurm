#!/bin/bash
#SBATCH --job-name=BRCA_NEO_B1
#SBATCH --partition=public
#SBATCH --qos=public
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=vnikoone@asu.edu
#SBATCH --output=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/07_BRCA_NextNEOpi/07_full_80_cohort/batch_1/logs/controller_%j.out
#SBATCH --error=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/07_BRCA_NextNEOpi/07_full_80_cohort/batch_1/logs/controller_%j.err

set -euo pipefail

ROOT=/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics
INSTALL="$ROOT/06_NextNEOpi"
REPO="$INSTALL/nextNEOpi"
CONFIG="$INSTALL/01_scripts/03_SOL_nextflow.config"
LICENSE=/scratch/vnikoone/999_envs/licenses/MiXCR_License.txt

BATCH_DIR="/scratch/vnikoone/01_ARPA/02_TCGA/03_BRCA/02_immunogenomics/07_BRCA_NextNEOpi/07_full_80_cohort/batch_1"
INPUT="$BATCH_DIR/batch_input.csv"
WORK="$BATCH_DIR/work"
RESULTS="$BATCH_DIR/results"
TMP="$BATCH_DIR/tmp"
LAUNCH="$BATCH_DIR/launch"

export PATH="/scratch/vnikoone/999_envs/nextneopi/bin:$PATH"
export NXF_VER=22.10.8
export NXF_CLASSPATH="$HOME/.nextflow/custom-libs/snappy-0.4.jar"
export NXF_OPTS="-Dleveldb.snappy=iq80"
export NXF_SYNTAX_PARSER=v1
export NXF_SINGULARITY_CACHEDIR="$INSTALL/containers"

test -s "$INPUT" || {
    echo "ERROR: Missing batch input: $INPUT"
    exit 1
}

test -s "$LICENSE" || {
    echo "ERROR: MiXCR license missing"
    exit 1
}

mkdir -p "$WORK" "$RESULTS" "$TMP" "$LAUNCH"

# Independent Nextflow launch directory and session/cache per batch.
cd "$LAUNCH"

# Fresh execution by default; use RESUME=1 to resume this batch.
RESUME_ARGS=()
if [[ "${RESUME:-0}" == "1" ]]; then
    RESUME_ARGS=(-resume)
fi

echo "Batch: batch_1"
echo "Launch directory: $LAUNCH"
echo "Work directory: $WORK"
echo "Resume: ${RESUME:-0}"
echo "Start: $(date)"

nextflow \
    -c "$REPO/conf/params.config" \
    -c "$CONFIG" \
    run "$REPO/nextNEOpi.nf" \
    "${RESUME_ARGS[@]}" \
    -profile singularity,cluster \
    -work-dir "$WORK" \
    --batchFile "$INPUT" \
    --outputDir "$RESULTS" \
    --tmpDir "$TMP" \
    --MIXCR_lic "$LICENSE" \
    --accept_license \
    --WES true \
    --TCR true \
    --disable_OptiType false \
    --exomeCaptureKit sureSelectV6 \
    -with-report "$BATCH_DIR/nextflow_report.html" \
    -with-timeline "$BATCH_DIR/nextflow_timeline.html" \
    -with-trace "$BATCH_DIR/nextflow_trace.tsv"

echo "Finished: $(date)"
