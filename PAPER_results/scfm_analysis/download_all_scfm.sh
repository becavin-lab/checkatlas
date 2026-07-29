#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# download_all_scfm.sh — reproduce all scFM-processed single-cell atlases
#
# Runs every download script in sequence to produce:
#   9 single-TF atlases    (*_TF_scfm.h5ad, from Census LTS 2025-11-08)
#   4 multi-scFM atlases   (*_multiscfm.h5ad, 2-Census merge, 6 scFMs each)
#   1 pancreas_scib        (original benchmark, from figshare)
#
# Prerequisites: pip install cellxgene-census tiledbsoma scanpy
# No GPU needed. No model weights downloaded.
#
# Usage:
#   bash download_all_scfm.sh                     # default output dir
#   bash download_all_scfm.sh /custom/output/dir  # custom output dir
# ---------------------------------------------------------------------------
set -euo pipefail

# ── Config ─────────────────────────────────────────────────────────────────

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${1:-/data/analysis/data_ganguly}"
START_TIME=$(date +%s)
TMPDIR="${TMPDIR:-/tmp}"
export SCFM_DATA_DIR="$OUTPUT_DIR"

RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'; NC='\033[0m'

# ── Helpers ────────────────────────────────────────────────────────────────

elapsed() { local now d; now=$(date +%s); d=$((now - START_TIME)); printf "%02d:%02d" $((d/60)) $((d%60)); }

step() {
    printf "\n${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
    printf "${GREEN}[%s] %s${NC}\n" "$(elapsed)" "$1"
    printf "${GREEN}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}\n"
}

say()   { printf "${CYAN}[%s] %s${NC}\n" "$(elapsed)" "$1"; }
warn()  { printf "${YELLOW}[%s] WARN:  %s${NC}\n" "$(elapsed)" "$1"; }
die()   { printf "${RED}[%s] FATAL: %s${NC}\n" "$(elapsed)" "$1"; exit 1; }

check_pkg() {
    python3 -c "import $1" 2>/dev/null || die "Missing Python module '$1' — run: pip install $2"
}

run_py() {
    local script="$1" desc="$2"
    say "Running ${desc}..."
    python3 "${SCRIPT_DIR}/${script}" || die "${desc} failed (see above)"
    say "${desc} OK"
}

# ── Preflight ──────────────────────────────────────────────────────────────

step "Preflight checks"

check_pkg cellxgene_census cellxgene-census
check_pkg tiledbsoma      tiledbsoma
check_pkg scanpy          scanpy

mkdir -p "$OUTPUT_DIR"
say "Output directory: ${OUTPUT_DIR}"

for s in download_census_tf.py download_tf_eval_datasets.py \
         download_pancreas_myeloid.py download_multiscfm.py \
         verify_all.py verify_multiscfm.py; do
    [[ -f "${SCRIPT_DIR}/${s}" ]] || die "Missing script: ${SCRIPT_DIR}/${s}"
done

# ── Tier 1 — Single-TF atlases ─────────────────────────────────────────────

step "Tier 1a — tissue-matched single-TF atlases (blood, bone_marrow, liver, lung)"
say "Source: Census LTS 2025-11-08 — embeddings: tf-sapiens + scvi"
run_py download_census_tf.py "Tier 1a (4 single-TF atlases)"

step "Tier 1b — TF evaluation datasets (pbmc_10k, covid19, tabula_sapiens_2)"
say "Source: Census LTS 2025-11-08 — embeddings: tf-sapiens + tf-exemplar + scvi"
run_py download_tf_eval_datasets.py "Tier 1b (3 evaluation atlases)"

step "Tier 1c — pancreas + myeloid analogues"
say "Source: Census LTS 2025-11-08 — embeddings: tf-sapiens + tf-exemplar + scvi"
run_py download_pancreas_myeloid.py "Tier 1c (hpancreas + myeloid)"

# ── Tier 2 — Multi-scFM atlases ────────────────────────────────────────────

step "Tier 2 — multi-scFM atlases (6 SOTA scFMs, gold-standard join)"
say "Source: Census 2025-01-30 (TF, Geneformer, scVI) + Census 2023-12-15 (scGPT, UCE)"
say "Join: row-position with 3-way verification (obs_names + donor_id + assay = 100%)"
run_py download_multiscfm.py "Tier 2 (4 multi-scFM atlases)"

# ── Verification ───────────────────────────────────────────────────────────

step "Verification"

run_py verify_all.py "verify single-TF atlases (9 files)"
run_py verify_multiscfm.py "verify multi-scFM atlases (4 files)"

# ── Tier 3 — pancreas_scib ─────────────────────────────────────────────────

step "Tier 3 — pancreas_scib (original scIB benchmark)"

PANCREAS_PATH="${OUTPUT_DIR}/pancreas_scib.h5ad"
if [[ -f "$PANCREAS_PATH" ]]; then
    size=$(stat -c%s "$PANCREAS_PATH" 2>/dev/null || echo 0)
    if [[ $size -gt 5000000 ]]; then
        say "pancreas_scib.h5ad already exists ($(du -h "$PANCREAS_PATH" | cut -f1), skipping)"
    fi
else
    PANCREAS_URL="https://figshare.com/ndownloader/files/24539828"
    say "Downloading from figshare ($PANCREAS_URL)..."
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "$PANCREAS_PATH" "$PANCREAS_URL" || die "wget failed"
    elif command -v curl &>/dev/null; then
        curl -L -o "$PANCREAS_PATH" "$PANCREAS_URL" || die "curl failed"
    else
        warn "No wget/curl — download manually: ${PANCREAS_URL}"
    fi
    [[ -f "$PANCREAS_PATH" ]] && say "Downloaded ($(du -h "$PANCREAS_PATH" | cut -f1))"
fi

# ── Summary ────────────────────────────────────────────────────────────────

step "Done — all scFM atlases produced"

printf "\n${GREEN}Files in %s:${NC}\n" "$OUTPUT_DIR"
ls -lh "${OUTPUT_DIR}"/*_TF_scfm.h5ad \
      "${OUTPUT_DIR}"/*_multiscfm.h5ad \
      "${OUTPUT_DIR}"/pancreas_scib.h5ad 2>/dev/null \
    | awk '{ printf "  %-45s %s\n", $9, $5 }' || true

printf "\n${GREEN}Wall-clock: %s${NC}\n" "$(elapsed)"
printf "${GREEN}Output:    %s${NC}\n" "$OUTPUT_DIR"
