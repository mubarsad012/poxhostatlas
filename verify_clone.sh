#!/usr/bin/env bash
# verify_clone.sh — simulate a judge cloning this repo from scratch.
#
# Checks, in order:
#   1. fresh clone from GitHub (nothing from your local disk)
#   2. every script run_all.sh invokes actually exists
#   3. all scripts compile
#   4. requirements.txt installs into a clean venv
#   5. the offline-runnable analyses regenerate their real figures/tables
#   6. regenerated numbers match the manuscript
#
# Usage:  bash verify_clone.sh [repo-url]
# Needs:  network (GitHub + PyPI). Steps 5-6 need no data download.

set -uo pipefail   # NOT -e: we want to report every failure, not stop at the first

REPO_URL="${1:-https://github.com/mubarsad012/poxhostatlas.git}"
WORK="$(mktemp -d)"; trap 'rm -rf "$WORK"' EXIT
PASS=0; FAIL=0
ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; PASS=$((PASS+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; FAIL=$((FAIL+1)); }
head_() { printf '\n=== %s ===\n' "$1"; }

head_ "1. fresh clone"
git clone --quiet "$REPO_URL" "$WORK/repo" && ok "cloned $REPO_URL" || { bad "clone failed"; exit 1; }
cd "$WORK/repo"
ok "HEAD $(git rev-parse --short HEAD)"

head_ "2. every script run_all.sh invokes exists"
grep -oE 'scripts/[0-9]+_[a-zA-Z0-9_]+\.py' run_all.sh | sort -u | while IFS= read -r s; do
  [ -f "$s" ] && printf '  PASS  %s\n' "$s" || printf '  FAIL  MISSING %s\n' "$s"
done
MISS=$(grep -oE 'scripts/[0-9]+_[a-zA-Z0-9_]+\.py' run_all.sh | sort -u | while IFS= read -r s; do
  [ -f "$s" ] || echo x; done | wc -l | tr -d ' ')
[ "$MISS" = "0" ] && ok "no missing scripts" || bad "$MISS script(s) referenced but absent"

head_ "3. all scripts compile"
python3 -m py_compile scripts/*.py 2>/dev/null && ok "py_compile scripts/*.py" || bad "syntax error in scripts/"

head_ "4. clean venv + requirements.txt"
python3 -m venv venv >/dev/null 2>&1
./venv/bin/pip install --quiet --upgrade pip >/dev/null 2>&1
if ./venv/bin/pip install --quiet -r requirements.txt >/dev/null 2>&1; then
  ok "requirements.txt installed"
else
  bad "pip install -r requirements.txt failed"
fi
./venv/bin/python - <<'PY' 2>/dev/null && ok "all imports resolve" || bad "an import failed"
import pandas, numpy, scipy, matplotlib, seaborn, sklearn, networkx, statsmodels, openpyxl, pydeseq2
PY
./venv/bin/python -c "import matplotlib,pandas,numpy;print('  versions: matplotlib',matplotlib.__version__,'pandas',pandas.__version__,'numpy',numpy.__version__)"

head_ "5. offline analyses regenerate their outputs"
# Delete the committed outputs FIRST, so the checks below prove real regeneration
# rather than just finding the copies that shipped with the clone.
rm -rf results/genomewide results/external_validation
# These three need no raw-data download: their inputs are committed.
for s in 22_genomewide_reproducibility 23_genomewide_network_ranking 24_matia_external_validation; do
  if MPLBACKEND=Agg ./venv/bin/python "scripts/$s.py" >"$WORK/$s.log" 2>&1; then
    ok "scripts/$s.py ran (exit 0)"
  else
    bad "scripts/$s.py FAILED — tail:"; tail -6 "$WORK/$s.log" | sed 's/^/        /'
  fi
done
for f in results/genomewide/figures/genomewide_volcano.png \
         results/genomewide/network/top_host_hubs.csv \
         results/genomewide/network/figures/genomewide_network.png \
         results/external_validation/figures/matia_singlecell_shutoff.png \
         results/external_validation/figures/matia_cross_modality.png; do
  [ -s "$f" ] && ok "regenerated $f" || bad "MISSING $f"
done

head_ "6. regenerated numbers match the manuscript"
./venv/bin/python - <<'PY'
import pandas as pd, sys
def chk(label, got, want, tol):
    good = got is not None and abs(got-want) <= tol
    print(f"  {'PASS' if good else 'FAIL'}  {label}: got {got}, manuscript {want}")
    return good
allgood = True
sc = pd.read_csv("results/external_validation/matia_singlecell_shutoff_summary.csv")
hf = sc[sc.metric.eq("human_frac")].set_index("cell_type")
allgood &= chk("HeLa host-read frac (uninf)",  round(float(hf.loc['HeLa','median_uninfected']),3), 0.998, 0.002)
allgood &= chk("HeLa host-read frac (inf)",    round(float(hf.loc['HeLa','median_infected']),2),   0.76,  0.02)
allgood &= chk("BJ5ta host-read frac (inf)",   round(float(hf.loc['BJ5ta','median_infected']),2),  0.80,  0.02)
g = pd.read_csv("results/genomewide/genomewide_reproducible_host_genes.csv")
allgood &= chk("reproducible host genes", len(g), 1834, 0)
c = pd.read_csv("results/genomewide/network/network_communities.csv")
allgood &= chk("largest co-regulated module", int(c.n_genes.max()), 128, 0)
sys.exit(0 if allgood else 1)
PY
[ $? -eq 0 ] && ok "manuscript numbers reproduce" || bad "a manuscript number did not reproduce"

head_ "SUMMARY"
printf '  passed: %s   failed: %s\n' "$PASS" "$FAIL"
if [ "$FAIL" -eq 0 ]; then
  printf '  \033[32mA judge cloning this repo can reproduce the offline analyses.\033[0m\n'
else
  printf '  \033[31mFix the FAIL lines above before submitting.\033[0m\n'; exit 1
fi
printf '\n  NOTE: scripts 01-21 additionally re-download raw counts from NCBI GEO;\n'
printf '        run `bash run_all.sh` on a networked machine for the full pipeline.\n'
