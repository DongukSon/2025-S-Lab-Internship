

#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob    # skip unmatched globs

# ── sanity check ────────────────────────────────────────────────────────────────
: "${SUBJECTS_DIR:?Error:  SUBJECTS_DIR is not set}"

# ── iterate over all subjects ───────────────────────────────────────────────────
for subj_path in "$SUBJECTS_DIR"/Subj*; do
  [[ -d "$subj_path" ]] || continue          # skip if not a directory

  export SUBJ=$(basename "$subj_path")       # e.g. Subj014
  export ID="${SUBJ#Subj}"                   # e.g. 014
  echo "===== Processing $SUBJ (ID=$ID) ====="

  if [[ -f "Stat_final2/Subj${ID}.xlsx" ]]; then
      echo "[skip] ${ID}"
      continue
  fi
  
  SCRIPTS=(
    #"/home/donguk/research/Code/initial_set.sh"
    #"/home/donguk/research/Code/generate_surface.sh"
    #"/home/donguk/research/Code/proj_frac.sh"
    "/home/donguk/research/Code/depth_profile.sh"
    "/home/donguk/research/Code/create_excel.sh"
  )

  for script in "${SCRIPTS[@]}"; do
    if [[ -f "$script" ]]; then
      echo "[info] Running $script ..."
      if bash "$script"; then
        echo "[info] done: $script"
      else
        echo "[warn] error in $script — continuing" >&2
      fi
    else
      echo "[skip] $script missing."
    fi
  done

  echo "----- Finished $SUBJ -----"
done

echo "[done] all subjects processed"


: << 'BLOCK'
set -euo pipefail
shopt -s nullglob
IDS=( "004" "005" "006" "007" "022" 
      "023" "039" "040" "041" "048" 
      "049" "053" "055" "057" "060"
      "061" "062" "068" "069" "079" )

SCRIPTS=(
  #"/home/donguk/research/Code/initial_set.sh"
  #"/home/donguk/research/Code/generate_surface.sh"
  #"/home/donguk/research/Code/proj_frac.sh"
  "/home/donguk/research/Code/depth_profile.sh"
  "/home/donguk/research/Code/create_excel.sh"
)

for ID in "${IDS[@]}"; do
  SUBJ="Subj${ID}"
  echo
  echo "===== Processing ${SUBJ} ====="
  export ID
  export SUBJ

  for script in "${SCRIPTS[@]}"; do
    if [ -f "$script" ]; then
      echo "[info] Running $script ..."
      if ! bash "$script"; then
        echo "[warn] $script error occurred for ${SUBJ}, skipping to next script."
        continue
      fi
    else
      echo "[skip] $script missing."
    fi
  done

  echo "[info] Finished all scripts for ${SUBJ}"
done

echo
echo "[done] Finished processing all IDs"
BLOCK