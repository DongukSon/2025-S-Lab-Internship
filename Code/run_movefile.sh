#!/bin/bash

SRC_ROOT="/mnt/e/IDV"
DST_ROOT="/mnt/e/research_final/Data"

for subj_dir in "${SRC_ROOT}"/S*; do
  [ -d "$subj_dir" ] || continue

  subj=$(basename "$subj_dir")      # e.g. "S045"
  num=${subj#S}                     # e.g. "045"

  # 10#$num 로 숫자 앞의 0을 무시하고 정수 비교
  if [[ ! -f $DST_ROOT/$num/gre_t1wSinc.mgz ]]; then
    dst="${DST_ROOT}/${num}"
    mkdir -p "$dst"

    for suffix in scn_chi_n scn_chi_p gre_t1wSinc; do
      src_file="${subj_dir}/${subj}_${suffix}.nii.gz"
      if [ -f "$src_file" ]; then
        out_suffix="${suffix#scn_}"
        mri_convert "$src_file" "$dst/${out_suffix}.mgz"
        echo "Copied ${src_file} → ${dst}/${out_suffix}.mgz"
      else
        echo "Warning: ${src_file} not found, skipped."
      fi
    done
  else
    echo "Skipping ${subj}"
  fi
done
