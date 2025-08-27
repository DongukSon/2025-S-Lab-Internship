#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob    # 매치되는 파일이 없으면 for-loops 건너뜀

# rm -rf "$SUBJECTS_DIR/$SUBJ/stats/profiles/chi_p"
# rm -rf "$SUBJECTS_DIR/$SUBJ/stats/profiles/chi_n"
mkdir -p $SUBJECTS_DIR/$SUBJ/stats/profiles/chi_p
mkdir -p $SUBJECTS_DIR/$SUBJ/stats/profiles/chi_n
depths=(  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0 \
         -0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)

for h in lh rh; do
    for l in S_pc S_fs; do
        for d in "${depths[@]}"; do
            if [[ ! -f $SUBJECTS_DIR/$SUBJ/mri/roi/constrained_mask_ctx_${h}_${l}.mgz ]]; then
                echo "[skip] no mask"
                continue
            fi

            # if [[ -f $SUBJECTS_DIR/$SUBJ/stats/profiles/chi_p/ctx_${h}_${l}_d${d}.txt ]]; then
            #     echo "[skip] already depth profiled"
            #     continue
            # fi

            if (( $(echo "$d >= 0.0" | bc -l) )); then
                mri_vol2surf \
                --src       $SUBJECTS_DIR/$SUBJ/mri/roi/constrained_mask_ctx_${h}_${l}.mgz \
                --out       $SUBJECTS_DIR/$SUBJ/tmp/roi.mgh \
                --hemi      $h \
                --surf      equi${d}.pial \
                --interp    nearest \
                --out_type  mgh \
                --srcreg    $SUBJECTS_DIR/$SUBJ/register.dat
            else 
                mri_vol2surf \
                --src      $SUBJECTS_DIR/$SUBJ/mri/roi/constrained_mask_ctx_${h}_${l}.mgz \
                --out      $SUBJECTS_DIR/$SUBJ/tmp/roi.mgh \
                --hemi     $h \
                --surf     white \
                --projfrac $d \
                --interp   nearest \
                --out_type mgh \
                --srcreg   $SUBJECTS_DIR/$SUBJ/register.dat
            fi

            mri_cor2label \
            --i    $SUBJECTS_DIR/$SUBJ/tmp/roi.mgh \
            --id   1 \
            --l    $SUBJECTS_DIR/$SUBJ/tmp/roi.label \
            --surf $SUBJ $h white

            for suffix in p n; do
                mris_anatomical_stats \
                -l     $SUBJECTS_DIR/$SUBJ/tmp/roi.label \
                -t     $SUBJECTS_DIR/$SUBJ/surf/profiles/chi_${suffix}/${h}.d${d}.mgh \
                -f     $SUBJECTS_DIR/$SUBJ/stats/profiles/chi_${suffix}/ctx_${h}_${l}_d${d}.txt \
                $SUBJ  $h
            done
        done
    done
done




