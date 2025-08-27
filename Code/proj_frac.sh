
depths=(  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0 \
         -0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)
hemis=( lh rh )
mkdir -p $SUBJECTS_DIR/$SUBJ/surf/profiles
mkdir -p $SUBJECTS_DIR/$SUBJ/surf/profiles/chi_p
mkdir -p $SUBJECTS_DIR/$SUBJ/surf/profiles/chi_n

for h in "${hemis[@]}"; do
    for suffix in p n; do
        for d in "${depths[@]}"; do
            if (( $(echo "$d >= 0.0" | bc -l) )); then
                mri_vol2surf \
                --src       $SUBJECTS_DIR/$SUBJ/mri/chi_${suffix}.mgz \
                --out       $SUBJECTS_DIR/$SUBJ/surf/profiles/chi_${suffix}/${h}.d${d}.mgh \
                --hemi      $h \
                --surf      equi${d}.pial \
                --interp    trilinear \
                --out_type  mgh \
                --srcreg    $SUBJECTS_DIR/$SUBJ/register.dat

            else
                mri_vol2surf \
                --src       $SUBJECTS_DIR/$SUBJ/mri/chi_${suffix}.mgz \
                --out       $SUBJECTS_DIR/$SUBJ/surf/profiles/chi_${suffix}/${h}.d${d}.mgh \
                --hemi      $h \
                --projfrac  $d \
                --surf      white \
                --interp    trilinear \
                --out_type  mgh \
                --srcreg    $SUBJECTS_DIR/$SUBJ/register.dat
            fi
        done
    done
done