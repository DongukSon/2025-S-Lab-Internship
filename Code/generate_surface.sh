if [ ! -f "$SUBJECTS_DIR/$SUBJ/surf/rh.equi1.0.pial" ]; then
    cd /home/donguk/surface_tools
    hemis=( lh rh )
    for h in "${hemis[@]}"; do
        python -m surface_tools.equivolumetric_surfaces.generate_equivolumetric_surfaces \
        --smoothing 0 \
        $SUBJECTS_DIR/$SUBJ/surf/${h}.pial \
        $SUBJECTS_DIR/$SUBJ/surf/${h}.white \
        11 \
        $SUBJECTS_DIR/$SUBJ/surf/${h}.equi \
        --software freesurfer \
        --subject_id $SUBJ
    done
else
    echo "[info] surfaces are already generated ..."
fi
