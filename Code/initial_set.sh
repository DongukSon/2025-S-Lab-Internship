#!/usr/bin/env bash
echo $SUBJ
  # mri_convert \
  #   $SUBJECTS_DIR/$SUBJ/mri/brain.mgz \
  #   $SUBJECTS_DIR/$SUBJ/mri/brain.nii.gz

  # flirt \
  #   -in    /mnt/e/IDV/S${ID}/S${ID}_scn_chi_n.nii.gz \
  #   -ref   $SUBJECTS_DIR/$SUBJ/mri/brain.nii.gz \
  #   -dof   6 \
  #   -cost  mutualinfo \
  #   -omat  $SUBJECTS_DIR/$SUBJ/chi2T1.mat \
  #   -interp sinc \
  #   -out   $SUBJECTS_DIR/$SUBJ/mri/chi_n.nii.gz
  # flirt \
  #   -in  /mnt/e/IDV/S${ID}/S${ID}_scn_chi_p.nii.gz \
  #   -ref $SUBJECTS_DIR/$SUBJ/mri/brain.nii.gz\
  #   -applyxfm \
  #   -interp sinc \
  #   -init $SUBJECTS_DIR/$SUBJ/chi2T1.mat \
  #   -out $SUBJECTS_DIR/$SUBJ/mri/chi_p.nii.gz
  # mri_convert $SUBJECTS_DIR/$SUBJ/mri/chi_p.nii.gz $SUBJECTS_DIR/$SUBJ/mri/chi_p2.mgz
  # mri_convert $SUBJECTS_DIR/$SUBJ/mri/chi_n.nii.gz $SUBJECTS_DIR/$SUBJ/mri/chi_n2.mgz

# if [ ! -f "$SUBJECTS_DIR/$SUBJ/mri/chi_p.mgz" ]; then
  mri_vol2vol \
    --mov       /mnt/e/research_final/Data/${ID}/chi_p.mgz \
    --targ      $SUBJECTS_DIR/$SUBJ/mri/T1.mgz \
    --regheader \
    --o         $SUBJECTS_DIR/$SUBJ/mri/chi_p.mgz \
    --interp    cubic
  mri_vol2vol \
    --mov       /mnt/e/research_final/Data/${ID}/chi_n.mgz \
    --targ      $SUBJECTS_DIR/$SUBJ/mri/T1.mgz \
    --regheader \
    --o         $SUBJECTS_DIR/$SUBJ/mri/chi_n.mgz \
    --interp    cubic
# fi

if [ ! -f "$SUBJECTS_DIR/$SUBJ/label/lh.S_front_middle.label" ]; then
  mri_annotation2label \
    --subject    $SUBJ \
    --hemi       lh \
    --annotation aparc.a2009s \
    --outdir     $SUBJECTS_DIR/$SUBJ/label
  mri_annotation2label \
    --subject    $SUBJ \
    --hemi       rh \
    --annotation aparc.a2009s \
    --outdir     $SUBJECTS_DIR/$SUBJ/label
else
  echo "[skip] labels are already generated ..."
fi

# if [ ! -f  "$SUBJECTS_DIR/$SUBJ/register.dat" ]; then
  tkregister2 \
    --mov  $SUBJECTS_DIR/$SUBJ/mri/chi_n.mgz \
    --targ $SUBJECTS_DIR/$SUBJ/mri/orig.mgz \
    --s    $SUBJ \
    --reg  $SUBJECTS_DIR/$SUBJ/register.dat \
    --noedit \
    --regheader
# else
#   echo "[skip] register.dat is already generated ..."
# fi