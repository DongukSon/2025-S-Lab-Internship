#!/usr/bin/env bash
set -euo pipefail
export ID=035
# 0) qform 비어 있으면 sform→qform 복사(무해)
# fslorient -copysform2qform "/mnt/e/IDV/S${ID}/S${ID}_gre_mag_e1_priored_corrected.nii.gz"
# fslorient -copysform2qform "/mnt/e/IDV/S${ID}/S${ID}_mpr_t1w.nii.gz"
# fslorient -copysform2qform "/mnt/e/IDV/S${ID}/S${ID}_scn_chi_n.nii.gz"
# fslorient -copysform2qform "/mnt/e/IDV/S${ID}/S${ID}_scn_chi_p.nii.gz"

# 1) 모두 fslreorient2std (정수 회전/flip만; 보간 없음)
fslreorient2std "/mnt/e/IDV/S${ID}/S${ID}_gre_mag_e1_priored_corrected.nii.gz" "/mnt/e/research/Data/${ID}/gre_std.nii.gz"
fslreorient2std "/mnt/e/IDV/S${ID}/S${ID}_mpr_t1w.nii.gz" "/mnt/e/research/Data/${ID}/t1_std.nii.gz"
fslreorient2std "/mnt/e/IDV/S${ID}/S${ID}_scn_chi_n.nii.gz" "/mnt/e/research/Data/${ID}/chi_n_std.nii.gz"
fslreorient2std "/mnt/e/IDV/S${ID}/S${ID}_scn_chi_p.nii.gz" "/mnt/e/research/Data/${ID}/chi_p_std.nii.gz"

# 2) GRE → T1 6-DOF 정합 (표준 방향끼리)
flirt \
  -in  "/mnt/e/research/Data/${ID}/gre_std.nii.gz" \
  -ref "/mnt/e/research/Data/${ID}/t1_std.nii.gz" \
  -dof 6 -cost normmi -bins 256 \
  -omat "/mnt/e/research/Data/${ID}/gre2t1_std.mat" \
  -out "/mnt/e/research/Data/${ID}/gre_on_t1_std_sinc.nii.gz" \
  -interp sinc

# 3) χ를 T1(std)에 직접 적용 (아이덴티티 리샘플 생략)
flirt \
  -in  "/mnt/e/research/Data/${ID}/chi_n_std.nii.gz" \
  -ref "/mnt/e/research/Data/${ID}/t1_std.nii.gz" \
  -applyxfm -init "/mnt/e/research/Data/${ID}/gre2t1_std.mat" \
  -out "/mnt/e/research/Data/${ID}/chi_n_on_t1_std.nii.gz" \
  -interp sinc

flirt \
  -in  "/mnt/e/research/Data/${ID}/chi_p_std.nii.gz" \
  -ref "/mnt/e/research/Data/${ID}/t1_std.nii.gz" \
  -applyxfm -init "/mnt/e/research/Data/${ID}/gre2t1_std.mat" \
  -out "/mnt/e/research/Data/${ID}/chi_p_on_t1_std.nii.gz" \
  -interp sinc

# 4) MGZ 변환
mri_convert "/mnt/e/research/Data/${ID}/chi_n_on_t1_std.nii.gz" "/mnt/e/research/Data/${ID}/chi_n_.mgz"
mri_convert "/mnt/e/research/Data/${ID}/chi_p_on_t1_std.nii.gz" "/mnt/e/research/Data/${ID}/chi_p_.mgz"


# export ID=036
# flirt \
#     -in  /mnt/e/IDV/S${ID}/S${ID}_gre_mag_e1_priored_corrected.nii.gz \
#     -ref /mnt/e/IDV/S${ID}/S${ID}_mpr_t1w.nii.gz \
#     -dof  6 -cost normmi -bins 256 \
#     -omat /mnt/e/research/Data/$ID/gre2t1.mat \
#     -out  /mnt/e/research/Data/$ID/gre_on_t1_sinc.nii.gz \
#     -interp sinc
    
# flirt \
#     -in /mnt/e/IDV/S${ID}/S${ID}_scn_chi_n.nii.gz \
#     -ref /mnt/e/IDV/S${ID}/S${ID}_mpr_t1w.nii.gz \
#     -applyxfm \
#     -init /mnt/e/research/Data/$ID/gre2t1.mat \
#     -out /mnt/e/research/Data/$ID/chi_n.nii.gz \
#     -interp sinc

# flirt \
#     -in /mnt/e/IDV/S${ID}/S${ID}_scn_chi_p.nii.gz \
#     -ref /mnt/e/IDV/S${ID}/S${ID}_mpr_t1w.nii.gz \
#     -applyxfm \
#     -init /mnt/e/research/Data/$ID/gre2t1.mat \
#     -out /mnt/e/research/Data/$ID/chi_p.nii.gz \
#     -interp sinc

# mri_convert /mnt/e/research/Data/$ID/chi_n.nii.gz /mnt/e/research/Data/$ID/chi_n_.mgz
# mri_convert /mnt/e/research/Data/$ID/chi_p.nii.gz /mnt/e/research/Data/$ID/chi_p_.mgz