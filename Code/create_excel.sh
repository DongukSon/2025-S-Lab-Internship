#!/usr/bin/env bash
set -euo pipefail

# ── 환경변수 SUBJ 확인 ───────────────────────────────
if [ -z "${SUBJ:-}" ]; then
  echo "Error: SUBJ 환경변수가 설정되어 있지 않습니다." >&2
  exit 1
fi

mkdir -p Stat_final2
OUTFILE="Stat_final2/${SUBJ}.xlsx"

# Depth 인덱스 0‒20 ↔ 파일 depth code 매핑

depths=(  0.0  0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1.0 \
         -0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7 -0.8 -0.9 -1.0)

# ROI prefix 매핑
declare -A PREFIX

# PREFIX[lh_psp]="ctx_lh_S_psp"
# PREFIX[rh_psp]="ctx_rh_S_psp"
# PREFIX[lh_fm]="ctx_lh_S_fm"
# PREFIX[rh_fm]="ctx_rh_S_fm"
# keys=( lh_psp_n lh_psp_p rh_psp_n rh_psp_p lh_fm_n lh_fm_p rh_fm_n rh_fm_p )

PREFIX[lh_pc]="ctx_lh_S_pc"
PREFIX[rh_pc]="ctx_rh_S_pc"
PREFIX[lh_fs]="ctx_lh_S_fs"
PREFIX[rh_fs]="ctx_rh_S_fs"
keys=( lh_pc_n lh_pc_p rh_pc_n rh_pc_p lh_fs_n lh_fs_p rh_fs_n rh_fs_p )

# ── CSV 임시파일 생성 ────────────────────────────────
TMP_AVG=$(mktemp)
TMP_STD=$(mktemp)
( IFS=','; echo "index,${keys[*]}" ) > "$TMP_AVG"
( IFS=','; echo "index,${keys[*]}" ) > "$TMP_STD"

for i in $(seq 0 20); do
  d=${depths[$i]}
  row_avg=("$i")
  row_std=("$i")

  for key in "${keys[@]}"; do
    IFS='_' read hemi roi np <<< "$key"
    base="${hemi}_${roi}"
    prefix="${PREFIX[$base]}"
    chi_dir="chi_${np}"
    txtdir="${SUBJECTS_DIR}/${SUBJ}/stats/profiles/${chi_dir}"
    txtfile="${txtdir}/${prefix}_d${d}.txt"

    if [ -f "$txtfile" ]; then
      # 5번째 필드(ThickAvg), 6번째 필드(ThickStd) 추출
      avg=$(awk '!/^#/ {last=$5} END{print last}' "$txtfile")
      std=$(awk '!/^#/ {last=$6} END{print last}' "$txtfile")
    else
      avg=""
      std=""
    fi
    row_avg+=("$avg")
    row_std+=("$std")
  done

  ( IFS=','; echo "${row_avg[*]}" ) >> "$TMP_AVG"
  ( IFS=','; echo "${row_std[*]}" ) >> "$TMP_STD"
done

# ── Python으로 Excel 작성 ───────────────────────────
python3 <<PYCODE
import csv
import xlsxwriter

# CSV 읽기
avg_rows = list(csv.reader(open(r"$TMP_AVG", newline='')))
std_rows = list(csv.reader(open(r"$TMP_STD", newline='')))

# 워크북 생성
wb = xlsxwriter.Workbook(r"$OUTFILE")
ws_avg = wb.add_worksheet("Avg")
ws_std = wb.add_worksheet("Std")

# Avg 시트 작성
for r, row in enumerate(avg_rows):
    for c, cell in enumerate(row):
        if r == 0:
            ws_avg.write(r, c, cell)
        else:
            if cell == "":
                ws_avg.write_blank(r, c, None)
            else:
                try:
                    ws_avg.write_number(r, c, float(cell))
                except ValueError:
                    ws_avg.write(r, c, cell)

# Std 시트 작성
for r, row in enumerate(std_rows):
    for c, cell in enumerate(row):
        if r == 0:
            ws_std.write(r, c, cell)
        else:
            if cell == "":
                ws_std.write_blank(r, c, None)
            else:
                try:
                    ws_std.write_number(r, c, float(cell))
                except ValueError:
                    ws_std.write(r, c, cell)

wb.close()
PYCODE

rm "$TMP_AVG" "$TMP_STD"
echo "Excel saved → $OUTFILE"