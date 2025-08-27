#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import nibabel as nib
import subprocess
from nibabel.freesurfer import read_geometry, read_annot
from nibabel.freesurfer.mghformat import MGHImage
from scipy.ndimage import distance_transform_edt
from scipy.spatial import cKDTree
from pathlib import Path
import os
from scipy.ndimage import binary_closing, binary_fill_holes
from nibabel.affines import apply_affine
from nibabel.freesurfer.io import read_morph_data

TXT_PATH     = '/home/donguk/research/ROI_tools/NewSliceNum.txt'
SUBJECTS_DIR = Path(os.environ['SUBJECTS_DIR'])
HEMI         = ['lh', 'rh']
ROI_TYPE     = 2
# TARGET_LABEL_IDS = {'S_psp': 70, 'S_fm' : 54}
# seed_radius   = {'S_psp': 1.0, 'S_fm': 1.0}
# expand_radius = {'S_psp': 6.0, 'S_fm': 6.0}
TARGET_LABEL_IDS = {'S_pc': 68, 'S_fs' : 55}
seed_radius   = {'S_pc': 1.0, 'S_fs': 1.0}
expand_radius = {'S_pc': 6.0, 'S_fs': 6.0}

def build_roi_cfg(slice_vals, type):
    if type == 1:
        return {
            'lh': {
                'S_psp': dict(slice_idx=slice_vals[0], label_id=11170),
                'S_fm' : dict(slice_idx=slice_vals[2], label_id=11154),
            },
            'rh': {
                'S_psp': dict(slice_idx=slice_vals[1], label_id=12170),
                'S_fm' : dict(slice_idx=slice_vals[3], label_id=12154),
            },
        }
    elif type == 2:
        return {
            'lh': {
                'S_pc': dict(slice_idx=slice_vals[0], label_id=11168),
                'S_fs' : dict(slice_idx=slice_vals[2], label_id=11155),
            },
            'rh': {
                'S_pc': dict(slice_idx=slice_vals[1], label_id=12168),
                'S_fs' : dict(slice_idx=slice_vals[3], label_id=12155),
            },
        }

def _to_int(tok: str) -> int:
    """'N/A', '', '-', '0' → 0, 그 외 정수 문자열은 int로 변환."""
    tok = tok.strip()
    if tok in {'', '-', 'N/A', 'NA', '0'}:
        return 0
    return int(tok)

def load_slice_table(path):
    rows = []
    with open(path, 'r', encoding='utf-8') as f:
        for ln in f:
            ln = ln.strip()
            if not ln:
                continue
            tokens = ln.split()
            subj  = tokens[0]
            idxs  = list(map(_to_int, tokens[1:5]))
            rows.append((subj, idxs))
    return rows

def make_coronal_mask(t1_img, slice_idx, out_path):
    m = np.zeros(t1_img.shape, dtype=np.uint8)
    m[:, :, slice_idx] = 1
    nib.save(nib.MGHImage(m, t1_img.affine, t1_img.header), str(out_path))

def make_binarize_mask(mask_path, label_id, out_path, aseg_path):
    subprocess.check_call([
        'mri_binarize', 
        '--i', str(aseg_path), 
        '--o', str(out_path),
        '--match', str(label_id), 
        '--mask', str(mask_path),
    ])

def make_within_mask(bin_path, out_path, r_seed, hemi):
    img       = nib.load(str(bin_path))
    data_bool = img.get_fdata().astype(bool)
    vox_idx   = np.column_stack(np.nonzero(data_bool))
    world_pts = nib.affines.apply_affine(img.affine, vox_idx)
    #   lh: x - z,  rh: -x - z
    if hemi == 'lh':
        scores = world_pts[:, 0] - 1.8*world_pts[:, 2]
    else:  # 'rh'
        scores = -world_pts[:, 0] - 1.8*world_pts[:, 2]
        
    max_score  = scores.max()
    candidates = np.where(scores == max_score)[0]
    idx_seed   = candidates[np.argmin(world_pts[candidates, 2])]
    seed_world = world_pts[idx_seed]
    within = np.linalg.norm(world_pts - seed_world, axis=1) <= r_seed
    mask = np.zeros_like(data_bool, dtype=np.uint8)
    mask[tuple(vox_idx[within].T)] = 1
    nib.save(nib.MGHImage(mask, img.affine, img.header), str(out_path))

def make_expanded_mask(within_path, out_path, r_expand, aseg_path, label_id):
    img       = nib.load(str(within_path))
    mask3d    = img.get_fdata().astype(bool)
    affine    = img.affine
    header    = img.header

    zooms     = header.get_zooms()[:3]             # (dx, dy, dz) in mm
    dt3d      = distance_transform_edt(~mask3d, sampling=zooms)

    close_enough = (dt3d <= r_expand)

    # load aseg volume to check label IDs
    aseg_img  = nib.load(str(aseg_path))
    aseg_data = aseg_img.get_fdata().astype(int)

    # build final mask: within distance AND matching label
    final_mask = np.zeros_like(mask3d, dtype=np.uint8)
    sel = close_enough & (aseg_data == label_id)
    final_mask[sel] = 1
    nib.save(MGHImage(final_mask, affine, header), str(out_path))

def make_depth_mask(expanded_path, out_path, surf_path, aseg_path):
    img = nib.load(str(expanded_path))
    affine      = img.affine
    header      = img.header
    vox2ras_tkr = img.header.get_vox2ras_tkr()
    inv_tkr     = np.linalg.inv(vox2ras_tkr)

    # aseg 파일을 로드해서 해부학적 레이블 정보를 가져옵니다.
    aseg_img = nib.load(str(aseg_path))
    aseg_data = aseg_img.get_fdata().astype(int)
    # --- 새로운 부분 끝 ---

    verts, faces = read_geometry(str(surf_path))
    # face normal
    face_normals = np.cross(
        verts[faces[:, 1]] - verts[faces[:, 0]],
        verts[faces[:, 2]] - verts[faces[:, 0]]
    )
    face_normals /= np.linalg.norm(face_normals, axis=1, keepdims=True)
    # vertex normal
    vertex_normals = np.zeros_like(verts)
    for i in range(faces.shape[0]):
        for j in range(3):
            vertex_normals[faces[i, j]] += face_normals[i]
    vertex_normals /= np.linalg.norm(vertex_normals, axis=1, keepdims=True)

    vert_ijk = np.round(nib.affines.apply_affine(inv_tkr, verts)).astype(int)

    exp_img  = nib.load(str(expanded_path))
    exp_data = exp_img.get_fdata().astype(bool)

    valid = (
        (vert_ijk[:,0]>=0)&(vert_ijk[:,0]<exp_data.shape[0]) &
        (vert_ijk[:,1]>=0)&(vert_ijk[:,1]<exp_data.shape[1]) &
        (vert_ijk[:,2]>=0)&(vert_ijk[:,2]<exp_data.shape[2])
    )

    inside = np.zeros(verts.shape[0], dtype=bool)
    inside[valid] = exp_data[
        vert_ijk[valid, 0],
        vert_ijk[valid, 1],
        vert_ijk[valid, 2]
    ]
    keep = np.where(inside)[0]

    step_mm = 0.5
    n_steps = 12
    mask = np.zeros(img.shape, dtype=np.uint8)

    neighbor_offsets = [
        (0,  0,  0), (0, 0, 1), (0, 0, -1)
    ]

    for d in range(0, n_steps + 1):
        ras_coords = verts[keep] - d * step_mm * vertex_normals[keep]
        ijk_coords = nib.affines.apply_affine(inv_tkr, ras_coords)
        ijk_coords = np.round(ijk_coords).astype(int)

        for i, j, k in ijk_coords:
            for di, dj, dk in neighbor_offsets:
                ii, jj, kk = i + di, j + dj, k + dk
                if (0 <= ii < mask.shape[0] and 0 <= jj < mask.shape[1] and 0 <= kk < mask.shape[2]):
                    # 해당 복셀이 화이트 매터(레이블 2 또는 41)일 경우에만 마스크에 기록합니다.
                    if aseg_data[ii, jj, kk] == 2 or aseg_data[ii, jj, kk] == 41:
                        mask[ii, jj, kk] = 1
                    # --- 수정된 부분 끝 ---

    nib.save(nib.MGHImage(mask, affine, header), str(out_path))
    print(f"out_img saved to {out_path}")

def make_constrained_mask(expanded_path, depth_path, slice_idx, out_path, close_radius):
    img_exp   = nib.load(str(expanded_path))
    data_exp  = img_exp.get_fdata().astype(bool)
    img_norm  = nib.load(str(depth_path))
    data_norm = img_norm.get_fdata().astype(bool)

    combined = data_exp | data_norm

    mask2d = combined[:, :, slice_idx]

    # define a square structuring element of size (2*close_radius+1)
    struct = np.ones((2*close_radius+1, 2*close_radius+1), dtype=bool)

    # 1) close small gaps
    closed = binary_closing(mask2d, structure=struct)
    # 2) fill interior holes
    filled = binary_fill_holes(closed)

    # build full-volume output but only that slice is non-zero
    out_mask = np.zeros_like(combined, dtype=np.uint8)
    out_mask[:, :, slice_idx] = filled.astype(np.uint8)

    # save as MGZ
    nib.save(
        MGHImage(out_mask, img_exp.affine, img_exp.header),
        str(out_path)
    )
    print(f"constrained mask saved to {out_path}")

def load_surface(subj_dir, subj, hemi):
    surf  = subj_dir / subj / "surf"  / f"{hemi}.white"
    annot = subj_dir / subj / "label" / f"{hemi}.aparc.a2009s.annot"
    verts, _         = read_geometry(surf)
    labels, _, names = read_annot(annot)
    names = [n.decode("utf-8") for n in names]
    return verts, labels, names

def get_thickness_thresholds(subj_path, hemi, vlabels, label_ids, scale=1.5):
    t_path = subj_path / 'surf' / f'{hemi}.thickness'
    thickness = read_morph_data(str(t_path))     # (V,) 단위: mm

    thresholds = {}
    for lb in label_ids:
        t_vals = thickness[vlabels == lb]
        if len(t_vals) == 0:
            # fallback: 3.5 mm (또는 다른 값)
            thresholds[lb] = 3.5
        else:
            thresholds[lb] = float(np.mean(t_vals) * scale)
    return thresholds

def get_n_ring_neighbors(vertex_adj, vid, n):
    """
    vertex_adj: list of sets, vertex_adj[v]는 v의 1-ring 이웃 set
    vid: 시작 정점 인덱스
    n: 확장할 ring 수 (1이면 1-ring, 2면 2-ring)
    """
    visited = {vid}
    frontier = {vid}
    for _ in range(n):
        # frontier에 있는 모든 정점의 1-ring 이웃을 모아서
        next_frontier = set().union(*(vertex_adj[v] for v in frontier))
        # 이미 방문한 정점 빼고
        frontier = next_frontier - visited
        visited |= frontier
    visited.remove(vid)
    return visited

def run_pipeline(subj, slice_vals):
    os.environ['SUBJ'] = subj
    roi_cfg   = build_roi_cfg(slice_vals, ROI_TYPE)
    subj_path = SUBJECTS_DIR / subj
    t1_path   = subj_path / 'mri' / 'T1.mgz'
    aseg_path = subj_path / 'mri' / 'aparc.a2009s+aseg.mgz'
    tmp_dir   = subj_path / 'tmp'
    roi_dir   = subj_path / 'mri' / 'roi'
    # if roi_dir.exists():
    #     shutil.rmtree(roi_dir)
    # if tmp_dir.exists():
    #     shutil.rmtree(tmp_dir)
    roi_dir.mkdir(exist_ok=True)
    tmp_dir.mkdir(exist_ok=True)

    t1_img   = nib.load(str(t1_path))
    aseg_img = nib.load(str(aseg_path))    

    for hemi, roi_dict in roi_cfg.items():
        for roi, cfg in roi_dict.items():
            slice_idx = cfg['slice_idx']
            if slice_idx == 0:
                print(f'[skip] {subj}:{hemi}_{roi} (slice_idx={slice_idx})')
                continue
            label_id  = cfg['label_id']
            r_seed    = seed_radius[roi]
            r_expand  = expand_radius[roi]
            prefix = f'ctx_{hemi}_{roi}'
            surf_path = subj_path / 'surf' / f'{hemi}.white'
            cor = tmp_dir / f'coronal_mask_{prefix}.mgz'
            bin = tmp_dir / f'binarized_mask_{prefix}.mgz'
            wit = tmp_dir / f'within_mask_{prefix}.mgz'
            exp = tmp_dir / f'expanded_mask_{prefix}.mgz'
            dep = tmp_dir / f'depth_mask_{prefix}.mgz'
            con = tmp_dir / f'constrained_mask_{prefix}.mgz'

            print(f'\n── [{subj}] {prefix}: slice {slice_idx}, label {label_id} ──')
            make_coronal_mask(t1_img, slice_idx, cor)
            make_binarize_mask(cor, label_id, bin, aseg_path)
            make_within_mask(bin, wit, r_seed, hemi)
            make_expanded_mask(wit, exp, r_expand, aseg_path, label_id)
            make_depth_mask(exp, dep, surf_path, aseg_path)
            make_constrained_mask(exp, dep, slice_idx, con, 3)
    vox2ras_tkr = aseg_img.header.get_vox2ras_tkr()
    for hemi in HEMI:
        surface = subj_path/'surf'/f'{hemi}.white'
        verts, faces = read_geometry(str(surface))
        vlabels, _, _ = read_annot(str(subj_path/'label'/f'{hemi}.aparc.a2009s.annot'))
        vertex_adj = [set() for _ in range(len(verts))]
        for a, b, c in faces:
            vertex_adj[a].update((b, c))
            vertex_adj[b].update((a, c))
            vertex_adj[c].update((a, b))
        n_verts = verts.shape[0]
        for target_label_name, target_label_id in TARGET_LABEL_IDS.items():
            mask_path = tmp_dir / f"constrained_mask_ctx_{hemi}_{target_label_name}.mgz"
            if not mask_path.exists():
                print(f"[skip] {mask_path} not found — skipped")
                continue
            t_path = subj_path/'surf'/f'{hemi}.thickness'
            thickness = read_morph_data(str(t_path))  # shape (n_verts,), 단위 mm
            # (1) mask voxels
            mask_img  = nib.load(mask_path)
            mask_data = mask_img.get_fdata()
            ijk       = np.column_stack(np.nonzero(mask_data == 1)).astype(int)  # (N,3)

            # (2) tkRAS coords
            world = apply_affine(vox2ras_tkr, ijk)                         # (N,3)

            # ── (3) slice‑based KD‑tree ───────────────────────────────────────────────────
            #  (A) coronal mask 불러오기 (make_coronal_mask 때 생성된 파일)
            cor_path = tmp_dir / f"coronal_mask_ctx_{hemi}_{target_label_name}.mgz"
            cor_data = nib.load(str(cor_path)).get_fdata().astype(bool)
            
            #  (B) 모든 표면 정점 → voxel index 변환
            ras2vox   = np.linalg.inv(vox2ras_tkr)
            vox_verts = np.round(apply_affine(ras2vox, verts)).astype(int)  # (V,3)

            #  (C) coronal slice 위 정점만 필터링
            on_slice    = cor_data[
                vox_verts[:, 0],
                vox_verts[:, 1],
                vox_verts[:, 2]
            ]  # boolean mask (V,)
            slice_verts = verts[on_slice]               # (V_slice, 3)
            slice_vids  = np.nonzero(on_slice)[0]       # 원래 verts 인덱스

            if len(slice_verts) == 0:
                raise RuntimeError(f"hemi={hemi}, roi={target_label_name}: no vertex on slice")

            #  (D) 슬라이스 전용 KD‑tree 생성
            slice_tree = cKDTree(slice_verts)
            dists_slice, idx_slice = slice_tree.query(world, workers=-1)
            # idx_slice 는 slice_verts 내의 인덱스 → 원래 verts 인덱스로 복원
            closest_vids  = slice_vids[idx_slice]
            region_ids    = vlabels[closest_vids]  # (N,)

            # ── (4) constrain mask ──────────────────────────────────────────────────
            n_ring = 2

            label_match = (vlabels == target_label_id)
            # n_ring까지 모두 같은 label이어야 통과
            mask_label = np.array([
                all(label_match[nbr] for nbr in get_n_ring_neighbors(vertex_adj, vid, n_ring))
                for vid in closest_vids
            ], dtype=bool)

            scale      = 1.8
            local_thresh = thickness[closest_vids] * scale  # shape (N,)
            mask_dist = (dists_slice <= local_thresh)
            # max_dist_mm = dist_lims[target_label_id]
            # mask_dist   = (dists_slice <= max_dist_mm)
            keep_mask  = mask_label & mask_dist
            kept_vox  = ijk[keep_mask]

            constrained = np.zeros(mask_data.shape, dtype=np.uint8)
            constrained[tuple(kept_vox.T)] = 1

            out_path = roi_dir / f"constrained_mask_ctx_{hemi}_{target_label_name}.mgz"
            nib.save(nib.MGHImage(constrained, affine=mask_img.affine, header=mask_img.header), out_path)

            # ── (5) report ──────────────────────────────────────────────────────────
            kept = np.count_nonzero(keep_mask)
            print(f"[info] total input voxels = {len(ijk):,d}")
            print(f"[info] kept voxels        = {kept:,d}")
            print(f"[info] output saved to    : {out_path}")

def main():
    table = load_slice_table(TXT_PATH)
    if not table:
        print('NewSliceNum.txt is empty')
        return

    for subj, slice_vals in table:
        print(f'\n===== {subj} start =====')
        try:
            run_pipeline(subj, slice_vals)
        except Exception as e:
            print(f'[error] {subj}: {e}')

if __name__ == '__main__':
    main()
