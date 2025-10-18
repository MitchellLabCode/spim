#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Mask and background subtract a sequence of TIFF volumetric images.

This removes the outside (ex beads) while retaining an iLastik-identified ROI (embryo)
"""

# ==== USER CONFIG ============================================================
H5_PATH        = r"D:\\mblBootcamp\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined\\hindgut_avgAcrossTime_Merged_downsampled4x_Probabilities_float32.h5"
H5_DATASET     = "/exported_data"     # typical iLastik dataset path (adjust if needed)
ILASTIK_CH     = 0                    # which class prob to use (0 or 1 for 2-class iLastik)
THRESH         = 0.5                  # threshold on probability to create mask
INPUT_GLOB     = r"D:\\mblBootcamp\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined\\deconvolved18iter_15it\\TP*.tif"  # wildcards OK
OUTPUT_DIR     = r"D:\\mblBootcamp\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined\\masked_deconvolved_18it_15it_bgsub\\"                # will be created if missing
OUTPUT_MIP_DIR     = r"D:\\mblBootcamp\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined\\masked_max_projection_bgsub\\"                # will be created if missing
PRESERVE_DTYPE = True                 # cast back to each original TIFF dtype
PAD_VALUE      = 0                    # value to pad with if sizes off by a pixel at edges
# ---- USER scaling for background (y0:y1, x0:x1) ----
fudge_factor = 0.01

# ============================================================================

import os
import sys
import glob
from pathlib import Path
import numpy as np
import h5py
import tifffile as tiff
from scipy.ndimage import binary_fill_holes, binary_dilation, generate_binary_structure

# Try SciPy for interpolation; fall back to NumPy repeat if missing
try:
    from scipy.ndimage import zoom as ndi_zoom
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

def read_ilastik_prob_volume(h5_path, dataset, ch_index):
    with h5py.File(h5_path, "r") as f:
        data = f[dataset]
        arr = data[()]  # load to RAM (float32)
    # Expect (Z, C, Y, X). If shapes differ, try to infer.
    if arr.ndim != 4:
        raise ValueError(f"Expected 4D iLastik data (Z,C,Y,X). Got shape {arr.shape}")
    Z, C, Y, X = arr.shape
    if ch_index < 0 or ch_index >= C:
        raise IndexError(f"Channel index {ch_index} out of range 0..{C-1}")
    prob = arr[:, ch_index, :, :].astype(np.float32, copy=False)  # (Z, Y, X)
    return prob

def first_tiff_shape(input_glob):
    files = sorted(glob.glob(input_glob))
    if not files:
        raise SystemExit(f"No files matched: {input_glob}")
    with tiff.TiffFile(files[0]) as tf:
        shp = tf.series[0].shape
    if len(shp) != 3:
        raise ValueError(f"Expected 3D TIFF stacks. Found shape {shp} in {files[0]}")
    return tuple(shp), files

def fit_to_shape(arr, target_shape, pad_value=0):
    """
    Center-pad or center-crop arr (Z,Y,X) to target_shape.
    Padding uses constant pad_value. Cropping is symmetric.
    """
    z, y, x = arr.shape
    tz, ty, tx = target_shape
    out = arr

    # pad or crop Z
    if z < tz:
        pad_before = (tz - z) // 2
        pad_after  = tz - z - pad_before
        out = np.pad(out, ((pad_before, pad_after), (0,0), (0,0)),
                     mode="constant", constant_values=pad_value)
    elif z > tz:
        cut = (z - tz) // 2
        out = out[cut:cut+tz, :, :]

    # pad or crop Y
    z, y, x = out.shape
    if y < ty:
        pad_before = (ty - y) // 2
        pad_after  = ty - y - pad_before
        out = np.pad(out, ((0,0), (pad_before, pad_after), (0,0)),
                     mode="constant", constant_values=pad_value)
    elif y > ty:
        cut = (y - ty) // 2
        out = out[:, cut:cut+ty, :]

    # pad or crop X
    z, y, x = out.shape
    if x < tx:
        pad_before = (tx - x) // 2
        pad_after  = tx - x - pad_before
        out = np.pad(out, ((0,0), (0,0), (pad_before, pad_after)),
                     mode="constant", constant_values=pad_value)
    elif x > tx:
        cut = (x - tx) // 2
        out = out[:, :, cut:cut+tx]

    return out

def upsample_to_shape(prob_vol, target_shape, order=1, pad_value=0):
    """
    Upsample prob_vol (Z,Y,X) to target_shape with interpolation.
    Uses SciPy's ndimage.zoom if available (preferred). If not,
    uses nearest-neighbor via repeat (rough but dependency-free).
    Then center pad/crop to target_shape to fix any rounding.
    """
    sz, sy, sx = prob_vol.shape
    tz, ty, tx = target_shape

    if _HAS_SCIPY:
        # zoom factors (float); output will be close but may be off by 1
        factors = (tz / sz, ty / sy, tx / sx)
        # scipy.ndimage.zoom returns size ≈ round(in * factor)
        up = ndi_zoom(prob_vol, zoom=factors, order=order, prefilter=(order > 1))
    else:
        IOError('Better to use scipy! Please install it')
        # Nearest neighbor fallback using integer repeats per axis
        # compute integer repeats (at least 1)
        rz = max(1, int(round(tz / sz)))
        ry = max(1, int(round(ty / sy)))
        rx = max(1, int(round(tx / sx)))
        up = np.repeat(np.repeat(np.repeat(prob_vol, rz, axis=0), ry, axis=1), rx, axis=2)

    # Final tidy pad/crop to exact target shape
    up = fit_to_shape(up, target_shape, pad_value=pad_value)
    return up

def binarize(prob, thresh):
    mask = (prob >= float(thresh)).astype(np.float32, copy=False)  # 0/1 float mask
    return mask

def apply_mask_and_cast(vol, mask, preserve_dtype=True):
    """
    vol: ndarray (Z,Y,X), any dtype. mask: float32 0/1 same shape.
    Multiply in float32, then cast back if requested.
    """
    out = vol.astype(np.float32, copy=False) * mask
    if preserve_dtype:
        # convert back to original dtype with rounding/clipping if integer
        dt = vol.dtype
        if np.issubdtype(dt, np.integer):
            info = np.iinfo(dt)
            out = np.clip(np.rint(out), info.min, info.max).astype(dt, copy=False)
        elif np.issubdtype(dt, np.floating):
            out = out.astype(dt, copy=False)
    return out

def main():
    # 1) Read iLastik probabilities
    print("[INFO] Loading iLastik probability volume...")
    prob = read_ilastik_prob_volume(H5_PATH, H5_DATASET, ILASTIK_CH)  # (Z,Y,X)
    print(f"[INFO] iLastik prob shape: {prob.shape} (Z,Y,X)")

    # 2) Determine target stack shape from first TIFF
    target_shape, files = first_tiff_shape(INPUT_GLOB)
    print(f"[INFO] Target TIFF shape: {target_shape} (Z,Y,X)")
    print(f"[INFO] Matched {len(files)} stack(s).")

    # 3) Upsample probs to target shape
    print("[INFO] Upsampling probability to target shape...")
    prob_up = upsample_to_shape(prob, target_shape, order=1, pad_value=PAD_VALUE)

    # 4) Threshold to binary mask
    print(f"[INFO] Thresholding at {THRESH} ...")
    mask = binarize(prob_up, THRESH)  # float32 0/1

    # Dilate a bit to catch small boundary holes/gaps (≈ Fiji imdilate)
    st = generate_binary_structure(3, 1)  # 6-connected; use (3,2) for thicker dilation
    mask = binary_dilation(mask, structure=st, iterations=1)

    # -- NEW: fill internal holes in the mask --
    mask = binary_fill_holes(mask).astype(np.float32)
    print("[INFO] Filled holes in mask using scipy.ndimage.binary_fill_holes")

    # 5) Ensure output directory
    outdir = Path(OUTPUT_DIR)
    outdir.mkdir(parents=True, exist_ok=True)

    # 5b) Ensure output directory
    mipdir = Path(OUTPUT_MIP_DIR)
    mipdir.mkdir(parents=True, exist_ok=True)

    # 6) Apply mask to each stack and write
    n = len(files)
    for i, fp in enumerate(files, 1):
        print(f"\r[INFO] Processing {i}/{n}: {fp}", end="", flush=True)
        vol = tiff.imread(fp)  # (Z,Y,X)

        # Ensure mask matches this stack's shape (rarely needed but safe)
        m = fit_to_shape(mask, vol.shape, pad_value=PAD_VALUE) if vol.shape != mask.shape else mask

        # ---- Background on Z-max projection of ORIGINAL volume ----
        # ---- compute bg_min from inside-mask only (robust to outside zeros) ----
        mip_orig = np.max(vol, axis=0)  # (Y, X)
        mip_mask = np.max(m, axis=0) > 0  # where object projects in MIP
        inside_vals = mip_orig[mip_mask]
        if inside_vals.size == 0:
            raise ValueError("Mask is empty over this stack; cannot compute baseline.")

        # FILTER FOR STRICTLY NONZERO VALUES
        nonzero_inside = inside_vals[inside_vals > 0]

        if nonzero_inside.size == 0:
            print("WARNING: No nonzero inside-mask values found; falling back to minimal positive constant")
            bg_min = 1e-3
            Exception("No nonzero inside-mask values found; This should be an error")
        else:
            bg_min = float(nonzero_inside.min())  # <-- smallest NONZERO value

        print('bg_min = ' + str(bg_min))

        # ---- zero masked voxels, scale whole volume by bg_min (so bg_min → ~1) ----
        work = vol.astype(np.float32, copy=False)
        work /= (bg_min / fudge_factor)
        work[m == 0] = fudge_factor

        # ---- Subtract the smallest nonzero value FROM INSIDE THE MASK ONLY ----
        # (this aligns the dimmest in-mask background to 0 while keeping outside masked at 0)
        in_mask = (m > 0)
        nz_in_mask = in_mask & (work > 0)
        if np.any(nz_in_mask):
            min_nz = float(work[nz_in_mask].min())
            # Subtract this baseline so the minimum nonzero becomes zero ----
            work -= min_nz
        else:
            min_nz = 0.0

        # ---- Clip any tiny negatives from float precision ----
        np.maximum(work, 0.0, out=work)

        # ---- Optional cast back to original dtype ----
        if PRESERVE_DTYPE:
            dt = vol.dtype
            if np.issubdtype(dt, np.integer):
                info = np.iinfo(dt)
                work = np.clip(np.rint(work), info.min, info.max).astype(dt, copy=False)
            else:
                work = work.astype(dt, copy=False)

        # ---- Write scaled, masked stack ----
        base = os.path.basename(fp)
        out_path = outdir / base
        tiff.imwrite(str(out_path), work, bigtiff=False, compression="zlib", metadata=None)

        # ---- Write Z-max MIP of the scaled stack (clear it's a MIP) ----
        mip = np.max(work, axis=0)  # (Y, X)
        mip_path = mipdir / (Path(base).stem + "_MIP.tiff")
        tiff.imwrite(str(mip_path), mip, bigtiff=False, compression="zlib", metadata=None)

    print("\n[OK] Done.")


if __name__ == "__main__":
    main()
