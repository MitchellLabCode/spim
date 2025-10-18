#!/usr/bin/env python3
"""
Average a series of 3D TIFF volumes and write the mean as TIFF or HDF5.

Usage examples:
  # Average all .tif in a folder and write a BigTIFF
  python avg_stacks.py /path/to/folder/*.tif -o /path/to/output/mean.tif

  # Average using z-slice chunking (lower RAM), ignore NaNs, and write HDF5
  python avg_stacks.py "/data/stacks/*.tif" -o /data/mean.h5 --h5 --dataset /Data --chunked --ignore-nan
"""


# ==== USER CONFIG (EDIT THIS ONLY) ==========================================
INPUT_GLOB = "D:\\mblBootcamp\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined\\deconvolved18iter_15it\\TP*Ch0*.tif"     # wildcards OK
OUTPUT_PATH = "D:\\mblBootcamp\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined\\ch0_avg_uint16.tiff"  # .tif or .h5
CHUNKED = True                            # True = memory safe z-slice summing
IGNORE_NAN = False                        # True = NaN-aware averaging
TARGET_DTYPE = None                       # e.g. "float32" or "uint16" or None
# ============================================================================

import argparse
from pathlib import Path
import numpy as np
import tifffile as tiff
import h5py
import glob

def read_shape_dtype_tiff(fn):
    with tiff.TiffFile(str(fn)) as tf:
        arr_shape = tf.series[0].shape  # expect (Z, Y, X) or (Y, X, Z) if weird—handled below
        dtype = tf.series[0].dtype
    return arr_shape, dtype

def load_slice_tiff(fn, z):
    # Read only one Z plane from a 3D stack
    return tiff.imread(str(fn), key=int(z))  # shape (Y, X)

def average_tiff_series(filepaths, chunked=False, ignore_nan=False, target_dtype=None, read_mode="by_file"):
    if len(filepaths) == 0:
        raise ValueError("No input files provided.")

    # Inspect shape/dtype from the first file
    shape0, _ = read_shape_dtype_tiff(filepaths[0])
    Z, Y, X = shape0

    # Verify shapes match
    for fp in filepaths[1:]:
        shp, _ = read_shape_dtype_tiff(fp)
        if shp != shape0:
            raise ValueError(f"Shape mismatch: {fp} has {shp}, expected {shape0}")

    # --- Prefer file-major reading for speed on large datasets ---
    if read_mode == "by_file":
        sum_vol = np.zeros(shape0, dtype=np.float64)
        count_vol = np.zeros(shape0, dtype=np.uint32) if ignore_nan else None

        n = len(filepaths)
        for i, fp in enumerate(filepaths, 1):
            print(f"\rReading stack {i}/{n}: {fp}", end="", flush=True)
            vol = tiff.imread(str(fp)).astype(np.float64, copy=False)

            if ignore_nan:
                mask = ~np.isnan(vol)
                sum_vol[mask] += vol[mask]
                count_vol[mask] += 1
            else:
                sum_vol += vol

        print()  # newline after progress
        if ignore_nan:
            with np.errstate(invalid="ignore", divide="ignore"):
                mean_vol = sum_vol / np.maximum(count_vol, 1)
                mean_vol[count_vol == 0] = np.nan
        else:
            mean_vol = sum_vol / float(len(filepaths))

    else:
        # Fallback: slice-major (older approach). Keeps peak RAM low but slower I/O.
        sum_vol = np.zeros((Z, Y, X), dtype=np.float64)
        count_vol = np.zeros((Z, Y, X), dtype=np.uint32) if ignore_nan else None

        for z in range(Z):
            print(f"\rAveraging slice {z+1}/{Z} ...", end="", flush=True)
            if ignore_nan:
                plane_sum = np.zeros((Y, X), dtype=np.float64)
                plane_cnt = np.zeros((Y, X), dtype=np.uint32)
                for fp in filepaths:
                    slc = load_slice_tiff(fp, z).astype(np.float64, copy=False)
                    mask = ~np.isnan(slc)
                    plane_sum[mask] += slc[mask]
                    plane_cnt[mask] += 1
                sum_vol[z] = plane_sum
                count_vol[z] = plane_cnt
            else:
                plane_sum = np.zeros((Y, X), dtype=np.float64)
                for fp in filepaths:
                    slc = load_slice_tiff(fp, z).astype(np.float64, copy=False)
                    plane_sum += slc
                sum_vol[z] = plane_sum

        print()
        if ignore_nan:
            with np.errstate(invalid="ignore", divide="ignore"):
                mean_vol = sum_vol / np.maximum(count_vol, 1)
                mean_vol[count_vol == 0] = np.nan
        else:
            mean_vol = sum_vol / float(len(filepaths))

    # Cast output
    if target_dtype is None:
        out = mean_vol.astype(np.float32, copy=False)
    else:
        td = np.dtype(target_dtype)
        out = np.rint(mean_vol).astype(td, copy=False) if np.issubdtype(td, np.integer) else mean_vol.astype(td, copy=False)
    return out


def write_tiff(out_path: Path, volume: np.ndarray):
    # BigTIFF automatically used for large data; use compression to save space
    tiff.imwrite(
        str(out_path),
        volume,
        bigtiff=True,
        compression="zlib",
        metadata=None,
    )

def write_h5(out_path: Path, volume: np.ndarray, dataset="/Data", gzip=4, chunks=True):
    out_path = Path(out_path)
    with h5py.File(out_path, "w") as h5:
        dset = h5.create_dataset(
            dataset,
            data=volume,
            compression="gzip" if gzip else None,
            compression_opts=gzip if gzip else None,
            chunks=chunks,
        )
        dset.attrs["axes"] = "zyx"
        dset.attrs["dtype_original"] = str(volume.dtype)

def scale_to_uint16(vol, vmin=None, vmax=None, percentiles=(0.5, 99.5)):
    """
    Scale a float/any-volume to uint16.
    - If percentiles is not None, compute (vmin, vmax) from those percentiles.
    - Else use provided vmin/vmax; if missing, fall back to nanmin/nanmax.
    - NaNs become 0.
    Returns (vol_u16, vmin, vmax).
    """
    import numpy as np
    x = vol.astype(np.float64, copy=False)

    if percentiles is not None:
        vmin, vmax = np.nanpercentile(x, percentiles)
    else:
        if vmin is None:
            vmin = np.nanmin(x)
        if vmax is None:
            vmax = np.nanmax(x)

    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        # Degenerate case: all equal/NaN → zeros
        return np.zeros_like(vol, dtype=np.uint16), vmin, vmax

    x = (x - vmin) * (65535.0 / (vmax - vmin))
    x = np.clip(np.rint(x), 0, 65535).astype(np.uint16)
    # Preserve mask: set NaNs to 0 in the output
    if np.isnan(vol).any():
        x[np.isnan(vol)] = 0
    return x, vmin, vmax


def parse_args():
    p = argparse.ArgumentParser(description="Average 3D TIFF volumes voxel-wise.")
    p.add_argument("inputs", nargs="+", help="Input TIFF(s). Glob patterns accepted (Put them in double quotes!).")
    p.add_argument("-o", "--output", required=True, help="Output file path (.tif/.tiff or .h5).")
    p.add_argument("--h5", action="store_true", help="Write HDF5 instead of TIFF (inferred from extension if omitted).")
    p.add_argument("--dataset", default="/Data", help="HDF5 dataset path (default: /Data).")
    p.add_argument("--chunked", action="store_true", help="Process one Z-slice at a time to reduce RAM.")
    p.add_argument("--ignore-nan", action="store_true", help="Average while ignoring NaNs (NaN where all inputs are NaN).")
    p.add_argument("--dtype", default=None, help="Output dtype (e.g., float32, uint16). Default: float32.")
    return p.parse_args()

def main():
    args = parse_args()

    # Expand globs and sort for reproducibility
    files = []
    for pattern in args.inputs:
        files.extend(sorted(Path().glob(pattern)))
    files = sorted({Path(f) for f in files})  # unique & sorted

    if not files:
        raise SystemExit("No input files matched.")

    target_dtype = np.dtype(args.dtype) if args.dtype is not None else None

    mean_vol = average_tiff_series(
        files,
        chunked=args.chunked,
        ignore_nan=args.ignore_nan,
        target_dtype=target_dtype,
    )

    out_path = Path(args.output)
    to_h5 = args.h5 or out_path.suffix.lower() in (".h5", ".hdf5")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    if to_h5:
        write_h5(out_path, mean_vol, dataset=args.dataset)
        print(f"Wrote HDF5 mean volume to: {out_path} (dataset: {args.dataset})")
    else:
        write_tiff(out_path, mean_vol)
        print(f"Wrote TIFF mean volume to: {out_path}")

if __name__ == "__main__":

    # --- Load file list ---
    files = sorted(glob.glob(INPUT_GLOB))
    if not files:
        raise SystemExit(f"No files matched: {INPUT_GLOB}")

    # --- Perform average ---
    mean_vol = average_tiff_series(
        files,
        chunked=CHUNKED,
        ignore_nan=IGNORE_NAN,
        target_dtype=np.dtype(TARGET_DTYPE) if TARGET_DTYPE else None,
    )

    # --- Scale to uint16 (robust stretch) ---
    mean_u16, vmin, vmax = scale_to_uint16(mean_vol, percentiles=(0, 100))
    print(f"[INFO] Scaling to uint16 with vmin={vmin:.6g}, vmax={vmax:.6g} (percentiles 0.5–99.5)")

    # --- Choose output format based on OUTPUT_PATH extension ---
    out_path = Path(OUTPUT_PATH)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if out_path.suffix.lower() in (".h5", ".hdf5"):
        write_h5(out_path, mean_u16, dataset=H5_DATASET)  # stores as uint16
        print(f"[OK] Wrote HDF5 mean volume: {out_path} (uint16, dataset {H5_DATASET})")
    else:
        write_tiff(out_path, mean_u16)  # tifffile infers dtype from array
        print(f"[OK] Wrote TIFF mean volume: {out_path} (uint16)")
