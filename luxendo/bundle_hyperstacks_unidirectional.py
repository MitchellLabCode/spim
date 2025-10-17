"""
Batch max-projection of XYT HDF5 files into Z-stacks.

- Each .h5 file must contain a 3D dataset at /Data (interpreted as XYT in any order).
- Groups every 'nzplanes' frames along the time axis into Z-blocks.
- Max-projects slices zstart..end of each block.
- Saves one multi-page TIFF per .h5 file.
"""

import sys
from pathlib import Path
import numpy as np
import h5py
import tifffile as tiff


def max_project_blocks_h5(in_h5, out_tif, nzplanes, zstart_1based=1, dset="/Data", xyt_axes=[1,2,0]):
    """Perform blockwise max projection on one HDF5 file."""
    z0 = zstart_1based - 1
    with h5py.File(in_h5, "r") as f:
        if dset not in f:
            raise KeyError(f"Dataset '{dset}' not found in {in_h5}")
        
        ds = f[dset]
        
        x_ax, y_ax, t_ax = xyt_axes
        
        T, Y, X = ds.shape[t_ax], ds.shape[y_ax], ds.shape[x_ax]
        nblocks = T // nzplanes
        if T % nzplanes:
            print(f"[WARN] {in_h5.name}: ignoring last {T % nzplanes} frame(s)")
        projs = np.empty((nblocks, X, Y), dtype=ds.dtype)
        base = [slice(None)] * 3
        for b in range(nblocks):
            t_start = b * nzplanes + z0
            t_stop = (b + 1) * nzplanes
            sl = list(base)
            sl[t_ax] = slice(t_start, t_stop)
            block = ds[tuple(sl)]
            zproj = block.max(axis=t_ax)
            # reorder to (Y, X)
            
            if (y_ax, x_ax) != (1, 2):
                axes = [0, 1]
                if y_ax > x_ax:
                    
                    zproj = np.moveaxis(zproj, (x_ax-1, y_ax-1), (0, 1))
                    
            projs[b] = zproj
    tiff.imwrite(out_tif, projs, imagej=True, metadata={'axes': 'TYX'})
    print(f"[OK] {in_h5.name} -> {out_tif}")

def batch_project_directory(in_dir, out_dir, nzplanes, zstart_1based=1, dset="/Data"):
    """Loop over all .h5 files in a directory and process each."""
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(in_dir.glob("*.h5"))
    if not files:
        print(f"[WARN] No .h5 files found in {in_dir}")
        return
    for f in files:
        out_tif = out_dir / f"{f.stem}_maxZ.tif"
        try:
            max_project_blocks_h5(f, out_tif, nzplanes, zstart_1based, dset)
        except Exception as e:
            print(f"[ERROR] {f.name}: {e}")

if __name__ == "__main__":
    # CLI usage:
    #   python batch_h5_to_zproj.py <input_dir> <output_dir> <nzplanes> [zstart_1based]
    if len(sys.argv) < 4:
        # print("Usage: python batch_h5_to_zproj.py <input_dir> <output_dir> <nzplanes> [zstart_1based]")
        # sys.exit(1)
        in_dir = 'F:\\PROJECTS\\LightMicroscopyBootcamp2025\\48YGAL4klar_UASGFPnlsUASmChCAAX\\2025-10-15_202209_overnight_timelapse_200reps\\raw\\stack_0_channel_0_obj_right'
        out_dir = 'F:\\PROJECTS\\LightMicroscopyBootcamp2025\\48YGAL4klar_UASGFPnlsUASmChCAAX\\2025-10-15_202209_overnight_timelapse_200reps\\mips\\'
        nzplanes = 23
        zstart = 6
    else:
        in_dir = sys.argv[1]
        out_dir = sys.argv[2]
        nzplanes = int(sys.argv[3])
        zstart = int(sys.argv[4]) if len(sys.argv) > 4 else 1
    
    batch_project_directory(in_dir, out_dir, nzplanes, zstart)
