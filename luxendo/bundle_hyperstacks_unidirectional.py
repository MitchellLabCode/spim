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
    """Perform blockwise max projection on one HDF5 file.

    The function groups consecutive frames along the time axis into blocks of
    size ``nzplanes`` and computes a max projection over the frames in each
    block. The result is saved as a multi-page TIFF where each page is the
    max-projection of one block.

    Parameters
    ----------
    in_h5 : str or path-like
        Path to the input HDF5 file. The file is opened read-only inside the
        function.
    out_tif : str or path-like
        Path for the output TIFF file to write. A multi-page TIFF is written
        with axis ordering compatible with ImageJ ('TYX').
    nzplanes : int
        Number of consecutive time frames that form one Z-block to be
        collapsed (max-projected) into a single Z-slice in the output.
    zstart_1based : int, optional
        1-based start offset applied to the first frame of each block. For
        example, ``zstart_1based=1`` (default) means blocks start at the
        nominal first frame; ``zstart_1based=2`` shifts the start by one frame.
    dset : str, optional
        HDF5 dataset path inside the file to read image data from (default
        ``"/Data"``). The dataset must be 3-dimensional.
    xyt_axes : list of three ints, optional
        Mapping of dataset axes to X, Y and T axes. Provide a list of three
        0-based axis indices in the order ``[x_axis, y_axis, t_axis]`` so the
        code can index into ``ds.shape`` correctly. Default ``[1, 2, 0]``
        assumes the dataset ordering where axis 0 is time, axis 1 is X and
        axis 2 is Y (adjust if your dataset uses a different layout).

    Returns
    -------
    None
        Writes the TIFF file at ``out_tif``. Raises ``KeyError`` if ``dset`` is
        not found in the HDF5 file or other exceptions on I/O errors.
    """
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
    """Process all .h5 files in a directory and write Z-projected TIFFs.

    Parameters
    ----------
    in_dir : str or path-like
        Directory containing input ``.h5`` files to process.
    out_dir : str or path-like
        Directory where output TIFF files will be written. The directory is
        created if it does not exist.
    nzplanes : int
        Number of frames to group into each block for max projection.
    zstart_1based : int, optional
        1-based offset to use as the start index within each block (default
        is 1). Passed through to ``max_project_blocks_h5``.
    dset : str, optional
        Dataset path inside each HDF5 file to read image data from
        (default ``"/Data"``). Passed through to ``max_project_blocks_h5``.

    Returns
    -------
    None
        Writes one TIFF per input HDF5 file. Errors for individual files are
        printed but do not stop processing of the directory.
    """
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
