import os
import re
import tifffile
import numpy as np
from PIL import Image
from pathlib import Path


def process_tiff_directory(
        input_dir,
        output_dir,
        clipZ=(0, 0),
        clipX=(0, 0),
        clipY=(0, 0),
        process_in_reverse=True
):
    """
    Take MIPs along z for all OME-TIFF files in a single directory.

    Assumes data is ZXY (Z, X, Y).

    Files are expected to be named like:
        tp{tttttt}_c{c}_a{a}.ome.tif
    and MIPs are written into angle-specific subfolders:
        output_dir/mips_z/a{a}/...
    """
    os.makedirs(output_dir, exist_ok=True)
    mip_dir = output_dir   # os.path.join(output_dir, 'mips_z')
    os.makedirs(mip_dir, exist_ok=True)

    # Get all tif/ome.tif files in the directory
    filenames = [f for f in os.listdir(input_dir)
                 if f.lower().endswith('.tif') or f.lower().endswith('.tiff')]
    filenames = sorted(filenames)

    print(filenames)

    if process_in_reverse:
        filenames = filenames[::-1]

    for filename in filenames:
        # Enforce naming pattern (tp... .ome.tif)
        if not (filename.startswith('t') and '.ome.tif' in filename):
            print('skipping file ' + filename)
            continue

        # Extract angle from filename: tpXXXXXX_cY_aZ.ome.tif
        m = re.search(r'_a(\d+)\.ome\.tif$', filename)
        if m is None:
            print(f"Could not parse angle from filename, skipping: {filename}")
            continue
        angle_str = m.group(1)  # e.g. "0", "1", "15", ...

        filepath = os.path.join(input_dir, filename)
        base_filename = os.path.splitext(filename)[0]  # tpXXXXXX_cY_aZ.ome

        # Angle-specific output subfolder
        angle_dir = os.path.join(mip_dir, f"a{angle_str}")
        os.makedirs(angle_dir, exist_ok=True)

        out_path = os.path.join(angle_dir, f"{base_filename}_mipZ.tif")

        if os.path.exists(out_path):
            print(f"--> Output file already exists, skipping: {out_path}")
            continue

        print('-> Reading', filepath)
        data = tifffile.imread(filepath)

        # Expect ZXY (Z, X, Y)
        print('   data shape:', data.shape)

        # Apply clipping
        z0, z1 = clipZ
        x0, x1 = clipX
        y0, y1 = clipY

        # Handle 0 clipping gracefully
        z1 = data.shape[0] - z1 if z1 > 0 else data.shape[0]
        x1 = data.shape[1] - x1 if x1 > 0 else data.shape[1]
        y1 = data.shape[2] - y1 if y1 > 0 else data.shape[2]

        data_clipped = data[z0:z1, x0:x1, y0:y1]

        # Max projection along z (axis 0)
        mip_z = np.max(data_clipped, axis=0)

        print('--> Saving MIP to', out_path)
        Image.fromarray(mip_z).save(out_path)

def build_hyperstacks_from_mips(mips_root_dir):
    """
    Scrape all MIP images in each angle subfolder and build
    ImageJ/Fiji-compatible hyperstacks (T x C x Y x X) per angle.

    Expects files like:
        tp{tttttt}_c{c}_a{a}.ome_mipZ.tif
    stored in:
        mips_root_dir/mips_z/a{a}/
    """
    mips_root_dir = str(mips_root_dir)
    mips_z_dir = mips_root_dir  # os.path.join(mips_root_dir, 'mips_z')
    if not os.path.isdir(mips_z_dir):
        print(f"No mips_z directory found at {mips_z_dir}")
        return

    # Angle subfolders: a0, a1, ...
    angle_dirs = [d for d in os.listdir(mips_z_dir)
                  if os.path.isdir(os.path.join(mips_z_dir, d)) and d.startswith('a')]

    for angle_dir in sorted(angle_dirs):
        angle_path = os.path.join(mips_z_dir, angle_dir)
        print(f"\nProcessing angle folder: {angle_path}")

        files = [f for f in os.listdir(angle_path)
                 if f.lower().endswith('.tif') or f.lower().endswith('.tiff')]
        if not files:
            print(f"  No TIFFs found in {angle_path}, skipping.")
            continue

        # Parse time and channel from filenames
        # Pattern: tp000000_c0_a0.ome_mipZ.tif
        pattern = re.compile(r'^(t(\d+)_c(\d+)_a(\d+)\.ome)_mipZ\.tif$')

        entries = []
        for f in files:
            m = pattern.match(f)
            if m is None:
                print(f"  Skipping file that does not match pattern: {f}")
                continue
            base_name = m.group(1)
            t_str = m.group(2)
            c_str = m.group(3)
            a_str = m.group(4)  # should be consistent within this folder

            t = int(t_str)
            c = int(c_str)
            entries.append((t, c, f))

        if not entries:
            print(files)
            print(pattern)
            RuntimeError(f"  No files matched pattern in {angle_path}")
            print(f"  No files matched pattern in {angle_path}, skipping.")
            continue

        # Sort and get unique times and channels
        entries.sort(key=lambda x: (x[0], x[1]))
        times = sorted({e[0] for e in entries})
        channels = sorted({e[1] for e in entries})

        t_index = {t: i for i, t in enumerate(times)}
        c_index = {c: i for i, c in enumerate(channels)}

        # Read one image to get spatial shape and dtype
        sample_file = os.path.join(angle_path, entries[0][2])
        sample_img = tifffile.imread(sample_file)
        height, width = sample_img.shape
        dtype = sample_img.dtype

        nt = len(times)
        nc = len(channels)

        print(f"  Found {nt} timepoints, {nc} channels, image size {width}x{height}")

        # Allocate hyperstack array: T x C x Y x X
        stack = np.zeros((nt, nc, height, width), dtype=dtype)

        # Fill array
        for t, c, fname in entries:
            ti = t_index[t]
            ci = c_index[c]
            img_path = os.path.join(angle_path, fname)
            img = tifffile.imread(img_path)
            print(np.size(img))
            stack[ti, ci, :, :] = img

        # Save hyperstack TIFF (ImageJ compatible)
        out_name = f"{angle_dir}_hyperstack_TC.tif"
        out_path = os.path.join(mips_z_dir, out_name)

        print(f"  Saving hyperstack to {out_path}")
        tifffile.imwrite(
            out_path,
            stack,
            imagej=True,
            metadata={'axes': 'TCYX'}
        )


if __name__ == "__main__":
    # Directory that directly contains the tpXXXXXX_cY_aZ.ome.tif files
    input_dir = Path(
        r'E:\avistrok\bapGAL4_UAShidUASstingerHiRFP\20250806120735_bapGAL4_UAShidUASStingerHiRFP_5mpf_22x\20250806121230_bapGAL4_UAShidUASStingerHiRFP_5mpf_22x\unpacked'
    )

    # Where to put the mips_z folder
    output_dir = os.path.join(input_dir, 'mips')

    # Clipping (Z, X, Y): [low, high]
    clipZ = (0, 0)
    clipX = (0, 0)
    clipY = (0, 0)

    process_tiff_directory(
        input_dir=input_dir,
        output_dir=output_dir,
        clipZ=clipZ,
        clipX=clipX,
        clipY=clipY,
        process_in_reverse=False
    )

    # Build Fiji-ready hyperstacks from the MIPs
    build_hyperstacks_from_mips(output_dir)

    print('done')
