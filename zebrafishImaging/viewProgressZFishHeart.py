import os
import numpy as np
import h5py
import tifffile
import re


def process_h5_files(base_directory, timestamp, output_tiff, frame_indices=[0, 50, 100, 150], leftRight="right"):
    """
    Process HDF5 files in subdirectories ending with the specified leftRight suffix and generate a 3D TIFF file
    containing max intensity projections (MIP) from specified frames.

    Parameters
    ----------
    base_directory : str
        Path to the base directory containing subdirectories with HDF5 files.
    timestamp : str
        Timestamp to identify the specific HDF5 file (e.g., "00000").
    output_tiff : str
        Path to save the output 3D TIFF file.
    frame_indices : list of int, optional
        List of frame indices to extract and compute the MIP (default is [0, 50, 100, 150]).
    leftRight : str, optional
        Specifies whether to process "left" or "right" images (default is "right").

    Returns
    -------
    None
    """
    mip_stack = []

    # Identify subdirectories with the expected pattern
    subdirectories = [d for d in os.listdir(base_directory) if
                      d.startswith("stack_") and d.endswith(f"_{leftRight}") and os.path.isdir(
                          os.path.join(base_directory, d))]

    for subdir in sorted(subdirectories):
        h5_path = os.path.join(base_directory, subdir, f"Cam_{leftRight}_{timestamp}.lux")

        if os.path.exists(h5_path):
            with h5py.File(h5_path, "r") as h5_file:
                # Assuming dataset is named "data"
                dataset_name = list(h5_file.keys())[0]  # Gets the first dataset
                data = h5_file[dataset_name]

                # Extract the specified frames and compute MIP without loading entire dataset
                frames = np.stack([data[i] for i in frame_indices if i < data.shape[0]], axis=0)
                mip = np.max(frames, axis=0)
                mip_stack.append(mip)
        else:
            print(f"Warning: {h5_path} not found.")

    # Convert list to numpy array and save as TIFF
    if mip_stack:
        mip_stack = np.array(mip_stack)
        tifffile.imwrite(output_tiff, mip_stack.astype(np.uint16))
        print(f"Saved MIP stack to {output_tiff}")
    else:
        print("No valid data found.")


# Identify unique timestamps from available HDF5 files
base_directory = "E:/sid/2025-02-27_183701/raw/"
timestamp_set = set()

for subdir in os.listdir(base_directory):
    subdir_path = os.path.join(base_directory, subdir)
    if os.path.isdir(subdir_path) and subdir.endswith("_right"):
        for file in os.listdir(subdir_path):
            match = re.match(r"Cam_right_(\d+).lux", file)
            if match:
                timestamp_set.add(int(match.group(1)))

# Convert to sorted numpy array
timestamps = np.array(sorted(timestamp_set), dtype=int)

# Process files for each timestamp
for timestamp in timestamps:
    output_tiff = f"output_mip_stack_{timestamp:05d}.tif"
    process_h5_files(base_directory, f"{timestamp:05d}", output_tiff)
