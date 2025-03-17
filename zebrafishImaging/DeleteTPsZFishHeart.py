import os
import numpy as np
import h5py
import tifffile
import re
import shutil


def del_timepoints(in_dir, out_dir, timepoint):
    file_name_left = f"Cam_left_{timepoint:05d}.lux.h5"  # e.g., "Cam_left_00001.lux.h5"
    file_name_right = f"Cam_right_{timepoint:05d}.lux.h5"  # e.g., "Cam_right_00001.lux.h5"
    for root, dirs, files in os.walk(in_dir):
        # Check if either of the files exist in the current directory
        if file_name_left in files:
            # Move the file to the output directory
            src_file = os.path.join(root, file_name_left)
            shutil.move(src_file, os.path.join(out_dir, file_name_left))
            print(f"Moved: {src_file} to {out_dir}")

        if file_name_right in files:
            # Move the file to the output directory
            src_file = os.path.join(root, file_name_right)
            shutil.move(src_file, os.path.join(out_dir, file_name_right))
            print(f"Moved: {src_file} to {out_dir}")

    print("Process completed.")


if __name__ == "__main__":
    in_dir = 'E:/sid/20250214_333fps_h2avGFPmemCh/2025-02-14_175827/raw'  # input directory path
    out_dir = 'E:/sid/20250214_333fps_h2avGFPmemCh/To_be_deleted'  # output directory path
    timepoint = 4  # The timepoint to be deleted
    del_timepoints(in_dir, out_dir, timepoint)


