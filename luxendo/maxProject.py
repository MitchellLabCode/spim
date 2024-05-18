import os
import h5py
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


# Function to save MIP images

def save_mip_images(mip, axis, filename, output_dir):
    output_path = os.path.join(output_dir, f"{filename}_mip_axis_{axis}.tiff")
    print('-->' + output_path)
    Image.fromarray(mip).save(output_path)


# Function to process a single directory

def process_directory(input_dir, output_dir, plane=-1):
    # Go through all directories in the input_dir directory, make MIPS of all subdirectories
    #
    # Parameters
    # ----------
    # input_dir : str
    # output_dir : str
    # plane : int (default=-1)
    #   if >0, make record of single images of this z-plane
    #   if <0, use middle of the zstack as the plane of interest
    #
    os.makedirs(output_dir, exist_ok=True)
    for filename in os.listdir(input_dir):
        if filename.endswith('.h5'):
            filepath = os.path.join(input_dir, filename)
            print('-> ' + filepath)
            with h5py.File(filepath, 'r') as f:
                data = f['/Data'][()]
            mip_x = np.max(data, axis=0)
            # mip_y = np.max(data, axis=1)
            # mip_z = np.max(data, axis=2)
            base_filename = os.path.splitext(filename)[0]
            save_mip_images(mip_x, 'x', base_filename, output_dir)
            # save_mip_images(mip_y, 'y', base_filename, output_dir)
            # save_mip_images(mip_z, 'z', base_filename, output_dir)

            # Subset of z-stack as mip (mid-saggital, for ex)
            if plane < 0:
                # take middle z plane
                plane = np.round(np.size(data, 0) * 0.5).astype(np.int64)

            mip_x = data[plane]
            os.makedirs(output_dir + '_midZ', exist_ok=True)
            save_mip_images(mip_x, 'x', base_filename, output_dir + '_midZ')


if __name__ == "__main__":

    # parent_dir = 'D:\\rlondo\\HandGFPbynGAL4klar_H2AGFPUASMyo1CRFP\\2024-05-15_230441'  # Replace with your parent directory path
    parent_dir = 'D:\\rlondo\\HandGFPbynGAL4klar_UASmChCAAXH2AviRFP\\2024-05-17_184833'  # Replace with your parent directory path

    data_dir = os.path.join(parent_dir, 'raw')
    outdir = os.path.join(parent_dir, 'mips')
    for dir_name in os.listdir(data_dir):
        input_dir = os.path.join(data_dir, dir_name)
        print(input_dir)
        if os.path.isdir(input_dir):
            output_dir = os.path.join(outdir, dir_name)
            print('processing ' + input_dir + " -> " + output_dir)
            process_directory(input_dir, output_dir)
            print('done with' + input_dir)

    print('done')
    print('Now drag all mips subdirectories into Fiji. Run maxProjectMovies_simple.ijm after changing path in line 5')

