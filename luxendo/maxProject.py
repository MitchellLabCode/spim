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

def process_directory(input_dir, output_dir):
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


if __name__ == "__main__":

    parent_dir = 'D:\\rlondo\\HandGFPbynGAL4klar_H2AGFPUASmChCAAX\\2024-05-15_230441'  # Replace with your parent directory path
    data_dir = os.path.join(parent_dir, 'raw')
    outdir = os.path.join(parent_dir, 'mips')
    for dir_name in os.listdir(data_dir):
        input_dir = os.path.join(data_dir, dir_name)
        print(input_dir)
        if os.path.isdir(input_dir):
            output_dir = os.path.join(outdir, dir_name)
            print('processing ' + input_dir + " -> " + output_dir)
            process_directory(input_dir, output_dir)
