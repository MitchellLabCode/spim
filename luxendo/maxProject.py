import os
import h5py
import tifffile
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image


# Function to save MIP images
# def save_mip_images(mip, axis, filename, output_dir):
#     output_path = os.path.join(output_dir, f"{filename}_mip_axis_{axis}.tiff")
#     print('--> Saving to ' + output_path)
#     Image.fromarray(mip).save(output_path)


# Function to process a single directory
def process_directory(input_dir, output_dir, dir_name, plane=-1):
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
    max_proj_dir = os.path.join(output_dir, 'mips')
    os.makedirs(max_proj_dir, exist_ok=True)
    midz_output_dir = os.path.join(output_dir, 'midZ')
    os.makedirs(midz_output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        if filename.endswith('.h5'):
            base_filename = os.path.splitext(filename)[0]

            os.makedirs(os.path.join(max_proj_dir, dir_name+"_mips"), exist_ok=True)
            os.makedirs(os.path.join(midz_output_dir, dir_name+"_midZ"), exist_ok=True)
            max_proj_path = os.path.join(max_proj_dir, dir_name+"_mips", f"{base_filename}_mip_axis_x.tiff")
            midz_path = os.path.join(midz_output_dir, dir_name+"_midZ", f"{base_filename}_midZ_axis_x.tiff")

            if not os.path.exists(max_proj_path) or not os.path.exists(midz_path):
                filepath = os.path.join(input_dir, filename)
                print('-> ' + filepath)
                with h5py.File(filepath, 'r') as f:
                    data = f['/Data'][()]
                mip_x = np.max(data, axis=0)
                # mip_y = np.max(data, axis=1)
                # mip_z = np.max(data, axis=2)
                Image.fromarray(mip_x).save(max_proj_path)
                # save_mip_images(mip_y, 'y', base_filename, output_dir)
                # save_mip_images(mip_z, 'z', base_filename, output_dir)

                # Subset of z-stack as mip (mid-saggital, for ex)
                if plane < 0:
                    # take middle z plane
                    plane = np.round(np.size(data, 0) * 0.5).astype(np.int64)

                midz = data[plane]
                Image.fromarray(midz).save(midz_path)

            else:
                print(f"--> Output file already exists, skipping: {max_proj_path}")


if __name__ == "__main__":

    # parent_dir = 'D:\\rlondo\\HandGFPbynGAL4klar_H2AGFPUASMyo1CRFP\\2024-05-15_230441'  # Replace with your parent directory path
    # parent_dir = 'D:\\rlondo\\HandGFPbynGAL4klar_UASmChCAAXH2AviRFP\\2024-05-17_184833'  # Replace with your parent directory path
    # parent_dir = 'D:\\rlondo\\HandGFPbynGAL4klar_UASmChCAAXH2AviRFP\\2024-05-17_184833'  # Replace with your parent directory path

    # parent_dir = 'E:\\haibei\\48YGAL4klar_UASmChCAAXHiRFP\\2024-05-23_183541'
    # parent_dir = 'E:\\rio\\HandGFPbynGAL4klar_UASmChCAAXHiFP\\2024-05-27_125722_22C_butFocusDrifts'
    # parent_dir = 'E:\\rio\\HandGFPbynGAL4klar_UASMyo1CRFPHiFP\\2024-05-28_150855'
    # parent_dir = 'E:\\rio\\HandGFPbynGAL4klar_UASMyo1CRFPHiFP\\2024-05-28_150855_excellent'
    # parent_dir = 'E:\\rio\\HandGFPbynGAL4klar_UASmChCAAXHiFP\\2024-05-31_142410'
    parent_dir = 'E:\\haibei\\48YGAL4klar_UASLamGFPUASmChCAAX\\2024-06-03_142916'

    data_dir = os.path.join(parent_dir, 'raw')
    for dir_name in os.listdir(data_dir):
        input_dir = os.path.join(data_dir, dir_name)
        print(input_dir)
        if os.path.isdir(input_dir):
            output_dir = parent_dir  # os.path.join(outdir, 'mips', dir_name)
            print('processing ' + input_dir + " -> " + output_dir + "\\mips or midZ")
            process_directory(input_dir, output_dir, dir_name)
            print('done with ' + input_dir)

    print('done')
    print('Now drag all mips subdirectories into Fiji. Run maxProjectMovies_simple.ijm after changing path in line 5')
