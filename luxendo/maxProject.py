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
def process_directory(input_dir, output_dir, dir_name, plane=-1, clipZ=[0,0], clipX=[0,0], clipY=[0,0],
                      planeX=[-1], do_mips=True, do_midZ=True, do_midX=False):
    # Go through all directories in the input_dir directory, make MIPS of all subdirectories
    # Note data is in format ZXY.
    #
    # Parameters
    # ----------
    # input_dir : str
    # output_dir : str
    # plane : int (default=-1)
    #   if >0, make record of single images of this z-plane for midZ
    #   if <0, use middle of the zstack as the plane of interest
    # clipZ : length-2 list of ints
    # clipX : length-2 list of ints
    # clipY : length-2 list of ints
    #
    os.makedirs(output_dir, exist_ok=True)
    max_proj_dir = os.path.join(output_dir, 'mips')
    os.makedirs(max_proj_dir, exist_ok=True)
    midz_output_dir = os.path.join(output_dir, 'midZ')
    os.makedirs(midz_output_dir, exist_ok=True)
    if do_midX:
        if len(planeX) == 1 and planeX[0] < 0:
            midx_output_dir = os.path.join(output_dir, 'midX')
            os.makedirs(midx_output_dir, exist_ok=True)
        else:
            for i in range(len(planeX)):
                midx_output_dir = os.path.join(output_dir, f'midX{planeX[i]:04d}')
                os.makedirs(midx_output_dir, exist_ok=True)

    for filename in sorted(os.listdir(input_dir)):
        if filename.endswith('.h5'):
            base_filename = os.path.splitext(filename)[0]

            # MIPs
            os.makedirs(os.path.join(max_proj_dir, dir_name+"_mips"), exist_ok=True)
            max_proj_path = os.path.join(max_proj_dir, dir_name+"_mips", f"{base_filename}_mip_axis_x.tiff")

            # Z saggital
            os.makedirs(os.path.join(midz_output_dir, dir_name+"_midZ"), exist_ok=True)
            midz_path = os.path.join(midz_output_dir, dir_name+"_midZ", f"{base_filename}_midZ_axis_x.tiff")

            # X cranial
            write_midX = False
            if do_midX:
                for xplane in planeX:
                    midx_output_dir = os.path.join(output_dir, f'midX{xplane:04d}')
                    os.makedirs(os.path.join(midx_output_dir, dir_name+f"_midX{xplane:04d}"), exist_ok=True)
                    midx_path = os.path.join(midx_output_dir, dir_name+f"_midX{xplane:04d}", f"{base_filename}_midX{xplane:04d}.tiff")
                    write_midX = write_midX or not os.path.exists(midx_path)

            if not os.path.exists(max_proj_path) or not os.path.exists(midz_path) or write_midX:
                filepath = os.path.join(input_dir, filename)
                print('-> ' + filepath)
                with h5py.File(filepath, 'r') as f:
                    data = f['/Data'][()]

                print(np.shape(data))

                if do_mips:
                    print(f"--> {max_proj_path}")

                    mip_x = np.max(data[clipZ[0]:data.shape[1]-clipZ[1],
                                   clipX[0]:data.shape[1] - clipX[1],
                                   clipY[0]:data.shape[2]-clipY[1]], axis=0)
                    # mip_y = np.max(data, axis=1)
                    # mip_z = np.max(data, axis=2)
                    Image.fromarray(mip_x).save(max_proj_path)
                    # save_mip_images(mip_y, 'y', base_filename, output_dir)
                    # save_mip_images(mip_z, 'z', base_filename, output_dir)

                if do_midZ:
                    print(f"--> {midz_path}")
                    # Subset of z-stack as mip (mid-saggital, for ex)
                    if plane < 0:
                        # take middle z plane
                        plane = np.round(np.size(data, 0) * 0.5).astype(np.int64)

                    midz = data[plane, clipX[0]:data.shape[1] - clipX[1],
                                   clipY[0]:data.shape[2]-clipY[1]]
                    Image.fromarray(midz).save(midz_path)

                # Subset of x-stack as mip (mid-cranial)
                if write_midX:
                    for xplane in planeX:
                        if xplane < 0:
                            # take middle z plane
                            xplane = np.round(np.size(data, 1) * 0.5).astype(np.int64)
                            midx_path = os.path.join(midx_output_dir, dir_name + "_midX", f"{base_filename}_midX.tiff")
                        else:
                            midx_output_dir = os.path.join(output_dir, f'midX{xplane:04d}')
                            midx_path = os.path.join(midx_output_dir, dir_name + f"_midX{xplane:04d}",
                                                     f"{base_filename}_midX{xplane:04d}.tiff")

                        print(f"--> {midx_path}")
                        midx = data[:, xplane, clipY[0]:data.shape[2]-clipY[1]]
                        Image.fromarray(midx).save(midx_path)

            else:
                print(f"--> Output file already exists, skipping: {max_proj_path}")


if __name__ == "__main__":
    parent_dir = 'E:\\avistrok\\bapGAL4_UAShidUASstingerHiRFP\\2025-08-06_121230' #bapGAL4_UAShidUASstingerHiRFP'
    parent_dir = '/project/npmitchell/avistrok/bapGAL4_UAShidUASstingerHiRFP/2025-08-04_131451/'

    #parent_dir = 'E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-17\\2025-06-17_165029'
    # "E:\boris\bynGAL4klar_UASmChCAAXHiRFP\2025-06-16\2025-06-16_153503"
    multiple_dirs = False
    #list_of_parent_dirs = ['E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-28\\2025-06-28_173043', 'E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-28\\2025-06-28_173646', 'E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-28\\2025-06-28_173852', 'E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-28\\2025-06-28_181654', 'E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-28\\2025-06-28_181857', 'E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-28\\2025-06-28_184516', 'E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-28\\2025-06-28_191308', 'E:\\boris\\bynGAL4klar_UASmChCAAXHiRFP\\2025-06-28\\2025-06-28_202253']

    do_midX = True
    do_midZ = True
    do_mips = True
    planeX = [620, 720, 820, 920, 1020]

    clipZ = [30,30]
    clipY = [00,00]
    clips = dict['z': clipZ, 'y': clipY]

    if multiple_dirs is False:
        data_dir = os.path.join(parent_dir, 'raw')
        for dir_name in os.listdir(data_dir):
            input_dir = os.path.join(data_dir, dir_name)
            print(input_dir)
            if os.path.isdir(input_dir):
                output_dir = parent_dir  # os.path.join(outdir, 'mips', dir_name)
                print('processing ' + input_dir + " -> " + output_dir + "\\mips or midZ")
                process_directory(input_dir, output_dir, dir_name, clipZ=clipZ, clipY=clipY, do_mips=do_mips,
                                  do_midZ=do_midZ, do_midX=do_midX, planeX=planeX)
                print('done with ' + input_dir)
    else:
        for parent_dir in list_of_parent_dirs:
            data_dir = os.path.join(parent_dir, 'raw')
            for dir_name in os.listdir(data_dir):
                input_dir = os.path.join(data_dir, dir_name)
                print(input_dir)
                if os.path.isdir(input_dir):
                    output_dir = parent_dir  # os.path.join(outdir, 'mips', dir_name)
                    print('processing ' + input_dir + " -> " + output_dir + "\\mips or midZ")
                    process_directory(input_dir, output_dir, dir_name, clipZ=clipZ, clipY=clipY)
                    print('done with ' + input_dir)




    print('done')
    print('Now drag all mips subdirectories into Fiji. Run maxProjectMovies_simple.ijm after changing path in line 5')
