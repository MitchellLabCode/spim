import os
import h5py
import numpy as np
from tifffile import imwrite

"""Unpack Luxendo datasets to OME-TIFF format readable by Preibisch plugins without limiting the LUT"""


def convert_hdf5_files_in_directory(directory, output_directory, overwrite):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".h5") or file.endswith(".hdf5"):
                file_path = os.path.join(root, file)
                convert_hdf5_data_to_ome_tiff(file_path, output_directory, overwrite)


def convert_hdf5_data_to_ome_tiff(filename, output_directory, overwrite):
    with h5py.File(filename, 'r') as f:
        for dataset_name in f:
            dataset = f[dataset_name]
            if is_unsigned_integer(dataset):
                uint16_data = dataset[:]
                save_as_ome_tiff(filename, dataset_name, uint16_data, output_directory, overwrite)


def is_unsigned_integer(dataset):
    dtype = dataset.dtype
    return np.issubdtype(dtype, np.uint16)


def save_as_ome_tiff(hdf5_filename, dataset_name, uint16_data, output_directory, overwrite):
    subdirectory_name = os.path.basename(os.path.dirname(hdf5_filename))
    original_filename = os.path.splitext(os.path.basename(hdf5_filename))[0]
    tiff_filename = f"{subdirectory_name}_{original_filename}.ome.tif"
    tiff_path = os.path.join(output_directory, tiff_filename)

    # Convert uint16 data to uint8 for OME-TIFF (assuming data is in the range [0, 65535])
    # uint8_data = (uint16_data / 65535 * 255).astype(np.uint8)
    imwrite(tiff_path, uint16_data,
            metadata={'axes': 'ZYX'})  # You may need to adjust metadata as per your requirement
    print(f"Converted {dataset_name} in {hdf5_filename} to {tiff_path}")

if __name__ == "__main__":
    directory_path = 'D:\\rlondo\\bynGAL4klar_UAS-mCherry-CAAX\\2024-04-25_155740\\raw'  # Replace with your parent directory path
    output_directory = 'D:\\rlondo\\bynGAL4klar_UAS-mCherry-CAAX\\2024-04-25_155740\\tif'  # Replace with your output directory path
    overwrite = True

    os.makedirs(output_directory, exist_ok=True)
    convert_hdf5_files_in_directory(directory_path, output_directory, overwrite)
    print("Conversion completed.")
