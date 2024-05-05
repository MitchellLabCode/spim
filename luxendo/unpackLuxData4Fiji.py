import os
import h5py
import numpy as np
from tifffile import imwrite

"""Unpack Luxendo datasets to OME-TIFF format readable by Preibisch plugins without limiting the LUT"""


def convert_hdf5_files_in_directory(directory, overwrite):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".h5") or file.endswith(".hdf5"):
                file_path = os.path.join(root, file)
                convert_hdf5_data_to_ome_tiff(file_path, overwrite)


def convert_hdf5_data_to_ome_tiff(filename, overwrite):
    with h5py.File(filename, 'r') as f:
        for dataset_name in f:
            dataset = f[dataset_name]
            if is_unsigned_integer(dataset):
                uint16_data = dataset[:]
                save_as_ome_tiff(filename, dataset_name, uint16_data, overwrite)
+_

def is_unsigned_integer(dataset):
    dtype = dataset.dtype
    return np.issubdtype(dtype, np.uint16)


def save_as_ome_tiff(hdf5_filename, dataset_name, uint16_data, overwrite):
    tiff_filename = f"{os.path.splitext(hdf5_filename)[0]}_{dataset_name}.ome.tif"
    # Convert uint16 data to uint8 for OME-TIFF (assuming data is in the range [0, 65535])
    # uint8_data = (uint16_data / 65535 * 255).astype(np.uint8)
    imwrite(tiff_filename, uint16_data,
            metadata={'axes': 'ZYX'})  # You may need to adjust metadata as per your requirement
    print(f"Converted {dataset_name} in {hdf5_filename} to {tiff_filename}")


if __name__ == "__main__":
    directory_path = "D:\\npmitchell\\202403271400_james_3color_16bit\\2024-03-27_135655\\raw"
    overwrite = True
    convert_hdf5_files_in_directory(directory_path, overwrite)
    print("Conversion completed.")
