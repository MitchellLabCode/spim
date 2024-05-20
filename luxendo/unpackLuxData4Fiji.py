import os
import h5py
import numpy as np
from tifffile import imwrite
import xml.etree.ElementTree as ET
import re

"""Unpack Luxendo datasets to OME-TIFF format readable by Preibisch plugins without limiting the LUT"""


def parse_bdv_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()

    angle_map = {}

    for setup in root.findall(".//ViewSetup"):
        setup_id = setup.find("id").text
        setup_name = setup.find("name").text

        # Extract angle value using regular expression
        match = re.search(r'ang:h(\d+)', setup_name)
        if match:
            angle_value = match.group(1)
            if "left" in setup_name:
                angle_map[f"stack_{setup_id}_left"] = f"Angle{angle_value}"
            elif "right" in setup_name:
                angle_map[f"stack_{setup_id}_right"] = f"Angle{angle_value}"

    return angle_map


def convert_hdf5_files_in_directory(directory, output_directory, bdv_xml_file, overwrite):
    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok=True)

    # Parse the bdv.xml file
    angle_map = parse_bdv_xml(bdv_xml_file)

    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".h5") or file.endswith(".hdf5"):
                file_path = os.path.join(root, file)
                convert_hdf5_data_to_ome_tiff(file_path, output_directory, angle_map, overwrite)


def convert_hdf5_data_to_ome_tiff(filename, output_directory, angle_map, overwrite):
    with h5py.File(filename, 'r') as f:
        for dataset_name in f:
            dataset = f[dataset_name]
            if is_unsigned_integer(dataset):
                uint16_data = dataset[:]
                save_as_ome_tiff(filename, dataset_name, uint16_data, output_directory, angle_map, overwrite)


def is_unsigned_integer(dataset):
    dtype = dataset.dtype
    return np.issubdtype(dtype, np.uint16)


def save_as_ome_tiff(hdf5_filename, dataset_name, uint16_data, output_directory, angle_map, overwrite):
    subdirectory_name = os.path.basename(os.path.dirname(hdf5_filename))
    original_filename = os.path.splitext(os.path.basename(hdf5_filename))[0]

    # Check if the dataset_name matches any key in the angle_map
    angle_suffix = angle_map.get(dataset_name, dataset_name)
    tiff_filename = f"{subdirectory_name}_{original_filename}_{angle_suffix}.ome.tif"
    tiff_filename = tiff_filename.replace('_left_Cam_left', '_')
    tiff_filename = tiff_filename.replace('_right_Cam_right', '_')
    tiff_filename = tiff_filename.replace('.lux.', '.')

    tiff_path = os.path.join(output_directory, tiff_filename)

    # Convert uint16 data to uint8 for OME-TIFF (assuming data is in the range [0, 65535])
    # uint8_data = (uint16_data / 65535 * 255).astype(np.uint8)
    imwrite(tiff_path, uint16_data,
            metadata={'axes': 'ZYX'})  # You may need to adjust metadata as per your requirement
    print(f"Converted {dataset_name} in {hdf5_filename} to {tiff_path}")


if __name__ == "__main__":
    directory_path = 'D:\\rlondo\\bynGAL4klar_UAS-mCherry-CAAX\\2024-04-25_155740\\raw'  # Replace with your parent directory path
    output_directory = 'D:\\output_directory'  # Replace with your output directory path
    bdv_xml_file = os.path.join(directory_path, "bdv.xml")  # Path to the bdv.xml file
    overwrite = True
    convert_hdf5_files_in_directory(directory_path, output_directory, bdv_xml_file, overwrite)
    print("Conversion completed.")
