import os
import h5py
import numpy as np
from tifffile import imwrite
import xml.etree.ElementTree as ET # no need to install
import re # no need to install#
from ome_types import from_xml
from ome_types.model import OME, Image, Pixels, Channel, TiffData

"""Unpack Luxendo datasets to OME-TIFF format readable by Preibisch plugins without limiting the LUT"""
def parse_bdv_xml(xml_file): # xml_file is bdv_xml_file, which is the path to the 'bdv.xml'
    # 'parse' means make a file in one format to another format, so that the file structure is clearer
    # xml_file = bdv_xml_file
    tree = ET.parse(xml_file) # <xml.etree.ElementTree.ElementTree object at 0x7dfd8bc54790>
    root = tree.getroot() # <Element 'SpimData' at 0x7dfd8b341300>
    # it seems that this memory address can change every time the code s executed
    # root.tag # 'SpimData'
    # root.attib # {'version': '0.2'}s

    angle_map = {}

    for setup in root.findall(".//ViewSetup"):
        # len(root.findall(".//ViewSetup")) # this is a list. Len = 12.
        # In every <ViewSetup>, there are <id>, <name>, etc
        # root.findall(".//ViewSetup")[1] # <Element 'ViewSetup' at 0x7dfd8b3410d0>
        setup_id = setup.find("id").text  # 1, 2, 3, ...
        setup_name = setup.find("name").text  # for example
        stack = setup_name[8]
        # print('This is stack ' + stack)
        # 'ch:0_st:0_ang:h45-v90_obj:right_cam:right'
        if 'left' in setup_name:
            start_idx = setup_name.find('ang:h')
            end_idx = setup_name.find('-v')
            angle_map['stack_' + setup_name[8] + '_channel_' + setup_name[3] + '_obj_left'] = setup_name[start_idx + 5:end_idx]

        if 'right' in setup_name:
            start_idx = setup_name.find('ang:h')
            end_idx = setup_name.find('-v')
            angle_map['stack_' + setup_name[8] + '_channel_' + setup_name[3] + '_obj_right'] = setup_name[start_idx + 5:end_idx]

        # print(angle_map)
        # {'stack_1_channel_0_obj_left': '240', 'stack_1_channel_0_obj_right': '60', 'stack_1_channel_2_obj_left': '240', 'stack_1_channel_2_obj_right': '60', 'stack_2_channel_0_obj_left': '300', 'stack_2_channel_0_obj_right': '120', 'stack_2_channel_2_obj_left': '300', 'stack_2_channel_2_obj_right': '120', 'stack_3_channel_0_obj_left': '0', 'stack_3_channel_0_obj_right': '180'}
    return angle_map

def convert_hdf5_files_in_directory(directory, output_directory, bdv_xml_file, overwrite):
    # Create output directory if it doesn't exist
    # directory = datadir
    os.makedirs(output_directory, exist_ok=True)

    # Parse the bdv.xml file
    angle_map = parse_bdv_xml(bdv_xml_file)

    for root, dirs, files in os.walk(directory):
        # os.walk("...") returns a generator   <generator object _walk at 0x7dfd8bc536f0>
        for file in sorted(files, key=extract_timepoint):
            if file.endswith(".h5") or file.endswith(".hdf5"):
                file_path = os.path.join(root, file)
                # file path (an example): '/mnt/data/midgut_tubulogenesis/48Y-GAL4-klar_UASmChCAAXSLamGFP/2024-05-06_143154/raw/stack_1_channel_1_obj_right/Cam_right_00077.lux.h5'
                # Use this file to test the function "convert_hdf5_data_to_ome_tiff"
                convert_hdf5_data_to_ome_tiff(file_path, output_directory, angle_map, overwrite)


def extract_timepoint(filename):
    # Adjust the regex to match the specific format of your filenames
    match = re.search(r'_(\d+)\.lux\.h5$', filename)
    return int(match.group(1)) if match else float('inf')


# Read dataset. We should convert the dataset to an inverted way,
# and also rotate the dataset after this step.
def convert_hdf5_data_to_ome_tiff(filename, output_directory, angle_map, overwrite):
    # filename = file_path
    # check if the output file already exists

    # hdf5_filename = filename
    subdirectory_name = os.path.basename(os.path.dirname(filename)) # stack_1_channel_1_obj_right
    original_filename = os.path.splitext(os.path.basename(filename))[0] # Cam_right_00077.lux.h5

    # Check if the dataset_name matches any key in the angle_map
    # angle_suffix = angle_map.get(dataset_name, dataset_name)
    for key in angle_map:
        if key in filename:
            # print('key is in filename')
            angle_suffix = angle_map[key]

    # angles=list(angle_map.values())
    # for i in angles:
    tiff_filename = f"{subdirectory_name}_{original_filename}_{angle_suffix}.ome.tif"
    tiff_filename = tiff_filename.replace('_left_Cam_left', '_tp')
    tiff_filename = tiff_filename.replace('_right_Cam_right', '_tp')
    tiff_filename = tiff_filename.replace('.lux', '_angle')
    tiff_filename = tiff_filename.replace('_obj', '')

    tiff_path = os.path.join(output_directory, tiff_filename)
    print('creating: ' + tiff_filename)

    # check if the file already exists on disk (tiff_path is full path)
    if os.path.isfile(tiff_path):
        print('file already on disk ' + tiff_filename)
    else:
        with h5py.File(filename, 'r') as f: # f is the name of the h5 file that is read by this code.
            # f acts like a dictionary in python
            # filename an example: '/mnt/data/.../raw/stack_1_channel_1_obj_right/Cam_right_00077.lux.h5'
            for dataset_name in f:
                dataset = f[dataset_name]
                # dataset_shape = dataset.shape
                if is_unsigned_integer(dataset):
                    uint16_data = dataset[:] # [::-1, :, :] does not work for hdf5 dataset
                    if 'left' in filename:
                        uint16_data_z_flipped = np.flip(uint16_data, axis=0) # Actually only the left one should be flipped.
                    elif 'right' in filename:
                        uint16_data_z_flipped = uint16_data
                        # print('From the right camera, not flipping')

                    # So later the right one should be done without this line.
                    # Or actually use if to judge if it is left / right
                    # we should invert the Z axis. The Z axis is the dim with shape == 266, which is the 0th axis here
                    # uint16_data_z_flipped[:, 0:1, 0:1]    AND ALSO   uint16_data[:, 0:1, 0:1]
                    # The above 2 do look inverted, and they do each has ~200 items. It seems that the z axis IS flipped
                    # save_as_ome_tiff(filename, dataset_name, uint16_data_z_flipped, output_directory, angle_map, overwrite)

                    # Create OME metadata
                    # UserWarning: Unrecognized fields for type <class 'ome_types._autogenerated.ome_2016_06.pixels.Pixels'>: {'tiff_data'}
                    #   coro = func()
                    ome_metadata = OME(  ## This is an OME object.
                        images=[
                            Image(  ## This OME object only has a single Image object, whose id is "image:0"
                                id="Image:0",
                                name=tiff_filename,
                                pixels=Pixels(
                                    dimension_order="XYZCT",  # orignially zyx, now converted to x,y,z,channel,time
                                    type="uint16",
                                    size_x=uint16_data_z_flipped.shape[2],
                                    size_y=uint16_data_z_flipped.shape[1],
                                    size_z=uint16_data_z_flipped.shape[0],
                                    size_c=1,
                                    size_t=1,
                                    channels=[Channel(id="Channel:0:0", samples_per_pixel=1)],
                                    # tiff_data=[
                                    #     TiffData(first_t=0, first_c=0, first_z=0, ifd=0)
                                    # ], # seems that the OME object does not have this attribute
                                    tiff_data_blocks=[
                                        TiffData(first_t=0, first_c=0, first_z=0, ifd=0)
                                    ]
                                )
                            )
                        ]
                    )

                    # Write the OME-TIFF file  - SHOULD BE MODIFIED
                    # imwrite(tiff_path, uint16_data_z_flipped, metadata={"axes": "ZYX", "OME": ome_metadata.to_xml()})
                    # NPM: This said XYZ until 2024-06-11, but we changed to XYZ above so flipping here in metadata. Is this right?
                    imwrite(tiff_path, uint16_data_z_flipped, metadata={'axes': 'XYZ'}, description=ome_metadata.to_xml())
                    print(f"Converted {dataset_name} in {filename} to {tiff_path}")
        """
        # with h5py.File(filename, 'r') as f:
            # f_keys = list(f.keys()) # ['Data', 'metadata']   # There are 2 datasets. (Also shown in Fiji when importing h5)
            # path: /Data       size: 283x2048x2048 type: uInt16    element size[um]: 1.0x0.195x0.195
            # path: /metadata   size: 1     type: uSTRING(-1):{}    element size[um]: unknown
            # dataset = f['Data']
            # dset_shape = dset.shape # {tuple:3}  (266,2048,2048)
            # dset_type = dset.dtype # dtype('uint16')
            # info on hdf5 object dataset: https://docs.h5py.org/en/stable/high/dataset.html
        """


def is_unsigned_integer(dataset):
    dtype = dataset.dtype
    return np.issubdtype(dtype, np.uint16)


# def save_as_ome_tiff(hdf5_filename, dataset_name, uint16_data_z_flipped, output_directory, angle_map, overwrite):


if __name__ == "__main__":

    directory_path = '/project/npmitchell/canto/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/2024-06-26_175116_crisp_120s_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP'
    directory_path = 'F:\\PROJECTS\\LightMicroscopyBootcamp2025\\PSF_33x_ALP_grp\\2025-10-14_121043'
    directory_path = 'F:\\PROJECTS\\LightMicroscopyBootcamp2025\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251014144825_60spf_combined\\2025-10-14_144825_tp0'

    datadir = os.path.join(directory_path, 'raw')
    # This is where the raw hdf5 files are. In linux use '/'
    output_directory = os.path.join(directory_path, 'unpacked')
    bdv_xml_file = os.path.join(directory_path, "bdv.xml")  # Path to the bdv.xml file
    overwrite = True
    convert_hdf5_files_in_directory(datadir, output_directory, bdv_xml_file, overwrite)
    print("Conversion completed.")
