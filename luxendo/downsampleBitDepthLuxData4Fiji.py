import os
import h5py
import numpy as np

"""Convert Luxendo datasets to format readable by Preibisch plugins by limiting LUT (for signed 16bit max).
This reduces the information in the data by a factor of two.
"""

def convert_hdf5_files_in_directory(directory, forceOverwrite):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".h5") or file.endswith(".hdf5"):
                file_path = os.path.join(root, file)

                # Check if converted.txt exists
                txt_file_path = os.path.join(root, "converted.txt")
                if not os.path.exists(txt_file_path) or forceOverwrite:
                    print('Considering ' + file_path)
                    convert_hdf5_data(file_path)

                    # mark as converted
                    with open(txt_file_path, "w") as txt_file:
                        txt_file.write("converted")
                else:
                    print('Already done: ' + file_path)


def convert_hdf5_data(filename):
    with h5py.File(filename, 'r+') as f:
        for dataset_name in f:
            dataset = f[dataset_name]

            # Check if dataset is unsigned integer
            if is_unsigned_integer(dataset):
                print('>>Converting ' + filename)
                # Convert unsigned integer to signed integer
                minval = np.min(dataset[:])
                maxval = np.max(dataset[:])
                signed_data = convert_data(dataset)
                print('>> max: ' + str(maxval) + ' -> ' + str(np.max(signed_data[:])))
                print('>> min: ' + str(minval) + ' -> ' + str(np.min(signed_data[:])))
                # assert (np.max(signed_data[:]) <= np.iinfo(np.int16).max)
                # Write back to dataset
                dataset[...] = signed_data


def is_unsigned_integer(dataset):
    dtype = dataset.dtype
    return np.issubdtype(dtype, np.uint16)


def convert_data(dataset):
    # divide the values by 2 to keep ints in range [-2^15, 2^15].
    unsigned_data = np.floor(dataset[:] * 0.999)
    signed_data = np.array(unsigned_data, dtype=np.uint32)
    return signed_data


if __name__ == "__main__":
    directory_path = "D:\\npmitchell\\202403271400_james_3color\\2024-03-27_135655\\raw"
    convert_hdf5_files_in_directory(directory_path, True)
    print("Conversion completed.")
