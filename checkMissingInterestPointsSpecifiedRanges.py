import os

import os


def check_missing_files(directory, tp_ids, view_setup_ids):
    """
    Checks for missing files in the specified directory based on expected naming conventions.

    Expected files follow the format:
        tpId_<tp_id>_viewSetupId_<view_setup_id>.<suffix>
    where:
        - tp_id is either 9 or 48
        - view_setup_id ranges from 1 to 17
        - suffix is either 'beads.corr.txt' or 'beads.ip.txt'

    Prints a list of missing files if any are not found.

    Parameters
    ----------
    directory: path-like str
        Path to the directory containing the files
    tp_ids : integer list
        timepoints to look for
    view_setup_ids : list of integers
        view_setup_ids to look for
    """
    if not os.path.isdir(directory):
        print(f"Error: Directory '{directory}' does not exist.")
        return

    suffixes = ["beads.corr.txt", "beads.ip.txt"]
    missing_files = []

    for tp_id in tp_ids:
        for view_setup_id in view_setup_ids:
            for suffix in suffixes:
                filename = f"tpId_{tp_id}_viewSetupId_{view_setup_id}.{suffix}"
                if not os.path.isfile(os.path.join(directory, filename)):
                    missing_files.append(filename)

    if missing_files:
        print("Missing files:")
        for file in missing_files:
            print(file)
    else:
        print("All expected files are present.")



if __name__ == "__main__":
    """
    Defines the directory to check and executes the file-checking function.
    """
    directory = "/project/npmitchell/canto/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/2024-06-26_175116_crisp_120s_HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/unpacked/interestpoints/"  # Define the directory here
    # directory = "/project/npmitchell/alex/HandGFP48YGAL4UASmChCAAXHiRFP_20250303_22p1C_excellent/2025-03-03_162051/unpacked/interestpoints"

    tp_ids = range(0,119)
    view_setup_ids = range(0, 18)
    check_missing_files(directory, tp_ids, view_setup_ids)