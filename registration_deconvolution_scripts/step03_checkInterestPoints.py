import os
import re
import matplotlib.pyplot as plt
import argparse


def extract_tp_and_view_setup_ids(directory, pattern):
    """
    Extracts tp_ids and view_setup_ids from filenames in the given directory.

    Parameters
    ----------
    directory : str
        Path to the directory containing the files.

    Returns
    -------
    tuple
        A sorted list of tp_ids and a dictionary mapping tp_ids to sets of view_setup_ids.
    """
    tp_ids = set()
    view_setup_ids = {}

    for filename in os.listdir(directory):
        match = pattern.match(filename)
        if match:
            channel, tp_id, angle = map(int, match.groups())
            tp_ids.add(tp_id)
            view_setup_id = channel * 6 + angle // 60
            view_setup_ids.setdefault(tp_id, set()).add(view_setup_id)

    return sorted(tp_ids), view_setup_ids


def check_missing_files(directory, tp_ids, view_setup_ids):
    """
    Checks for missing files in the interestpoints subdirectory based on extracted tp_ids and view_setup_ids.
    Additionally, counts the number of lines in each found "beads.corr.txt" file and plots the results.

    Parameters
    ----------
    directory : str
        Path to the main directory containing the 'interestpoints' subdirectory.
    tp_ids : list
        List of extracted time point IDs.
    view_setup_ids : dict
        Dictionary mapping tp_ids to sets of view_setup_ids.
    """
    interestpoints_dir = os.path.join(directory, 'interestpoints')
    if not os.path.isdir(interestpoints_dir):
        print(f"Error: Directory '{interestpoints_dir}' does not exist.")
        return

    suffixes = ["beads.corr.txt", "beads.ip.txt"]
    missing_files = []
    line_counts = {}

    for tp_id in tp_ids:
        for view_setup_id in view_setup_ids.get(tp_id, []):
            for suffix in suffixes:
                filename = f"tpId_{tp_id}_viewSetupId_{view_setup_id}.{suffix}"
                filepath = os.path.join(interestpoints_dir, filename)
                if not os.path.isfile(filepath):
                    missing_files.append(filename)
                else:
                    if suffix == "beads.ip.txt":
                        with open(filepath, 'r') as file:
                            line_counts.setdefault(view_setup_id, []).append((tp_id, sum(1 for _ in file)-1))

    if missing_files:
        print("Missing files:")
        for file in missing_files:
            print(file)
    else:
        print("All expected files are present.")

    plt.figure(figsize=(10, 5))
    for view_setup_id, data in line_counts.items():
        data.sort()
        x_vals = [tp for tp, _ in data]
        y_vals = [count for _, count in data]
        plt.plot(x_vals, y_vals, label=f"View Setup {view_setup_id}")

    plt.xlabel("Timepoint")
    plt.ylabel("Interest points")
    plt.title("Interest points summary")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    plt.tight_layout()
    plot_path = os.path.join(directory, 'interestpoints_summary.png')
    plt.savefig(plot_path)
    print(f"Plot saved to {plot_path}")
    plt.show()


if __name__ == "__main__":
    """
    Main function to define directory and execute file checking and plotting.
    
    Example Usage
    -------------
    
    """
    parser = argparse.ArgumentParser(description="Check for missing interest points and generate summary plot.")
    parser.add_argument(
        "-dir", "--directory",
        type=str,
        default="./",
        help="Path to the directory containing the data files."
    )

    args = parser.parse_args()
    directory = args.directory

    pattern = re.compile(r'c(\d+)_t(\d+)_a(\d+)\.ome\.tif')
    tp_ids, view_setup_ids = extract_tp_and_view_setup_ids(directory, pattern)
    check_missing_files(directory, tp_ids, view_setup_ids)