import os
import argparse

"""
Compare two hash files and check for different hashes or missing data files. Report missing files on command line and 
in output file. 
Written to compare output of md5deep-4.3 software.
Hash files are to be in the format:
hash1  /path/to/data/file1.exten
hash2  /path/to/data/file2.exten
...

Example Usage
-------------
python compare_hashes.py examples/hashlist1.txt examples/Atlas_Data_Runt_Mobile5drive.txt examples/hash_comparison.txt
python compare_hashes.py examples/hashlist1.txt examples/hashlist_tweak.txt examples/hash_comparison.txt


"""


def parse_hashlist(input_filename):
    """Load dictionary from txt file in outdir. This is equivalent to lattice_elasticity.load_params()

    Parameters
    ----------
    input_filename: str
        The path to the txt file with hashes and filenames.

    Returns
    -------
    hashd: dict
        A dictionary, with all hashes associated with filenames

    """
    hashd = {}
    # If outdir is the entire path, it has .txt at end, and so use this to split it up into dir and filename

    with open(input_filename) as f:
        # for line in f:
        #     print line
        for line in f:
            if '# ' not in line:
                (hash, fn) = line.split('  ')
                hashd[fn] = hash

    return hashd


def compare_hash_dicts(hd1, hd2):
    """Compare two hash dictionaries for similarity

    Parameters
    ----------
    hd1: dict
        hash dictionary 1
    hd2: dict
        hash dictionary 2

    Returns
    -------
    comp : dict
        dictionary with keys:
            filename : str, filename
            missing : dict with (keys, values)
                (filename, hash) : str, str
            extra : dict with (keys, values)
                (filename, hash) : str, str
            different : dict with (keys, values)
                (filename, [hash, hash_master]) : str, list of strings
            filename_master : str, filename
    """
    comp = {'missing': {}, 'extra': {}, 'different': {}}
    for key in hd1:
        if key in hd2:
            if hd1[key] != hd2[key]:
                comp['different'][key] = [hd1[key], hd2[key]]
        else:
            comp['extra'][key] = hd1[key]

    for key in hd2:
        if key not in hd1:
            comp['missing'][key] = hd2[key]

    return comp


def write_comparison_dict(dd, outfn, fn1, fn2, padding_var=7):
    """Write the comparison to a text file with easy-to-read formatting

    Parameters
    ----------
    dd : dict
        the comparison dictionary
    outfn : str
        file to write the comparison to disk

    Returns
    -------
    """
    with open(outfn, 'w') as myfile:
        myfile.write('# Comparison of ' + fn1 + ' with ' + fn2 + '\n')

    for mainkey in dd:
        with open(outfn, 'a') as myfile:
            myfile.write('\n\n# '+mainkey+'\n')

        subd = dd[mainkey]
        # print('mainkey = ' + mainkey)
        # print('subd = ')
        # print(subd)

        if len(subd) == 0:
            with open(outfn, 'a') as myfile:
                myfile.write('none\n')
        else:
            for key in subd:
                # print('key = ', key)
                # print('value = ', subd[key])
                with open(outfn, 'a') as myfile:
                    if isinstance(subd[key], str) or isinstance(subd[key], list):
                        myfile.write('{{0: <{}}}'.format(padding_var).format(key))
                        # + '= ' + subd[key] + '\n')
                    else:
                        raise RuntimeError('Could not write key,value pair to disk, since not a string.')


if __name__ == "__main__":
    # PARSE ARGUMENTS
    parser = argparse.ArgumentParser(description='Compare two hash lists and report missing files')
    # Build two positional arguments, one for each file
    parser.add_argument('file1', type=str, nargs='?',
                        help='Full or relative path to first hash file to compare',
                        default='check_string_for_empty')
    parser.add_argument('file2', type=str, nargs='?',
                        help='Full or relative path to second hash file to compare',
                        default='check_string_for_empty')
    # Build final optional positional argument for output file
    parser.add_argument('outfn', nargs='?', help='Path to output file with info on the comparison',
                        type=str, default='./hash_comparison.txt')
    args = parser.parse_args()

    if args.file1 == 'check_string_for_empty':
        raise RuntimeError('Must supply filepath for hash list to check')
    if args.file2 == 'check_string_for_empty':
        raise RuntimeError('Must supply filepath for master hash list to check against')

    print('Parsing file 1...')
    hashd1 = parse_hashlist(args.file1)
    print('Parsing file 2...')
    hashd2 = parse_hashlist(args.file2)
    print('Comparing hashes and filenames...')
    comp = compare_hash_dicts(hashd1, hashd2)
    # print(comp)
    print('Writing to disk')
    write_comparison_dict(comp, args.outfn, args.file1, args.file2, padding_var=7)
    print('Wrote comparison to file: ' + args.outfn)

