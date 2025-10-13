import os
import glob

# Directory and search pattern
datadir = '/path/to/data/'
searchstr = 'data*Time_0*_Angle_*.ome.tif'
removestr = 'data\\'   # note: backslash must be escaped
replacementstr = ''

# Get matching files
files = glob.glob(os.path.join(datadir, searchstr))

# Iterate through files and rename
for oldfn in files:
    pathn, filename = os.path.split(oldfn)
    newfilename = filename.replace(removestr, replacementstr)
    newfn = os.path.join(pathn, newfilename)
    
    print(f"Replacing {filename}\n  with -> {newfilename}")
    os.rename(oldfn, newfn)
