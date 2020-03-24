# spimCode
code for SPIM microscopy

extractStacks...m files are MATLAB scripts for unpacking data from "as fast as possible" format from lightsheet computer into separate files for each view, timepoint, & channel. 

hash_md5_parsing is a collection of codes for checking the output of an md5deep checksum against another instance of the same checksum. For ex, say you have a dataset on one computer and transferred it to another computer. Call these two copies now data1 and data2. Are they actually identical? Hard to say a priori. What to do? Run a checksum on each (separately, one checksum on one computer, the other on the other computer). Now we wish to compare the long txt files created by each to compare hashcodes. To do so, first change the path to the file (since these are different machines, they have different paths). So for example, hashlist2.txt is the same (scrambled) list as hashlist1.txt in the examples directory, but the path of the data is E:/Runt/ instead of D:/Atlas_Data/Runt/, so the parser would think these are different files. To prep hashlist2.txt for comparison, you can run
$ cd /path/to/this/code/hash_md5_parsing/
$ python replace_substring_in_hashfile.py examples/hashlist2.txt examples/hashlist2_fixed.txt "E:\Runt" "D:\Atlas_Data\Runt"
Then compare the string-replaced hash with the original:
$ python compare_hashes.py examples/hashlist1.txt examples/hashlist_tweak.txt examples/hash_comparison.txt
which has the syntax:
$ python program_name.py <hashlist_file> <master_hashlist_file> <output_filename>

