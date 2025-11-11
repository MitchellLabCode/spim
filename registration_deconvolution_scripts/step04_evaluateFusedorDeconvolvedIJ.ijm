

tp_array = newArray(0,1,50,90,130,170,210,226);
for (i = 0; i <= lengthOf(tp_array) ; i++){
	run("Merge Channels...", "c2=TP" + tp_array[i] +"_Ch0_Ill0_Ang48,108,168,228,288,348.tif c6=TP" + tp_array[i] + "_Ch1_Ill0_Ang48,108,168,228,288,348.tif create");
}

