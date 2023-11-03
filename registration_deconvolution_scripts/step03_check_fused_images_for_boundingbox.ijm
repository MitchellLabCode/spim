// datdir="/mnt/crunch/48YGAL4klarLamGFPCAAXmCh/202308122137_180s_1p4um_2mW2mW_48YG4kLamGFPCAAXmCh_4views_0p25_0p5ms/dataLamGFP_resavedViaMATLAB/fused3x_LamGFP/";
datdir="/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/data_resavedViaMATLAB/";

	
// Define the range of angles you want to iterate over
tps = newArray(0, 1, 10, 30, 45, 65, 75, 90, 101);

for (i = 0; i < tps.length; i++) {
    tp = tps[i];
    
	open(datdir+"TP"+tp+"_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	open(datdir+"TP"+tp+"_Ch2_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	
	run("Merge Channels...", "c2=TP"+tp+"_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif c6=TP"+tp+"_Ch2_Ill0_Ang0,45,90,135,180,225,270,315.tif create");
	selectWindow("Composite");
	saveAs("Tiff", datdir+"fused3x_GFPnlsCAAXmCh/composite_TP"+tp+".tif");

	// GFPnls starting guess:
	makeRectangle(159, 246, 200, 326);
	// LamGFP starting guess: makeRectangle(213, 142, 205, 472);
	// close();


}