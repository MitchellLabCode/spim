
datdir = "/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/data_resavedViaMATLAB/fused2x/";
tps = newArray(0,10,20,30,40,50,60,70,80,90,100);


//loop over the various camera view indices
for (i = 0; i < tps.length; i++){
	tp = tps[i] ;

	// ch1
	open(datdir + "TP" + tp + "_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	selectWindow("TP" + tp + "_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir + "mips/c1_TP" + tp + ".tif");
	close();

	// reslice laterally
	selectWindow("TP" + tp + "_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	run("Reslice [/]...", "output=2.000 start=Left");
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir + "mips_left/c1_TP" + tp + ".tif");
	close();
	close() ;
	
	// ch2
	open(datdir + "TP" + tp + "_Ch2_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	selectWindow("TP" + tp + "_Ch2_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir + "mips/c2_TP" + tp + ".tif");
	close();
	
	// reslice laterally
	selectWindow("TP" + tp + "_Ch2_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	run("Reslice [/]...", "output=2.000 start=Left");
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir + "mips_left/c2_TP" + tp + ".tif");
	close();
	close();
	close() ;
	
	// Merge
	//run("Merge Channels...", "c2=[MAX_TP" + tp + "_Ch2_Ill0_Ang0,45,90,135,180,225,270,315.tif] c6=[MAX_TP" + tp + "_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif] create");
	//saveAs("Tiff", "/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/data_resavedViaMATLAB/fused2x/mips/c2_TP" + tp + ".tif");
	//close();
	
}

run("Image Sequence...", "dir=" + datdir + "mips/ filter=c1 sort");
saveAs("Tiff", datdir + "mips_c1.tif");
run("Image Sequence...", "dir=" + datdir + "mips/ filter=c2 sort");
saveAs("Tiff", datdir + "mips_c2.tif");
run("Merge Channels...", "c2=[mips_c2.tif] c6=[mips_c1.tif] create");
// Labels and scalebar
run("Label...", "format=00:00 starting=0 interval=3 x=5 y=20 font=24 text=[] range=1-11 use");
run("Set Scale...", "distance=1 known=0.5238 unit=um");
run("Scale Bar...", "width=50 height=50 thickness=6 font=24 color=White background=None location=[Lower Right] horizontal overlay");
saveAs("Tiff", datdir + "mips_Composite.tif");


run("Image Sequence...", "dir=" + datdir + "mips_left/ filter=c1 sort");
saveAs("Tiff", datdir + "mips_left_c1.tif");
run("Image Sequence...", "dir=" + datdir + "mips_left/ filter=c2 sort");
saveAs("Tiff", datdir + "mips_left_c2.tif");
run("Merge Channels...", "c2=[mips_left_c2.tif] c6=[mips_left_c1.tif] create");
// Labels and scalebar
run("Label...", "format=00:00 starting=0 interval=3 x=5 y=20 font=24 text=[] range=1-11 use");
run("Set Scale...", "distance=1 known=0.5238 unit=um");
run("Scale Bar...", "width=50 height=50 thickness=6 font=24 color=White background=None location=[Lower Right] horizontal overlay");
saveAs("Tiff", datdir + "mips_left_Composite.tif");


