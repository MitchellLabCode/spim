
datdir = "/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/data_resavedViaMATLAB/fused2x/";
tps = newArray(0,10,20,30,40,50,60,70,80,90,100);
slice0 = "165-175" ;
slice1 = "265-275" ;
slice2 = "365-375" ;
rotAngle = "122" ;

//loop over the various camera view indices
for (i = 0; i < tps.length; i++){
	tp = tps[i] ;

	//ch1
	open(datdir + "TP" + tp + "_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	selectWindow("TP" + tp + "_Ch1_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	run("Reslice [/]...", "output=2.000 start=Top");

	// MIP
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir+"reslice/mips/MAX_reslice_c1_TP" + tp + ".tif");
	close();
	
	// Get MIPs of thin slices of the restack
	run("Make Substack...", "slices="+slice0);
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir+"reslice/slice0/MAX_slice0_c1_TP" + tp + ".tif");
	close();
	close();
	run("Make Substack...", "slices="+slice1);
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir+"reslice/slice1/MAX_slice1_c1_TP" + tp + ".tif");
	close();
	close();
	run("Make Substack...", "slices="+slice2);
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir+"reslice/slice2/MAX_slice2_c1_TP" + tp + ".tif");
	close();
	close();

	//ch2
	open(datdir + "TP" + tp + "_Ch2_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	selectWindow("TP" + tp + "_Ch2_Ill0_Ang0,45,90,135,180,225,270,315.tif");
	run("Reslice [/]...", "output=2.000 start=Top");
	
	// MIP
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir+"reslice/mips/MAX_reslice_c2_TP" + tp + ".tif");
	close();
	
	// Get MIPs of thin slices of the restack
	run("Make Substack...", "slices="+slice0);
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir+"reslice/slice0/MAX_slice0_c2_TP" + tp + ".tif");
	close();
	close();
	run("Make Substack...", "slices="+slice1);
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir+"reslice/slice1/MAX_slice1_c2_TP" + tp + ".tif");
	close();
	close();
	run("Make Substack...", "slices="+slice2);
	run("Z Project...", "projection=[Max Intensity]");
	saveAs("Tiff", datdir+"reslice/slice2/MAX_slice2_c2_TP" + tp + ".tif");
	close();
	close();
	
	// Optional: save merged 2color stack
	run("Merge Channels...", "c2=[Reslice of TP" + tp + "_Ch2_Ill0_Ang0,45,90,135,180,225,270,315] c6=[Reslice of TP" + tp + "_Ch1_Ill0_Ang0,45,90,135,180,225,270,315] create");
	saveAs("Tiff", datdir+"reslice/Composite_TP" + tp + ".tif");
		
	close();
	close();
	close();
	
}


// Assemble MIPs into movies
run("Image Sequence...", "dir=" + datdir + "reslice/mips/ filter=c1 sort");
saveAs("Tiff", datdir + "mips_reslice_c1.tif");

run("Image Sequence...", "dir=" + datdir + "reslice/mips/ filter=c2 sort");
saveAs("Tiff", datdir + "mips_reslice_c2.tif");

run("Merge Channels...", "c2=[mips_reslice_c2.tif] c6=[mips_reslice_c1.tif] create");
run("Rotate... ", "angle="+rotAngle+" grid=1 interpolation=Bilinear");
// Labels and scalebar
run("Label...", "format=00:00 starting=0 interval=3 x=5 y=20 font=24 text=[] range=1-11 use");
run("Set Scale...", "distance=1 known=0.5238 unit=um");
run("Scale Bar...", "width=50 height=50 thickness=6 font=24 color=White background=None location=[Lower Right] horizontal overlay");
saveAs("Tiff", datdir + "mips_reslice_Composite.tif");


// Assemble slice0 MIPs into movies
run("Image Sequence...", "dir=" + datdir + "reslice/slice0/ filter=c1 sort");
saveAs("Tiff", datdir + "mips_slice0_c1.tif");
run("Image Sequence...", "dir=" + datdir + "reslice/slice0/ filter=c2 sort");
saveAs("Tiff", datdir + "mips_slice0_c2.tif");
run("Merge Channels...", "c2=[mips_slice0_c2.tif] c6=[mips_slice0_c1.tif] create");
run("Rotate... ", "angle="+rotAngle+" grid=1 interpolation=Bilinear");
// Labels and scalebar
run("Label...", "format=00:00 starting=0 interval=3 x=5 y=20 font=24 text=[] range=1-11 use");
run("Set Scale...", "distance=1 known=0.5238 unit=um");
run("Scale Bar...", "width=50 height=50 thickness=6 font=24 color=White background=None location=[Lower Right] horizontal overlay");
saveAs("Tiff", datdir + "mips_slice0.tif");

//slice1
run("Image Sequence...", "dir=" + datdir + "reslice/slice1/ filter=c1 sort");
saveAs("Tiff", datdir + "mips_slice1_c1.tif");
run("Image Sequence...", "dir=" + datdir + "reslice/slice1/ filter=c2 sort");
saveAs("Tiff", datdir + "mips_slice1_c2.tif");
run("Merge Channels...", "c2=[mips_slice1_c2.tif] c6=[mips_slice1_c1.tif] create");
run("Rotate... ", "angle="+rotAngle+" grid=1 interpolation=Bilinear");
// Labels and scalebar
run("Label...", "format=00:00 starting=0 interval=3 x=5 y=20 font=24 text=[] range=1-11 use");
run("Set Scale...", "distance=1 known=0.5238 unit=um");
run("Scale Bar...", "width=50 height=50 thickness=6 font=24 color=White background=None location=[Lower Right] horizontal overlay");
saveAs("Tiff", datdir + "mips_slice1.tif");

//slice2
run("Image Sequence...", "dir=" + datdir + "reslice/slice2/ filter=c1 sort");
saveAs("Tiff", datdir + "mips_slice2_c1.tif");
run("Image Sequence...", "dir=" + datdir + "reslice/slice2/ filter=c2 sort");
saveAs("Tiff", datdir + "mips_slice2_c2.tif");
run("Merge Channels...", "c2=[mips_slice2_c2.tif] c6=[mips_slice2_c1.tif] create");
run("Rotate... ", "angle="+rotAngle+" grid=1 interpolation=Bilinear");
// Labels and scalebar
run("Label...", "format=00:00 starting=0 interval=3 x=5 y=20 font=24 text=[] range=1-11 use");
run("Set Scale...", "distance=1 known=0.5238 unit=um");
run("Scale Bar...", "width=50 height=50 thickness=6 font=24 color=White background=None location=[Lower Right] horizontal overlay");
saveAs("Tiff", datdir + "mips_slice2.tif");


close();
close();
close();

close();
close();
close();
	
