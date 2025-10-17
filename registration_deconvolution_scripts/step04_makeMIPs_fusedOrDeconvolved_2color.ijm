//////////////////////////////////////////////////////
// Make MIPS of fusions for select timepoints in tps
//////////////////////////////////////////////////////

datdir = "F:\\PROJECTS\\LightMicroscopyBootcamp2025\\bynGAL4klar_UASmChCAAXUASGFPnls\\20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined\\deconvolved18iter_15it\\";
//datdir = "/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/data_resavedViaMATLAB/fused2x/";
n = 119;
//tps = newArray(0,10,20,30,40,50,60,70,80,90,100,110,118);
tps = newArray(n);
for (i = 0; i < n; i++) {
	tps[i] = i;
}

// Options
dt = "2" ; // dt in minutes or seconds
pix_per_um = "0.195" ;
ch0='0' ;
ch1='1' ;
//Angs = "0,45,90,135,180,225,270,315" ;
Angs = "21,81,141,201,261,321" ;
dash= '\\';



//
////loop over the various camera view indices
//for (i = 0; i < tps.length; i++){
//	tp = tps[i] ;
//
//	// ch1
//	open(datdir + "TP" + tp + "_Ch"+ch0+"_Ill0_Ang"+Angs+".tif");
//	selectWindow("TP" + tp + "_Ch"+ch0+"_Ill0_Ang"+Angs+".tif");
//	run("Z Project...", "projection=[Max Intensity]");
//	saveAs("Tiff", datdir + "mips"+dash+"c"+ch0+"_TP" + tp + ".tif");
//	close();
//
//	// reslice laterally
//	selectWindow("TP" + tp + "_Ch"+ch0+"_Ill0_Ang"+Angs+".tif");
//	run("Reslice [/]...", "output=2.000 start=Left");
//	run("Z Project...", "projection=[Max Intensity]");
//	saveAs("Tiff", datdir + "mips_left"+dash+"c"+ch0+"_TP" + tp + ".tif");
//	close();
//	close();
//	close();
//	
//	// ch2
//	open(datdir + "TP" + tp + "_Ch"+ch1+"_Ill0_Ang"+Angs+".tif");
//	selectWindow("TP" + tp + "_Ch"+ch1+"_Ill0_Ang"+Angs+".tif");
//	run("Z Project...", "projection=[Max Intensity]");
//	saveAs("Tiff", datdir + "mips"+dash+"c"+ch1+"_TP" + tp + ".tif");
//	close();
//	
//	// reslice laterally
//	selectWindow("TP" + tp + "_Ch"+ch1+"_Ill0_Ang"+Angs+".tif");
//	run("Reslice [/]...", "output=2.000 start=Left");
//	run("Z Project...", "projection=[Max Intensity]");
//	saveAs("Tiff", datdir + "mips_left"+dash+"c"+ch1+"_TP" + tp + ".tif");
//	close();
//	close();
//	close() ;
//	
//	// Merge
//	//run("Merge Channels...", "c2=[MAX_TP" + tp + "_Ch"+ch1+"_Ill0_Ang"+Angs+".tif] c6=[MAX_TP" + tp + "_Ch"+ch0+"_Ill0_Ang"+Angs+".tif] create");
//	//saveAs("Tiff", "/mnt"+dash+"crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/data_resavedViaMATLAB/fused2x/mips"+dash+"c"+ch1+"_TP" + tp + ".tif");
//	//close();
//	
//}


run("Image Sequence...", "dir=" + datdir + "mips"+dash+" filter=c"+ch0+" sort");
saveAs("Tiff", datdir + "mips_c"+ch0+".tif");
run("Image Sequence...", "dir=" + datdir + "mips"+dash+" filter=c"+ch1+" sort");
saveAs("Tiff", datdir + "mips_c"+ch1+".tif");
run("Merge Channels...", "c2=[mips_c"+ch1+".tif] c6=[mips_c"+ch0+".tif] create");
// Labels and scalebar
run("Label...", "format=00:00 starting=0 interval="+dt+" x=5 y=20 font=24 text=[] range=1-200 use");
run("Set Scale...", "distance=1 known="+pix_per_um+" unit=um");
run("Scale Bar...", "width=50 height=50 thickness=6 font=24 color=White background=None location=[Lower Right] horizontal overlay");
saveAs("Tiff", datdir + "mips_Composite.tif");


run("Image Sequence...", "dir=" + datdir + "mips_left"+dash+" filter=c"+ch0+" sort");
saveAs("Tiff", datdir + "mips_left_c"+ch0+".tif");
run("Image Sequence...", "dir=" + datdir + "mips_left"+dash+" filter=c"+ch1+" sort");
saveAs("Tiff", datdir + "mips_left_c"+ch1+".tif");
run("Merge Channels...", "c2=[mips_left_c"+ch1+".tif] c6=[mips_left_c"+ch0+".tif] create");
// Labels and scalebar
run("Label...", "format=00:00 starting=0 interval="+dt+" x=5 y=20 font=24 text=[] range=1-200 use");
run("Set Scale...", "distance=1 known="+pix_per_um+" unit=um");
run("Scale Bar...", "width=50 height=50 thickness=6 font=24 color=White background=None location=[Lower Right] horizontal overlay");
saveAs("Tiff", datdir + "mips_left_Composite.tif");


