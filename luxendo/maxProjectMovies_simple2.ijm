// Set the scale in microns per pixel
scale = 0.19500001;
interval = 2 ;
nStacks = 3 ;
fps = 7;

// Specify the output directory
outputDir = "D:/rlondo/HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/2024-05-17_184833/mips/";  // Change this to your desired output directory
outputDir = "E:/haibei/48YGAL4klar_UASmChCAAXHiRFP/2024-05-23_183541/mips/" ;
outputDir = "E:/rio/HandGFPbynGAL4klar_UASmChCAAXHiFP/2024-05-27_165010_22C_excellent/" ;
outputDir = "E:/rio/HandGFPbynGAL4klar_UASmChCAAXHiFP/2024-05-27_125722_22C_butFocusDrifts/" ;

// 3 color [cyan red white]
for (stck = 1; stck < nStacks+1; stck++) {

	run("Merge Channels...", "c1=stack_"+stck+"_channel_1_obj_left c4=stack_"+stck+"_channel_2_obj_left c5=stack_"+stck+"_channel_0_obj_left create");
	run("Enhance Contrast", "saturated=0.35");
	run("Label...", "format=00:00 starting=0 interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
	run("Set Scale...", "distance=1 known=" + scale + " unit=um");
	run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
	saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite.tif");

	run("Merge Channels...", "c1=stack_"+stck+"_channel_1_obj_right c4=stack_"+stck+"_channel_2_obj_right c5=stack_"+stck+"_channel_0_obj_right create");
	run("Enhance Contrast", "saturated=0.35");
	run("Label...", "format=00:00 starting=0 interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
	run("Set Scale...", "distance=1 known=" + scale + " unit=um");
	run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
	saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite.tif");
	
	
	run("Merge Channels...", "c1=stack_"+stck+"_channel_1_obj_left_midZ c4=stack_"+stck+"_channel_2_obj_left_midZ c5=stack_"+stck+"_channel_0_obj_left_midZ create");
	run("Enhance Contrast", "saturated=0.35");
	run("Label...", "format=00:00 starting=0 interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
	run("Set Scale...", "distance=1 known=" + scale + " unit=um");
	run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
	saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midZ.tif");


	run("Merge Channels...", "c1=stack_"+stck+"_channel_1_obj_right_midZ c4=stack_"+stck+"_channel_2_obj_right_midZ c5=stack_"+stck+"_channel_0_obj_right_midZ create");
	run("Enhance Contrast", "saturated=0.35");
	run("Label...", "format=00:00 starting=0 interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
	run("Set Scale...", "distance=1 known=" + scale + " unit=um");
	run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
	saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midZ.tif");
}

//Adjust LUTs for all MIP stacks
//
//// Now resave and save MIPs as AVIs
//for (stck = 1; d < 4; d++) {
//
//    // Get the name of the current image sequence
//	selectImage("Stack"+stck+"_left_composite.tif");
//    baseName = getTitle();
//    // Construct the filenames for the AVI and Tiff files
//    aviPath = outputDir + baseName + ".avi";
//
//    // Save the files in the designated directory with the constructed filenames
//	saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midZ.tif");
//    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
//    close();
//    
//    
//	selectImage("Stack"+stck+"_right_composite.tif");
//    baseName = getTitle();
//    // Construct the filenames for the AVI and Tiff files
//    aviPath = outputDir + baseName + ".avi";
//
//    // Save the files in the designated directory with the constructed filenames
//	saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midZ.tif");
//    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
//    close();
//}

// Readjust LUTs for midsection movies

// Now resave and save as AVIs
//for (stck = 1; stck < 4; stck++) {
//
//    // Get the name of the current image sequence
//    baseName="Stack"+stck+"_left_composite_midZ";
//	selectImage(baseName);
//    // Construct the filenames for the AVI and Tiff files
//    aviPath = outputDir + baseName + ".avi";
//
//    // Save the files in the designated directory with the constructed filenames
//	saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midZ.tif");
//    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
//    close();
//    
//    
//    // Get the name of the current image sequence
//    baseName="Stack"+stck+"_right_composite_midZ";
//	selectImage(baseName);
//    // Construct the filenames for the AVI and Tiff files
//    aviPath = outputDir + baseName + ".avi";
//
//    // Save the files in the designated directory with the constructed filenames
//	saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midZ.tif");
//    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
//    close();
//}
