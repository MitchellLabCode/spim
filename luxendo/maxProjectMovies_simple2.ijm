//// SET PARAMETERS ///////////////////////////////////////////////////////////////////
// Set the scale in microns per pixel
scale = 0.19500001; // 33x magnification has 0.195 microns / pixel
// scale = 0.2925 ; // scale for 22.2x
interval = 2 ; // minutes per dt 
nStacks = 3 ; //
nChannels = 3 ; // number of Channels
fps = 5; // frames per second do you want the output movies to play at
stack0 = 0 ; // what is the index of the first stack? Usually zero.
time0 = 22 ; // first frame's timestamp in minutes
//minval_mip = newArray(100, 100, 100); // Example min LUT values for 3 channels
//maxval_mip = newArray(375, 9000, 6000); // Example max LUT values for 3 channels
//minval_midZ = newArray(100, 100, 100); // Example min LUT values for 3 channels
//maxval_midZ = newArray(300, 230, 600); // Example max LUT values for 3 channels

if (nChannels == 2) {
	// 2 color [green magenta]
	// grnCh = 0 ;
	// magCh = 1 ;
	grnCh = 1 ;
	magCh = 2 ;
}
if (nChannels == 3) {
	// 3 color [cyan red white]
	cyanCh = 0 ;
	redCh = 1 ;
	whiteCh = 2 ;
}

// Specify the output directory
// Make sure to have a trailing forward slash at the end!
// outputDir = "E:/Chris/HandGFPbynGAL4klar_UASmChCAAXHiRFP/2024-07-02_174410/";
outputDir = "E:/wenjie/HandGFP48YGAL4klar_UASmChCAAXHiRFP/2024-12-22_201937/"; 
// outputDir = "/mnt/data/midgut_chirality/HandGFPbynGAL4klar_UASMyo1CRFPHiRFP/2024-05-29_182916_22C_excellent_120s/"



//// END OF PARAMETERS SECTION //////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
//// SECTION 0: Load individual TIFF stacks /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

// /*
for (stck = stack0; stck < nStacks+stack0; stck++) {
	for (ch = 0; ch < nChannels; ch++) {
		File.openSequence(outputDir + "mips/stack_" + stck + "_channel_" + ch + "_obj_left_mips/");
		File.openSequence(outputDir + "midZ/stack_" + stck + "_channel_" + ch + "_obj_left_midZ/");
		File.openSequence(outputDir + "mips/stack_" + stck + "_channel_" + ch + "_obj_right_mips/");
		File.openSequence(outputDir + "midZ/stack_" + stck + "_channel_" + ch + "_obj_right_midZ/");
	}
}


/////////////////////////////////////////////////////////////////////////////////////
//// SECTION I: Collate the different colors together into multi-color TIFF stacks //
/////////////////////////////////////////////////////////////////////////////////////


print(nChannels); 
if (nChannels == 2) {
	for (stck = stack0; stck < nStacks + stack0; stck++) {
	
		run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_left_mips c6=stack_"+stck+"_channel_"+magCh+"_obj_left_mips create");
		run("Enhance Contrast", "saturated=0.35");
		run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
		run("Set Scale...", "distance=1 known=" + scale + " unit=um");
		run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
		saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_mips.tif");
	
	
		run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_right_mips c6=stack_"+stck+"_channel_"+magCh+"_obj_right_mips create");
		run("Enhance Contrast", "saturated=0.35");
		run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
		run("Set Scale...", "distance=1 known=" + scale + " unit=um");
		run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
		saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_mips.tif");
		
	
		run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_left_midZ c6=stack_"+stck+"_channel_"+magCh+"_obj_left_midZ create");	
		run("Enhance Contrast", "saturated=0.35");
		run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
		run("Set Scale...", "distance=1 known=" + scale + " unit=um");
		run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
		saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midZ.tif");
	
	
		run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_right_midZ c6=stack_"+stck+"_channel_"+magCh+"_obj_right_midZ create");
		run("Enhance Contrast", "saturated=0.35");
		run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
		run("Set Scale...", "distance=1 known=" + scale + " unit=um");
		run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
		saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midZ.tif");
	}
}
if (nChannels == 3) {
	
	for (stck = 0+stack0; stck < nStacks+stack0; stck++) {
	
		run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_left_mips c4=stack_"+stck+"_channel_"+whiteCh+"_obj_left_mips c5=stack_"+stck+"_channel_"+cyanCh+"_obj_left_mips create");
		run("Enhance Contrast", "saturated=0.35");
		run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
		run("Set Scale...", "distance=1 known=" + scale + " unit=um");
		run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
		
		//// Loop through each channel and set the min and max values
		//for (i = 1; i <= nChannels; i++) {
		//    // Select the channel
		//    setChannel(i);
		//    
		//    // Set the min and max for the current channel
		//    setMinAndMax(minval_mip[i-1], maxval_mip[i-1]);
		//    
		//    // Apply the LUT to make the changes effective
		//    run("Apply LUT");
		//}		

		saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_mips.tif");
	
		run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_right_mips c4=stack_"+stck+"_channel_"+whiteCh+"_obj_right_mips c5=stack_"+stck+"_channel_"+cyanCh+"_obj_right_mips create");
		run("Enhance Contrast", "saturated=0.35");
		run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
		run("Set Scale...", "distance=1 known=" + scale + " unit=um");
		run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
		
		//// Loop through each channel and set the min and max values
		//for (i = 1; i <= nChannels; i++) {
		//    // Select the channel
		//    setChannel(i);
		//    
		//    // Set the min and max for the current channel
		//    setMinAndMax(minval_mip[i-1], maxval_mip[i-1]);
		//    
		//    // Apply the LUT to make the changes effective
		//    run("Apply LUT");
		//}		
		
		saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_mips.tif");
		
		
		run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_left_midZ c4=stack_"+stck+"_channel_"+whiteCh+"_obj_left_midZ c5=stack_"+stck+"_channel_"+cyanCh+"_obj_left_midZ create");
		run("Enhance Contrast", "saturated=0.35");
		run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
		run("Set Scale...", "distance=1 known=" + scale + " unit=um");
		run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
		
		//// Loop through each channel and set the min and max values
		//
		//for (i = 1; i <= nChannels; i++) {
		//    // Select the channel
		//    setChannel(i);
		//    
		//    // Set the min and max for the current channel
		//    setMinAndMax(minval_midZ[i-1], maxval_midZ[i-1]);
		//    
		//    // Apply the LUT to make the changes effective
		//    run("Apply LUT");
		//}		
		
		saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midZ.tif");
	
		run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_right_midZ c4=stack_"+stck+"_channel_"+whiteCh+"_obj_right_midZ c5=stack_"+stck+"_channel_"+cyanCh+"_obj_right_midZ create");
		run("Enhance Contrast", "saturated=0.35");
		run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-120 use");
		run("Set Scale...", "distance=1 known=" + scale + " unit=um");
		run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
		
		//// Loop through each channel and set the min and max values		
		//for (i = 1; i <= nChannels; i++) {
		//    // Select the channel
		//    setChannel(i);
		//    
		//    // Set the min and max for the current channel
		//    setMinAndMax(minval_midZ[i-1], maxval_midZ[i-1]);
		//    
		//    // Apply the LUT to make the changes effective
		//    run("Apply LUT");
		//}		
		
		saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midZ.tif");
	}
}
// */

//////////////////////////////////////////////////////////
//////// SECTION II //////////////////////////////////////
//////////////////////////////////////////////////////////

// // Adjust LUTs for all MIP stacks
/*

// // Now resave and save MIPs as AVIs
for (stck = stack0; stck < nStacks+stack0; stck++) {
   
	// // CAMERA LEFT
    // Get the name of the current image sequence
	selectImage("Stack"+stck+"_left_composite_mips.tif");
	//// Loop through each channel and set the min and max values
	//for (i = 1; i <= nChannels; i++) {
	//    // Select the channel
	//    Stack.setChannel(i);
	//    
	//    // Set the min and max for the current channel
	//    setMinAndMax(minval_mip[i-1], maxval_mip[i-1]);
	//    
	//    // Apply the LUT to make the changes effective
	//    run("Apply LUT");
	//}
	
    // Construct the filenames for the AVI and Tiff files
    aviPath = outputDir + "Stack"+stck+"_left_composite_mips.avi";

    // Save the files in the designated directory with the constructed filenames
	saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_mips.tif");
    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
    close();
    
    // // CAMERA RIGHT
	selectImage("Stack"+stck+"_right_composite_mips.tif");
	//// Loop through each channel and set the min and max values
	//for (i = 1; i <= nChannels; i++) {
	//    // Select the channel
	//    Stack.setChannel(i);
	//    
	//    // Set the min and max for the current channel
	//    setMinAndMax(minval_mip[i-1], maxval_mip[i-1]);
	//    
	//    // Apply the LUT to make the changes effective
	//    run("Apply LUT");
	//}
    // Construct the filenames for the AVI and Tiff files
    aviPath = outputDir + "Stack"+stck+"_right_composite_mips.avi";

    // Save the files in the designated directory with the constructed filenames
	  saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_mips.tif");
    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
    close();
}

*/

////////////////////////////////////////////////////////////
////////// SECTION III /////////////////////////////////////
////////////////////////////////////////////////////////////
// Readjust LUTs for midsection movies
/* 
// Now resave and save as AVIs
for (stck = stack0; stck < nStacks+stack0; stck++) {

    // Get the name of the current image sequence
    baseName="Stack"+stck+"_left_composite_midZ.tif";
	selectImage(baseName);
    // Construct the filenames for the AVI and Tiff files
    aviPath = outputDir + "Stack"+stck+"_left_composite_midZ.avi";

    // Save the files in the designated directory with the constructed filenames
	saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midZ.tif");
    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
    close();
    
    
    // Get the name of the current image sequence
    baseName="Stack"+stck+"_right_composite_midZ.tif";
	selectImage(baseName);
    // Construct the filenames for the AVI and Tiff files
    aviPath = outputDir + "Stack"+stck+"_right_composite_midZ.avi";

    // Save the files in the designated directory with the constructed filenames
	saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midZ.tif");
    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
    close();
}

*/

