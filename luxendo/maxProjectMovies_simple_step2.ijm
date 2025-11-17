//// SET PARAMETERS ///////////////////////////////////////////////////////////////////
//
nStacks = 3 ; // Number of stacks total
stack0 = 0 ; // what is the index of the first stack? Usually zero.
fps = 5; // frames per second do you want the output movies to play at

// Make sure to have a trailing forward slash at the end!
outputDir = "E:/avistrok/bapGAL4_UAShidUASstingerHiRFP/20251113121832_bapGAL4UAShidUASStingerHiRFP_3mpf_22x/";


//////////////////////////////////////////////////////////
//////// SECTION II //////////////////////////////////////
//////////////////////////////////////////////////////////

// // Adjust LUTs for all MIP stacks
// /*

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

// */


// Now run step3


