//// SET PARAMETERS ///////////////////////////////////////////////////////////////////

nStacks = 3 ; // Number of stacks total
stack0 = 0 ; // what is the index of the first stack? Usually zero.
fps = 5; // frames per second do you want the output movies to play at

midX_list = newArray("0620", "0720", "0820", "0920", "1020");

// Make sure to have a trailing forward slash at the end!
outputDir = "E:/avistrok/bapGAL4_UAShidUASstingerHiRFP/20251113121832_bapGAL4UAShidUASStingerHiRFP_3mpf_22x/";


// Now run step3

////////////////////////////////////////////////////////////
////////// SECTION III /////////////////////////////////////
////////////////////////////////////////////////////////////
// Readjust LUTs for midsection movies

// Now resave and save as AVIs

for (stck = stack0; stck < nStacks+stack0; stck++) {

	// Load midXs
	for (mx = 0; mx < midX_list.length; mx++) {
	    midX = midX_list[mx];
	    // Get the name of the current image sequence
	    baseName="Stack"+stck+"_left_composite_midX"+midX+".tif";
		selectImage(baseName);
	    // Construct the filenames for the AVI and Tiff files
	    aviPath = outputDir + "Stack"+stck+"_left_composite_midX"+midX+".avi";
	
	    // Save the files in the designated directory with the constructed filenames
		saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midX"+midX+".tif");
	    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
	    close();
	    
	    
	    // Get the name of the current image sequence
	    baseName="Stack"+stck+"_right_composite_midX"+midX+".tif";
		selectImage(baseName);
	    // Construct the filenames for the AVI and Tiff files
	    aviPath = outputDir + "Stack"+stck+"_right_composite_midX"+midX+".avi";
	
	    // Save the files in the designated directory with the constructed filenames
		saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midX"+midX+".tif");
	    run("AVI... ", "compression=JPEG frame="+fps+" save=" + aviPath);
	    close();
	}
}


