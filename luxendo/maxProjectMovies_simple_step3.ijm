//// SET PARAMETERS ///////////////////////////////////////////////////////////////////

nStacks = 3 ; // Number of stacks total
stack0 = 0 ; // what is the index of the first stack? Usually zero.
fps = 5; // frames per second do you want the output movies to play at

// Make sure to have a trailing forward slash at the end!
outputDir = "E:/Alex/HandGFP48YGAL4UASmChCAAXHiRFP/2025-03-12_144525/" ; 
outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/20250327_tiltedEmbryo_shiftsInY/2025-03-27_180348_part4_adjustY/";


outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/2025-06-13_164028/";
outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/2025-06-13_153744/";
//outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/2025-06-13_164244/";
outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/2025-06-13_172111/";
outputDir = "E:/Chris/bapGAL4UAShidUASStingerH2aviRFP/2025-06-22_171054/";


// Now run step3

////////////////////////////////////////////////////////////
////////// SECTION III /////////////////////////////////////
////////////////////////////////////////////////////////////
// Readjust LUTs for midsection movies

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


