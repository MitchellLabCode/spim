// Set the scale in microns per pixel
scale = 0.19500001;

// Specify the output directory
outputDir = "D:/rlondo/HandGFPbynGAL4klar_UASmChCAAXH2AviRFP/2024-05-17_184833/mips/";  // Change this to your desired output directory
outputDir = "E:/haibei/48YGAL4klar_UASmChCAAXHiRFP/2024-05-23_183541/mips/" ;

for (d = 0; d < 36; d++) {

    // Get the name of the current image sequence
    baseName = getTitle();
    // baseName = substring(imgName, 0, lengthOf(imgName) - 4);  // Remove the file extension

    run("Enhance Contrast", "saturated=0.35");
    run("Label...", "format=00:00 starting=0 interval=3 x=50 y=50 font=100 text=[] range=1-120 use");
    run("Set Scale...", "distance=1 known=" + scale + " unit=um");
    run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");

    // Construct the filenames for the AVI and Tiff files
    aviPath = outputDir + baseName + ".avi";
    //tiffPath = outputDir + baseName + "_processed.tif";

    // Save the files in the designated directory with the constructed filenames
    run("AVI... ", "compression=JPEG frame=5 save=" + aviPath);
    //saveAs("Tiff", tiffPath);
    close();
}
