// Set the scale in microns per pixel
scale = 0.19500001;

// parentDir = getDirectory("Choose the Parent Directory");
parentDir = "D:/rlondo/HandGFPbynGAL4klar_H2AGFPUASMyo1CRFP/2024-05-15_230441/mips/";
dirList = getFileList(parentDir);
nDirs = dirList.length;
print("Number of directories: " + nDirs);

for (d = 0; d < nDirs; d++) {
    dir = parentDir + dirList[d];
    print("Directory: " + dir);
    if (File.isDirectory(dir)) {
        list = getFileList(dir);
        print("Number of files: " + list.length);
        nIms = list.length;
        
        // Open all .tiff files as an image sequence
		print("Opening all .tiff files as an image sequence: " + dir);
		run("Image Sequence...", "open=[" + dir + "] sort use");
		
		
		//for (i = 0; i < nIms; i++) {
		//filePath = dir + "/" + list[i]; // Use forward slash for paths
		//print("File: " + filePath);
		//
		//// Check if file ends with .tif or .tiff
		//if (filePath.length() > 4) {
		//    // Get the last 4 and 5 characters of the file path
		//    lastFour = filePath.substring(filePath.length() - 4);
		//    lastFive = filePath.substring(filePath.length() - 5);
		//
		//    print("Last 4 characters: " + lastFour);
		//    print("Last 5 characters: " + lastFive);
		//
		//    // Check if the file ends with .tif or .tiff
		//    if (lastFour == ".tif" || lastFive == ".tiff") {
		//        print("Opening: " + filePath);
		//        open(filePath); // Use open instead of File.open
        
        // Process the image
        run("Enhance Contrast", "saturated=0.35");
        run("Label...", "format=00:00 starting=0 interval=2 x=50 y=50 font=36 text=[] range=1-120 use");
        run("Set Scale...", "distance=1 known=" + scale + " unit=um");
        run("Scale Bar...", "width=50 height=50 thickness=20 font=100 bold overlay");
        
        // Save the processed file paths
        aviPath = dir + ".avi";
        tiffPath = dir + "_processed.tif";

        run("AVI... ", "compression=JPEG frame=5 save=" + aviPath);
        saveAs("Tiff", tiffPath);
       
    }
}
