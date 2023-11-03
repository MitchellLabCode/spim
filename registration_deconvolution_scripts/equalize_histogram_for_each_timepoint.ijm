// Get the number of slices (timepoints) in the stack
nn=nSlices();

outdir="/mnt/crunch/48YGAL4klarLamGFPCAAXmCh/202308122137_180s_1p4um_2mW2mW_48YG4kLamGFPCAAXmCh_4views_0p25_0p5ms/dataLamGFP_resavedViaMATLAB/deconvolution_test/" ;

// Loop through each timepoint (slice)
for (slice = 1; slice <= nn; slice++) {
    // Set the current slice (timepoint)
    // selectImage(slice);

    // Get the current image
    run("Make Substack...", "slices="+slice);

    // Calculate the histogram
    //run("Histogram");
	//
    // Get the histogram data
    //getHistogram(lower, upper, min, max);
    //
	// Calculate the cumulative distribution function (CDF)
	//cdf = newArray(256);
	//cdf[0] = min;
	//for (i = 0; i < 255; i++) {
	//    cdf[i + 1] = cdf[i] + (upper[i] - lower[i]);
	//}
	//
	// Calculate the 1% and 99% intensity thresholds
	//totalPixels = getWidth() * getHeight();
	//lowerThreshold = cdf[0] + (0.01 * (totalPixels - 1));
	//upperThreshold = cdf[0] + (0.99 * (totalPixels - 1));

    // Apply histogram equalization with custom thresholds
    run("Enhance Contrast...", "saturated=0.35");
    
    // Save the adjusted image with a name based on the slice number
    saveAs("Tiff", outdir+"Slice_"+slice+".png");

    // Close the duplicate image
    close();

    // Update the progress bar
    // setProgress(slice / nn);
}

// Reset the progress bar
//setProgressBar(100);

// Select the original stack
//selectImage(id, 1);
