// Define the data directory path
//datadir = "/mnt/crunch/48YGAL4klarLamGFPCAAXmCh/202208122137_180s_1p4um_2mW2mW_48YG4kLamGFPCAAXmCh_4views_0p25_0p5ms/";
datadir = "/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202208101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh/";

// Define the range of angles you want to iterate over
angles = newArray(0, 45, 90, 135, 180, 225, 270, 315);

for (i = 0; i < angles.length; i++) {
    angle = angles[i];

    // ch1
    run("Image Sequence...", "open=" + datadir + "mips/MAX_Time_000000_Angle_0_c1_ls_1.ome.tif file=Angle_"+angle+"_c1 sort");
    run("Rotate 90 Degrees Right");
    run("Enhance Contrast", "saturated=0.6");
    run("Flip Horizontally", "stack");
    saveAs("Tiff", datadir + "mips_Angle_"+angle+"_GFP.tif");
    run("AVI... ", "compression=JPEG frame=10 save=" + datadir + "mips_Angle"+angle+"_GFP.avi");

    // ch2
    run("Image Sequence...", "open=" + datadir + "mips/MAX_Time_000000_Angle_180_c2_ls_1.ome.tif file=Angle_"+angle+"_c2 sort");
    run("Rotate 90 Degrees Right");
    run("Enhance Contrast", "saturated=5");
    run("Flip Horizontally", "stack");
    saveAs("Tiff", datadir + "mips_Angle_"+angle+"_RFP.tif");
    run("AVI... ", "compression=JPEG frame=10 save=" + datadir + "mips_Angle"+angle+"_RFP.avi");

    run("Merge Channels...", "c2=mips_Angle_"+angle+"_GFP.tif c6=mips_Angle_"+angle+"_RFP.tif create keep");
    run("Set Scale...", "distance=1 known=0.2619 unit=um");
    run("Scale Bar...", "width=100 height=8 font=36 color=White background=None location=[Lower Right] bold overlay");
    run("Label...", "format=00:00 starting=0 interval=3 x=5 y=41 font=36 text=[] range=1-80 use");
    saveAs("Tiff", datadir + "/mips_Angle"+angle+"_Composite.tif");
    run("AVI... ", "compression=JPEG frame=10 save=" + datadir + "mips_Angle"+angle+"_Composite.avi");

    close();
    close();
    close();

    open(datadir+"mips_Angle"+angle+"_GFP.avi");
    open(datadir+"mips_Angle"+angle+"_RFP.avi");
    open(datadir+"mips_Angle"+angle+"_Composite.avi");

	run("Combine...", "stack1=mips_Angle"+angle+"_GFP.avi stack2=mips_Angle"+angle+"_RFP.avi combine");
	run("Combine...", "stack1=[Combined Stacks] stack2=mips_Angle"+angle+"_Composite.avi combine");
	run("Size...", "width=1024 height=1539 depth=100 constrain average interpolation=Bilinear");
    run("AVI... ", "compression=JPEG frame=10 save=" + datadir + "mips_Angle"+angle+"_Final.avi");
    
    close();
}



// selectWindow("ch2_"+angle+".tif");
// run("AVI... ", "compression=JPEG frame=10 save=" + datadir + "mips_Angle"+angle+"_RFP.avi");
// close();

// selectWindow("ch1_"+angle+".tif");
// run("AVI... ", "compression=JPEG frame=10 save=" + datadir + "mips_Angle"+angle+"_GFP.avi");
// close();