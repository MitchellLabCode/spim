
rootdir = "D:/mblBootcamp/bynGAL4klar_UASmChCAAXUASGFPnls/20251011181532_bynGAL4klar_UASmChCAAX_UASGFPnls_combined/";
datdir = rootdir + "masked_deconvolved_18it_15it_bgsub/";
pix_per_um = 3.41880342 // 22x --> 3.41880342; 33x --> 5.1282;                  // pixels per micrometer
tps       = newArray(0, 60, 118);     // timepoints to process
anglestr  = "21,81,141,201,261,321";  // list of angles to process
dt = 2;
// --------------------------

// for (i=0; i < lengthOf(tps); i++){
i = 2 ;
tp = tps[i] ;

// STEP 1
//open(datdir+"TP"+tp+"_Ch0_Ill0_Ang"+anglestr+".tif");
//open(datdir+"TP"+tp+"_Ch1_Ill0_Ang"+anglestr+".tif");
//run("Merge Channels...", "c2=TP"+tp+"_Ch1_Ill0_Ang"+anglestr+".tif c6=TP"+tp+"_Ch0_Ill0_Ang"+anglestr+".tif create");

//
////run("Channels Tool...");
//Stack.setActiveChannels("10");
//setMinAndMax("0.015", "0.06");
//Stack.setActiveChannels("01");
//setMinAndMax("0.0199", "0.215");
//Stack.setActiveChannels("11");
//


// STEP 2
//run("3D Project...", "projection=[Brightest Point] axis=Y-Axis slice=1 initial=0 total=360 rotation=5 lower=0 upper=255 opacity=0 surface=100 interior=0 interpolate");
//selectImage("Projections of Composite");


// STEP 3
run("Remove Overlay") ;
run("Set Scale...", "distance="+pix_per_um+" known=1 unit=µm");
run("Scale Bar...", "width=50 height=8 thickness=12 font=60 location=[Upper Right] overlay");
run("Label...", "format=00:00 starting="+dt*tp+" interval=0 x=5 y=65 font=60 text=[] range=1-72 use");
saveAs("Tiff", rootdir+"3Projections_dtheta5_tp"+tp+".tif");

// }