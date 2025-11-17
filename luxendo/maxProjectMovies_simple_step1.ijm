//// SET PARAMETERS ///////////////////////////////////////////////////////////////////
// Set the scale in microns per pixel
// scale = 0.19500001; // 33x magnification has 0.195 microns / pixel
scale = 0.2925 ; // scale for 22.2x
interval = 3 ; // minutes per dt 
nStacks = 3 ; //
nChannels = 2 ; // number of Channels
stack0 = 0 ; // what is the index of the first stack? Usually zero.
time0 = 0 ; // first frame's timestamp in minutes
doMIPs = 0 ;
doMidZ = 0 ;
doMidX = 1 ;
midX_list = newArray("0620", "0720", "0820", "0920", "1020");

minCh = 1;
flipRight = true ;

if (nChannels == 2) {
	// 2 color [green magenta]
	grnCh = minCh ;
	magCh = minCh+1 ;
	// grnCh = 1 ;
	// magCh = 2 ;
}
if (nChannels == 3) {
	// 3 color [cyan red white]
	cyanCh = 0 ;
	redCh = 1 ;
	whiteCh = 2 ;
}

// Specify the output directory
// Make sure to have a trailing forward slash at the end!
//outputDir = "E:/Alex/HandGFP48YGAL4UASmChCAAXHiRFP/2025-03-12_144525/" ; 
//outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/20250327_tiltedEmbryo_shiftsInY/20250326_died/2025-03-26_215303_part4/";


//outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/2025-06-13_164028/";
//outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/2025-06-13_153744/";
//outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/2025-06-13_164244/";
//outputDir = "E:/boris/bynGAL4klar_UASmChCAAXHiRFP/2025-06-17_135805/";

//outputDir = "E:/Chris/bapGAL4UAShidUASStingerH2aviRFP/2025-06-22_171054/";
outputDir = "E:/avistrok/bapGAL4_UAShidUASstingerHiRFP/20251113121832_bapGAL4UAShidUASStingerHiRFP_3mpf_22x/";

//// END OF PARAMETERS SECTION //////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////
//// SECTION 0: Load individual TIFF stacks /////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

for (stck = stack0; stck < nStacks+stack0; stck++) {
	for (ch = minCh; ch < nChannels+minCh; ch++) {
		if (doMIPs) {
			File.openSequence(outputDir + "mips/stack_" + stck + "_channel_" + ch + "_obj_left_mips/");
			File.openSequence(outputDir + "mips/stack_" + stck + "_channel_" + ch + "_obj_right_mips/");
		}
		if (doMidZ) {
			File.openSequence(outputDir + "midZ/stack_" + stck + "_channel_" + ch + "_obj_left_midZ/");
			File.openSequence(outputDir + "midZ/stack_" + stck + "_channel_" + ch + "_obj_right_midZ/");
		}
				
		// Load midXs
		for (mx = 0; mx < midX_list.length; mx++) {
		    midX = midX_list[mx];
		
		    for (stck = stack0; stck < nStacks+stack0; stck++) {
		        for (ch = minCh; ch < nChannels+minCh; ch++) {
		
		            if (doMIPs) {
		                File.openSequence(outputDir + "mips/stack_" + stck + "_channel_" + ch + "_obj_left_mips/");
		                File.openSequence(outputDir + "mips/stack_" + stck + "_channel_" + ch + "_obj_right_mips/");
		            }
		
		            if (doMidZ) {
		                File.openSequence(outputDir + "midZ/stack_" + stck + "_channel_" + ch + "_obj_left_midZ/");
		                File.openSequence(outputDir + "midZ/stack_" + stck + "_channel_" + ch + "_obj_right_midZ/");
		            }
		
		            if (doMidX) {
		                File.openSequence(outputDir + "midX" + midX + "/stack_" + stck + "_channel_" + ch + "_obj_left_midX" + midX + "/");
		                File.openSequence(outputDir + "midX" + midX + "/stack_" + stck + "_channel_" + ch + "_obj_right_midX" + midX + "/");
		            }
		        }
		    }
		}



	}
}


/////////////////////////////////////////////////////////////////////////////////////
//// SECTION I: Collate the different colors together into multi-color TIFF stacks //
/////////////////////////////////////////////////////////////////////////////////////

// /*
print(nChannels); 
if (nChannels == 2) {
	for (stck = stack0; stck < nStacks + stack0; stck++) {
		if (doMIPs) {
			run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_left_mips c6=stack_"+stck+"_channel_"+magCh+"_obj_left_mips create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_mips.tif");
		
		
			run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_right_mips c6=stack_"+stck+"_channel_"+magCh+"_obj_right_mips create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			if (flipRight) {
				run("Flip Horizontally");
			}
			saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_mips.tif");
			
		}
		if (doMidZ) {
			run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_left_midZ c6=stack_"+stck+"_channel_"+magCh+"_obj_left_midZ create");	
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midZ.tif");
		
		
			run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_right_midZ c6=stack_"+stck+"_channel_"+magCh+"_obj_right_midZ create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			if (flipRight) {
				run("Flip Horizontally");
			}
			saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midZ.tif");
		}
		
		if (doMidX) {
			for (mx = 0; mx < midX_list.length; mx++) {
			    midX = midX_list[mx];
				
				run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_left_midX"+midX+" c6=stack_"+stck+"_channel_"+magCh+"_obj_left_midX"+midX+" create");	
				run("Enhance Contrast", "saturated=0.35");
				run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
				run("Set Scale...", "distance=1 known=" + scale + " unit=um");
				run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
				saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midX"+midX+".tif");
			
			
				run("Merge Channels...", "c2=stack_"+stck+"_channel_"+grnCh+"_obj_right_midX"+midX+" c6=stack_"+stck+"_channel_"+magCh+"_obj_right_midX"+midX+" create");
				run("Enhance Contrast", "saturated=0.35");
				run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
				run("Set Scale...", "distance=1 known=" + scale + " unit=um");
				run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
				if (flipRight) {
					run("Flip Horizontally");
				}
				saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midX"+midX+".tif");
			}
		}
	}
}
if (nChannels == 3) {
	
	for (stck = 0+stack0; stck < nStacks+stack0; stck++) {
			
		if (doMIPs) {
			run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_left_mips c4=stack_"+stck+"_channel_"+whiteCh+"_obj_left_mips c5=stack_"+stck+"_channel_"+cyanCh+"_obj_left_mips create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_mips.tif");
	
			run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_right_mips c4=stack_"+stck+"_channel_"+whiteCh+"_obj_right_mips c5=stack_"+stck+"_channel_"+cyanCh+"_obj_right_mips create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			if (flipRight) {
				run("Flip Horizontally");
			}	
			saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_mips.tif");
		
		}

		
		if (doMidZ) {
			run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_left_midZ c4=stack_"+stck+"_channel_"+whiteCh+"_obj_left_midZ c5=stack_"+stck+"_channel_"+cyanCh+"_obj_left_midZ create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midZ.tif");
		
			run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_right_midZ c4=stack_"+stck+"_channel_"+whiteCh+"_obj_right_midZ c5=stack_"+stck+"_channel_"+cyanCh+"_obj_right_midZ create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			if (flipRight) {
				run("Flip Horizontally");
			}
			saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midZ.tif");
		}
		
		if (doMidX) {
			run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_left_midX"+midX+" c4=stack_"+stck+"_channel_"+whiteCh+"_obj_left_midX"+midX+" c5=stack_"+stck+"_channel_"+cyanCh+"_obj_left_midX"+midX+" create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			saveAs("Tiff", outputDir + "Stack"+stck+"_left_composite_midX"+midX+".tif");
		
			run("Merge Channels...", "c1=stack_"+stck+"_channel_"+redCh+"_obj_right_midX"+midX+" c4=stack_"+stck+"_channel_"+whiteCh+"_obj_right_midX"+midX+" c5=stack_"+stck+"_channel_"+cyanCh+"_obj_right_midX"+midX+" create");
			run("Enhance Contrast", "saturated=0.35");
			run("Label...", "format=00:00 starting="+time0+" interval="+interval+" x=50 y=50 font=100 text=[] range=1-1200 use");
			run("Set Scale...", "distance=1 known=" + scale + " unit=um");
			run("Scale Bar...", "width=100 height=50 thickness=20 font=100 overlay");
			if (flipRight) {
				run("Flip Horizontally");
			}
			saveAs("Tiff", outputDir + "Stack"+stck+"_right_composite_midX"+midX+".tif");
		}

}
// */

// Now run step2

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


