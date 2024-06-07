datdir="/mnt/data/antibodies_lightsheet/abda_1to100_568/202101241145_e2t2_w48YGal4klar_ABDA_C11_568/view4_1p4um_25x_1mW_exp1p5/data/";

process_acit="process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_timepoint=[All Timepoints]";

// DETECTION
xml="select_xml="+datdir+"dataset.xml ";
detectOpts=" type_of_interest_point_detection=Difference-of-Gaussian label_interest_points=beads downsample_images set_minimal_and_maximal_intensity subpixel_localization=[3-dimensional quadratic fit] interest_point_specification=[Advanced ...] downsample_xy=[Match Z Resolution (less downsampling)] downsample_z=1x minimal_intensity=200 maximal_intensity=800 sigma=1.80150 threshold=0.2020 find_maxima compute_on=[CPU (Java)]";
run("Detect Interest Points for Registration", xml+process_acit+detectOpts);

// FAST REGISTRATION
fastOptions=process_acit+"registration_algorithm=[Fast 3d geometric hashing (rotation invariant)] type_of_registration=[Register timepoints individually] interest_points_channel_0=beads fix_tiles=[Fix first tile] map_back_tiles=[Do not map back (use this if tiles are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 allowed_error_for_ransac=5 significance=10";
run("Register Dataset based on Interest Points", xml+fastOptions);
run("Register Dataset based on Interest Points", xml+fastOptions);
run("Register Dataset based on Interest Points", xml+fastOptions);

//check it
//run("Fuse/Deconvolve Dataset", "select_xml="+datdir+"dataset.xml "+process_acit+" type_of_image_fusion=[Weighted-average fusion] bounding_box=[Define manually] fused_image=[Save as TIFF stack] minimal_x=0 minimal_y=0 minimal_z=0 maximal_x=1400 maximal_y=1400 maximal_z=2000 downsample=2 pixel_type=[16-bit unsigned integer] imglib2_container=[CellImg (large images)] process_views_in_paralell=All blend interpolation=[Linear Interpolation] output_file_directory="+datdir);
//open(datdir+"TP0_Ch0_Ill0_Ang0,45,90,135,180,225,270,315.tif")
//run("Brightness/Contrast...");
//run("Enhance Contrast", "saturated=0.35");

// REDUNDANT REGISTRATION
redOptions=process_acit+" registration_algorithm=[Redundant geometric local descriptor matching (translation invariant)] type_of_registration=[Register timepoints individually] interest_points_channel_0=beads fix_tiles=[Fix first tile] map_back_tiles=[Do not map back (use this if tiles are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 number_of_neighbors=3 redundancy=1 significance=3 allowed_error_for_ransac=5";
redOptions=process_acit+" registration_algorithm=[Iterative closest-point (ICP, no invariance)] type_of_registration=[Register timepoints individually] interest_points_channel_0=beads fix_tiles=[Fix first tile] map_back_tiles=[Do not map back (use this if tiles are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 maximal_distance=5 maximal_number=100";
run("Register Dataset based on Interest Points", xml+redOptions);

// ICP REGISTRATION
icpOptions=process_acit+" registration_algorithm=[Iterative closest-point (ICP, no invariance)] type_of_registration=[Register timepoints individually] interest_points_channel_0=beads fix_tiles=[Fix first tile] map_back_tiles=[Do not map back (use this if tiles are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 maximal_distance=5 maximal_number=100";
run("Register Dataset based on Interest Points", xml+icpOptions);

// FUSE THE RESULT
run("Fuse/Deconvolve Dataset", "select_xml="+datdir+"dataset.xml "+process_acit+" type_of_image_fusion=[Weighted-average fusion] bounding_box=[Define manually] fused_image=[Save as TIFF stack] minimal_x=0 minimal_y=0 minimal_z=-500 maximal_x=1400 maximal_y=1400 maximal_z=2000 downsample=1 pixel_type=[16-bit unsigned integer] imglib2_container=[CellImg (large images)] process_views_in_paralell=All blend interpolation=[Linear Interpolation] output_file_directory="+datdir);
open(datdir+"TP0_Ch0_Ill0_Ang0,45,90,135,180,225,270,315.tif")
//run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");
