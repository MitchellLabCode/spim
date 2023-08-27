// Run multi-pass registration on a single timepoint using the Multiview Reconstruction Plugin
// NPMitchell 2020-2023

// define dataset and path
datadir="/mnt/crunch/actb2-mem-cherry_h2afva-gfp/202307111346_4views_1p4um_1ms/";
dataset=datadir+"dataset.xml";
tp="74";


// basic registration line common to all passes
runline="select_xml="+dataset+" process_angle=[All angles] process_illumination=[All illuminations] process_timepoint=[Single Timepoint (Select from List)] ";
runline+="processing_timepoint=[Timepoint "+tp+"] ";

// customize the registration line for multiple passes with different methods
runline0=runline+"registration_algorithm=[Fast 3d geometric hashing (rotation invariant)] ";
runline0+="type_of_registration=[Register timepoints individually] interest_points_channel_1=beads interest_points_channel_2=beads fix_tiles=[Fix first tile] ";
runline0+="map_back_tiles=[Do not map back (use this if tiles are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 allowed_error_for_ransac=5 significance=10";

runline1=runline+"registration_algorithm=[Redundant geometric local descriptor matching (translation invariant)] ";
runline1+="type_of_registration=[Register timepoints individually] interest_points_channel_1=beads interest_points_channel_2=beads fix_tiles=[Fix first tile] ";
runline1+="map_back_tiles=[Do not map back (use this if tiles are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 number_of_neighbors=3 redundancy=1 significance=3 allowed_error_for_ransac=5";

runline2=runline+"registration_algorithm=[Iterative closest-point (ICP, no invariance)] ";
runline2+="type_of_registration=[Register timepoints individually] interest_points_channel_1=beads interest_points_channel_2=beads fix_tiles=[Fix first tile] ";
runline2+="map_back_tiles=[Do not map back (use this if tiles are fixed)] transformation=Affine regularize_model model_to_regularize_with=Rigid lamba=0.10 maximal_distance=5 maximal_number=100";

// run a series of registrations. repeating should improve the result
run("Register Dataset based on Interest Points", runline0);
run("Register Dataset based on Interest Points", runline0);
run("Register Dataset based on Interest Points", runline0);
run("Register Dataset based on Interest Points", runline0);
run("Register Dataset based on Interest Points", runline1);
run("Register Dataset based on Interest Points", runline2);
run("Register Dataset based on Interest Points", runline0);

// finally, fuse the result with some downsampling to inspect the quality of the result
run("Fuse/Deconvolve Dataset", "select_xml="+dataset+" process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_timepoint=[Single Timepoint (Select from List)] processing_timepoint=[Timepoint "+tp+"] type_of_image_fusion=[Weighted-average fusion] bounding_box=[Define manually] fused_image=[Save as TIFF stack] minimal_x=-674 minimal_y=-83 minimal_z=-519 maximal_x=2883 maximal_y=2068 maximal_z=2990 downsample=4 pixel_type=[16-bit unsigned integer] imglib2_container=[CellImg (large images)] process_views_in_paralell=All blend interpolation=[Linear Interpolation] output_file_directory="+datadir);



