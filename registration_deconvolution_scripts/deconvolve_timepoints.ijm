
datdir = "" ;
tps = "" ;

command = "browse="+datdir+"dataset.xml select_xml="+datdir+"dataset.xml process_angle=[All angles] process_channel=[All channels] process_illumination=[All illuminations] process_timepoint=[Range of Timepoints (Specify by Name)] process_following_timepoints=";
command+=tps+" type_of_image_fusion=[Multi-view deconvolution] bounding_box=[Define manually] fused_image=[Save as TIFF stack] ";
command+="minimal_x=128 minimal_y=674 minimal_z=156 maximal_x=838 maximal_y=1950 maximal_z=784 ";
command+="imglib2_container=[CellImg (large images)] imglib2_container_ffts=ArrayImg ";
command+="type_of_iteration=[Efficient Bayesian - Optimization I (fast, precise)] ";
command+="image_weights=[Virtual weights (less memory, slower)] osem_acceleration=[1 (balanced)] number_of_iterations=12 ";
command+="use_tikhonov_regularization tikhonov_parameter=0.0060 compute=[specify maximal blocksize manually] ";
command+="compute_on=[GPU (Nvidia CUDA via JNA)] psf_estimation=[Provide file with PSF] psf_display=[Do not show PSFs] ";
command+="output_file_directory="+datdir+" ";
command+="blocksize_x=800 blocksize_y=1400 blocksize_z=800 cuda_directory=/opt/Fiji.app/lib/linux64";

run("Fuse/Deconvolve Dataset", command);