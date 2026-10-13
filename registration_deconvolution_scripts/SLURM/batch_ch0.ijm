dataDir="/project/npmitchell/path/to/your/data_with_dataset.xml_in_it/";
minx="448";
miny="232";
minz="350";
maxx="1176";
maxy="792";
maxz="976";
niterations="40";
outdir=dataDir+"deconvolved32bit/";
mag="22x";
excitation="488";

passedArgument = call("java.lang.System.getenv", "SLURM_ARRAY_TASK_ID");

opts0="select_xml="+dataDir+"dataset.xml process_angle=[All angles] process_channel=[Single channel (Select from List)] " ;
opts1="process_illumination=[All illuminations] process_timepoint=[Multiple Timepoints (Select from List)] processing_channel=[channel 0] ";
opts1="timepoint_"+passedArgument+" type_of_image_fusion=[Multi-view deconvolution] bounding_box=[Define manually] fused_image=[Save as TIFF stack] ";
opts2="minimal_x="+minx+" minimal_y="+miny+" minimal_z="+minz+" maximal_x="+maxx+" maximal_y="+maxy+" maximal_z="+maxz+" ";
opts3="imglib2_container=[CellImg (large images)] imglib2_container_ffts=[CellImg (large images)] ";
opts4="save_memory type_of_iteration=[Efficient Bayesian - Optimization I (fast, precise)] image_weights=[Virtual weights (less memory, slower)] ";
opts5="osem_acceleration=[1 (balanced)] number_of_iterations="+niterations+" use_tikhonov_regularization tikhonov_parameter=0.0060 compute=[Entire image at once] ";
opts6="compute_on=[CPU (Java)] psf_estimation=[Provide file with PSF] psf_display=[Do not show PSFs] output_file_directory="+outdir+" ";
opts7=use_same_psf_for_all_angles/illuminations browse=/project/npmitchell/PSFs/"+mag+"/2um_waits_202407/mean_of_PSFs_"+excitation+".tif ";
opts8="transform_psfs psf_file=/project/npmitchell/PSFs/"+mag+"/2um_waits_202407/mean_of_PSFs_"+excitation+".tif "; 
options=opts0+opts1+opts2+opts3+opts4+opts5+opts6+opts7+opts8;

run("Fuse/Deconvolve Dataset", options);
eval("script", "System.exit(0);");
