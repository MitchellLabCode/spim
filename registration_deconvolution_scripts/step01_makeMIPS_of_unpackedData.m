%% Make MIPS of unpacked data to preview the dataset.
% Here, we Split each half volume into its own MIP
% NPMitchell 2023

% datdir = '/mnt/crunch/actb2-mem-cherry_h2afva-gfp/202307111346_4views_1p4um_1ms/' ;
datdir = '/mnt/crunch/48YGAL4klarGFPnlsCAAXmCh/202308101028_180s_1p4um_2mW2mW_48YG4knlsGFPCAAXmCh_0p25_3p0msexposure/data/';
channels = [1,2] ;
angles = 0:45:359 ;
times = 0:101 ;

% Go through all files
for ang = angles 
    for ch =channels
        for tt = times 
            % Search for a timestamp
            searchfn = ['Time_' sprintf('%06d', tt) '_Angle_' num2str(ang) '_c' num2str(ch) '_ls_1.ome.tif'] ;
            fns = dir(fullfile(datdir, ...
                searchfn)) ;
            
            % Check that there is only one matching file
            assert(length(fns) == 1)
            ii = 1 ;

            % output filenames
            outdir1 = fullfile(datdir, 'mips', ['a' num2str(ang) '_c' num2str(ch) '_view11']) ;
            outdir2 = fullfile(datdir, 'mips', ['a' num2str(ang) '_c' num2str(ch) '_view12']) ;
            if ~exist(outdir1, 'dir')
                mkdir(outdir1)
                mkdir(outdir2)
            end

            % output filename
            nm = ['a' num2str(ang) '_c' num2str(ch) '_' sprintf('%06d', tt) '_0.png'] ;
            outfn1 = fullfile(outdir1, nm) ;
            outfn2 = fullfile(outdir2, nm) ;

            % Generate the mips & write to disk as pngs
            if exist(outfn1, 'file') 
                disp(['Warning: output MIP already exists: ' outfn1])
            elseif exist(outfn2, 'file')
                disp(['Warning: output MIP already exists: ' outfn2])
            else
                % Read the data for this timepoint
                disp(['reading ' fns(ii).name])
                fn = fullfile(fns(ii).folder, fns(ii).name) ;
                dat = tiffreadVolume(fn) ;
                half = size(dat, 3) ;
                mip1 = max(dat(:, :, 1:half), [], 3) ;
                mip2 = max(dat(:, :, half:end), [], 3) ;
                mip1 = imadjust(mip1) ;
                mip2 = imadjust(mip2) ;
        
                % Write the mips to disk as pngs
                disp(['writing ' outfn1])
                imwrite(mip1, outfn1)
                imwrite(mip2, outfn2)   
            end
            
        end
    end
end