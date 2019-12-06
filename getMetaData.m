function meta = getMetaData(nameDummy,range)
    % getMetaData obtain the meta data from micro manager ome.tif stacks 
    % Micro manager stores files in a stack, that has no relation to the
    % actual stack size. 
    % 
    % name : name of file
    % range: Number of files in directory. Default 1 

    
    if nargin < 2
        range = 1; 
    end
    meta = struct('snum',[],'zpos',[],'etime',[],'chan',[],'stacksize',[],...
        'tempres',[],'spatialres',[],'nchan',[],'full',[]);
    name = [nameDummy,'.ome.tif'];
    
    for time = 0 : range-1
        
        if time == 0
        
        else   
            name = [nameDummy,'_',num2str(time),'.ome.tif'];
        end
        
        
        ncam  = 2; % specify the number of cameras; Could be optional!
        
        fi    = imfinfo(name);
        snum  = zeros(1,length(fi));
        zpos  = zeros(1,length(fi));
        etime = zeros(1,length(fi));
        chan  = zeros(1,length(fi));
        full  = cell(1,length(fi));
        
        for k = 1 : length(fi)

            % the micro manager specific metadata is stored in the
            % UnkownTags field of the file information structure. 
            v   = fi(k).UnknownTags(end).Value;
            
            full{k} = v;
            % extract the SliceIndex from the file information
            key = 'SliceIndex';
            sf  = strfind(v,key);
            snum(k) = sscanf(v(sf(1)+length(key)+2 : end) ,'%g',1);
            
            % extract the ZPosition in mum from the file information
%            key = 'ZPositionUm';
%            sf = strfind(v,key);
%            zpos(k) = sscanf(v(sf(1)+length(key)+2 : end) ,'%g',1);

            % extract the elapsed time from the file information
            key =  'ElapsedTime-ms';
            sf = strfind(v,key);
            etime(k) = sscanf(v(sf(1)+length(key)+2 : end) ,'%g',1);
        
            % extract the channel index
            key = '"ChannelIndex';
            sf = strfind(v,key);
            chan(k) = sscanf(v(sf(1)+length(key)+2 : end) ,'%g',1);
        end
    
        meta(time+1).snum = snum;
        meta(time+1).zpos = zpos;
        meta(time+1).etime = etime;
        meta(time+1).chan = chan;
        meta(time+1).full = full;
    end
    
    meta(1).stacksize = max(cat(2,meta.snum))+1;
    % We still don't know how many channels are stored. Below is a quick
    % and dirty way of figuring this out: 
    % for now we only have two channel data, so we only implement this
    % case. General will have a screen for more colours.
    meta(1).nchan = max(cat(2,meta.chan)+1)/ncam;
    % determine temporal and spatial resolution
    times = cat(2,meta.etime);
    meta(1).tempres    = round((times(meta(1).stacksize*ncam*meta(1).nchan)-times(1))/1000);
    meta(1).spatialres = (max(cat(2,meta.zpos))-min(cat(2,meta.zpos)))/max(cat(2,meta.snum));