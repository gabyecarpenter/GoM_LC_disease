%%End to end processing script for FGB disease backtracking

%STEP 1: Seascape tidying script updated January 2024
%%This processes raw trajectory file outputs from CMS into a usable
%%structure - "bigstruct".

%STEP 2:Generate particle density maps with 15/30/60 day non-convex hulls
%for all release periods (seasons) and annual aggregation. This code will
%save .mat files for data and .png images for maps generated.

%STEP 3:Generates the matrix imported into R to develop Figure 3-
%%maximum distance traveled by particles. Options for total distance
%%traveled and euclidian distance (use euclidian in manuscipt).

%STEP 4: Calculate connection of particles to reef areas. Reef areas are
%calculated in connectivity.m file (creates buffed polygons)

%%
% *%This script was run for each year - find and replace the year to modify.
% * 

%STEP 1: Tidying

%%Convert netcdf data to structure
trajlist = dir(fullfile(pwd,'/traj*'));

% Initialize the structure
bigstruct_0 = struct();

% Loop over the range of NetCDF files
for i = 1:length(trajlist)
    % Construct the filename
    filename = sprintf('traj_file_%02d.nc', i);  % Adjust format for two digits
    
    % Check if the file exists
    if isfile(filename)
        % Read the variables from the NetCDF file
        bigstruct_0(i).time = ncread(filename, 'time');
        bigstruct_0(i).location = ncread(filename, 'location');
        bigstruct_0(i).lon = ncread(filename, 'lon');
        bigstruct_0(i).lat = ncread(filename, 'lat');
        bigstruct_0(i).depth = ncread(filename, 'depth');
        bigstruct_0(i).distance = ncread(filename, 'distance');
        bigstruct_0(i).exitcode = ncread(filename, 'exitcode');
        bigstruct_0(i).releasedate = ncread(filename, 'releasedate');
    else
        warning('File %s does not exist.', filename);
    end
end

% bigstruct = bigstruct_0 
%Can create bigstruct immediatly if the number of release lines = number of
%traj files (not the case most of the time, so use the next loop)

%Depending on how many nodes/cores on the HPC, may need to adjust
%structure. This makes it so that each release line is a line in the
%structure. 

bigstruct = struct();
% bigstruct = struct('time', {}, 'location', {}, 'lon', {},'lat', {},'depth', {}, 'distance', {}, 'exitcode', {}, 'releasedate', {});

row_no = 0;
for i = 1:length(bigstruct_0)
    locs_in_file = unique(bigstruct_0(i).location);
    for j = 1:length(locs_in_file)
        row_no = row_no + 1
        r_index = bigstruct_0(i).location == locs_in_file(j);
        bigstruct(row_no).time = bigstruct_0(i).time;
        bigstruct(row_no).location = bigstruct_0(i).location(r_index);
        bigstruct(row_no).lon = bigstruct_0(i).lon(:,r_index);
        bigstruct(row_no).lat = bigstruct_0(i).lat(:,r_index);
        bigstruct(row_no).depth = bigstruct_0(i).depth(:,r_index);
        bigstruct(row_no).distance = bigstruct_0(i).distance(:,r_index);
        bigstruct(row_no).exitcode = bigstruct_0(i).exitcode(r_index);
        bigstruct(row_no).releasedate = bigstruct_0(i).releasedate(r_index);
    end
end



clear bigstruct_0

save traj_all.mat bigstruct -v7.3 %just saves the structure




%% Process ASCII output - SKIP IF YOUR OUTPUT IS IN NETCDF FORMAT
% The goal is to understand how each particle moves from its starting
% location. Leaving this code here in case there is a need to ouput in
% ASCII

% trajlist = dir('traj_file*'); %create list of all trajectory file in the directory, pwd = what is my current location
%     %can change to /output/traj if want to work in parent folder
% tic %big chunk is going to run through all trajectory files and bring all of them into one structure
% bigstruct = struct('loc',{},'lon',{},'lat',{},'dep',{},'flg',{},'tim',{});
% for k = 1:length(trajlist)
%     traj = trajlist(k).name;
%     trajfile = load(strcat(trajlist(k).folder,'/',traj),'-ascii');
% 
%     locs = unique(trajfile(:,1)); % Would help to have the vector of poly IDs/release lines from the release file here
%     outputrearrange = struct();
%     for l = 1:length(locs)
%         disp(k)
%         disp(l)
%         subtraj = trajfile(trajfile(:,1) == locs(l),:); %subset traj file by location
%         outputrearrange(l).loc = locs(l);
%         partcls = unique(subtraj(:,2));
%         times = unique(subtraj(:,3));
%         extstep = max(diff(times)); % This is likely the intended external time step
%         baddies = [];
%         baddies = find(rem(subtraj(:,3),extstep) > 0);
%         % Check flags
%         if (isempty(baddies))
%             disp("No bad timesteps")
%         elseif (any(subtraj(baddies,8) < 0))
%             disp('There are flags at uneven timesteps')
%             disp('Removing flagged lines from subtraj')
%             disp('Remaking particles and unique times')
%             subtraj(baddies,:) = [];
%             partcls = unique(subtraj(:,2)); % Just being safe
%             times = unique(subtraj(:,3));
%             disp('Finishing removing baddies')
%             if ~all(diff(times) == extstep)
%                 disp('Times still unequal STOP!')
%             end
%         else 
%             disp('You have a bigger problem times are bad flags are good')
%         end
%         newlon = nan(length(times),length(partcls));
%         newlat = nan(length(times),length(partcls));
%         newdep = nan(length(times),length(partcls));
%         newflg = nan(length(times),length(partcls));
%         newtim = nan(length(times),length(partcls));
%         for p = 1:length(partcls)
%             subp = subtraj(subtraj(:,2) == partcls(p),:); %subset again by particle
%             [~, I] = sort(subp(:,3));
%             subp = subp(I,:);
%             newlon(1:size(subp,1),p) = subp(:,4);
%             newlat(1:size(subp,1),p) = subp(:,5);
%             newdep(1:size(subp,1),p) = subp(:,6);
%             newflg(1:size(subp,1),p) = subp(:,8); % DH added 1/31/24
%             newtim(1:size(subp,1),p) = subp(:,3); % DH added 1/31/24
%         end
%         outputrearrange(l).lon = newlon;
%         outputrearrange(l).lat = newlat;
%         outputrearrange(l).dep = newdep;
%         outputrearrange(l).flg = newflg;
%         outputrearrange(l).tim = newtim;
%     end
%     bigstruct = [bigstruct,outputrearrange];
% end
% toc %bigstruct is the important one - have each of the locations (500) and a matrix for lat, long, and depth, (241 time stemps). 

% save traj_all.mat bigstruct -v7.3 %just saves the structure

%Quick figure to visualize trajectory files

figure()
hold on;
for i = 1:length(bigstruct) %count through each row of the structure and plot location
    plot([bigstruct(i).lon]-360, [bigstruct(i).lat]); 
end
axis equal
hold off;

% Let's bring in some land so we know what we're looking at.

S = shaperead('coastL1.shp');
for i = 1:length(S) %this for loop extracts the americas from a larger file
    Xloc = S(i).X;
    Yloc = S(i).Y;
    keepIndex = ~isnan(Xloc) & ~isnan(Yloc);
    Xloc = Xloc(keepIndex);
    Yloc = Yloc(keepIndex);
    fill(Xloc,Yloc,[.3 .3 .3]); hold on
end
for i = 1:length(bigstruct) %this plots our data
    plot([bigstruct(i).lon]-360, [bigstruct(i).lat]); 
end
axis equal
axis([-100,-80,14,32]) %this cuts it down


%STEP 2: Maps

% Initialize an empty cell array to store concatenated matrices
concatenatedData = cell(1, numel(bigstruct));
 
% Loop through each element of the struct array
for idx = 1:numel(bigstruct)
    % Get the matrix from the current element
    currentMatrix_lat_all = bigstruct(idx).lat; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlat_all{idx} = currentMatrix_lat_all(:);
end
 
% Concatenate all the column vectors into a single array
lat_all = cat(1, concatenatedlat_all{:});

 
% Initialize an empty cell array to store concatenated matrices
concatenatedData_all = cell(1, numel(bigstruct));
 
% Loop through each element of the struct array
for idx = 1:numel(bigstruct)
    % Get the matrix from the current element
    currentMatrix_lon_all = bigstruct(idx).lon; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlon_all{idx} = currentMatrix_lon_all(:);
end
 
% Concatenate all the column vectors into a single array
lon_all = cat(1, concatenatedlon_all{:});
data_all = [lon_all, lat_all];

save lat_all_YYYY.mat lat_all -v7.3
save lon_all_YYYY.mat lon_all -v7.3

%%Ok, this is where the loops is going to get annoying - I am going to do
%%seperate loops for each season, so that I can visualize seasonality

%DECEMBER RELEASE (lines 1-20)
% Initialize an empty cell array to store concatenated matrices
concatenatedData_lat_dec = cell(1, 20);
 
% Loop through each element of the struct array
for idx = 1:20
    % Get the matrix from the current element
    currentMatrix_lat_dec = bigstruct(idx).lat; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlat_dec{idx} = currentMatrix_lat_dec(:);
end
 
% Concatenate all the column vectors into a single array
lat_dec = cat(1, concatenatedlat_dec{:});

 
% Initialize an empty cell array to store concatenated matrices
concatenatedData_lon_dec = cell(1, 20);
 
% Loop through each element of the struct array
for idx = 1:20
    % Get the matrix from the current element
    currentMatrix_lon_dec = bigstruct(idx).lon; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlon_dec{idx} = currentMatrix_lon_dec(:);
end
 
% Concatenate all the column vectors into a single array
lon_dec = cat(1, concatenatedlon_dec{:});
data_dec = [lon_dec, lat_dec];

save data_dec_YYYY.mat data_dec -v7.3


%OCTOBER RELEASE (lines 21-40)
% Initialize an empty cell array to store concatenated matrices
concatenatedData_lat_oct = cell(1, 20);
 
% Loop through each element of the struct array
for idx = 21:40
    % Get the matrix from the current element
    currentMatrix_lat_oct = bigstruct(idx).lat; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlat_oct{idx} = currentMatrix_lat_oct(:);
end
 
% Concatenate all the column vectors into a single array
lat_oct = cat(1, concatenatedlat_oct{:});

 
% Initialize an empty cell array to store concatenated matrices
concatenatedData_lon_oct = cell(1, 20);
 
% Loop through each element of the struct array
for idx = 21:40
    % Get the matrix from the current element
    currentMatrix_lon_oct = bigstruct(idx).lon; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlon_oct{idx} = currentMatrix_lon_oct(:);
end
 
% Concatenate all the column vectors into a single array
lon_oct = cat(1, concatenatedlon_oct{:});
data_oct = [lon_oct, lat_oct];

save data_oct_YYYY.mat data_oct -v7.3

%JULY RELEASE (lines 41-60)
% Initialize an empty cell array to store concatenated matrices
concatenatedData_lat_jul = cell(1, 20);
 
% Loop through each element of the struct array
for idx = 41:60
    % Get the matrix from the current element
    currentMatrix_lat_jul = bigstruct(idx).lat; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlat_jul{idx} = currentMatrix_lat_jul(:);
end
 
% Concatenate all the column vectors into a single array
lat_jul = cat(1, concatenatedlat_jul{:});

 
% Initialize an empty cell array to store concatenated matrices
concatenatedData_lon_jul = cell(1, 20);
 
% Loop through each element of the struct array
for idx = 41:60
    % Get the matrix from the current element
    currentMatrix_lon_jul = bigstruct(idx).lon; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlon_jul{idx} = currentMatrix_lon_jul(:);
end
 
% Concatenate all the column vectors into a single array
lon_jul = cat(1, concatenatedlon_jul{:});
data_jul = [lon_jul, lat_jul];

save data_jul_YYYY.mat data_jul -v7.3

%APRIL RELEASE (lines 61-80)
% Initialize an empty cell array to store concatenated matrices
concatenatedData_lat_apr = cell(1, 20);
 
% Loop through each element of the struct array
for idx = 61:80
    % Get the matrix from the current element
    currentMatrix_lat_apr = bigstruct(idx).lat; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlat_apr{idx} = currentMatrix_lat_apr(:);
end
 
% Concatenate all the column vectors into a single array
lat_apr = cat(1, concatenatedlat_apr{:});

 
% Initialize an empty cell array to store concatenated matrices
concatenatedData_lon_apr = cell(1, 20);
 
% Loop through each element of the struct array
for idx = 61:80
    % Get the matrix from the current element
    currentMatrix_lon_apr = bigstruct(idx).lon; % Assuming the field name is 'matrix'
    % Concatenate the matrix into a column vector
    concatenatedlon_apr{idx} = currentMatrix_lon_apr(:);
end
 
% Concatenate all the column vectors into a single array
lon_apr = cat(1, concatenatedlon_apr{:});
data_apr = [lon_apr, lat_apr];

save data_apr_YYYY.mat data_apr -v7.3

%To generate 15/30/60 day non-convex hulls + maps

%ANNUAL
%lats and lons for 15 day non-convex hull - timestep 120 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_120_all = [bigstruct(:).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_120_all = lons_120_all(1:120,:);
lons_120_all = lons_120_all(:);

lats_120_all = [bigstruct(:).lat];
lats_120_all = lats_120_all(1:120,:);
lats_120_all = lats_120_all(:);

%lats and lons for 30 day non-convex hull - timestep 240 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_240_all = [bigstruct(:).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_240_all = lons_240_all(1:240,:);
lons_240_all = lons_240_all(:);

lats_240_all = [bigstruct(:).lat];
lats_240_all = lats_240_all(1:240,:);
lats_240_all = lats_240_all(:);

%lats and lons for 60 day non-convex hulls - timestep 480 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_480_all = [bigstruct(:).lon];
lons_480_all = lons_480_all(1:480,:);
lons_480_all = lons_480_all(:);

lats_480_all = [bigstruct(:).lat];
lats_480_all = lats_480_all(1:480,:);
lats_480_all = lats_480_all(:);

%%Create non-covex hulls with boundary function - remove nan's and make
%coefficient 0.7 for intermediate enveloping of points (1 = most compact, 0
%= least compact) - https://www.mathworks.com/help/matlab/ref/boundary.html

%15 DAY
%need to subest first so the index location in boundary line up
lons_120_all = lons_120_all(~isnan(lons_120_all));
lats_120_all = lats_120_all(~isnan(lats_120_all));

%now we can run boundary - 15 day
tic
k120_all = boundary(lons_120_all,lats_120_all,.7);
toc

save k120_all.mat k120_all -v7.3

%30 DAY
%need to subest first so the index location in boundary line up
lons_240_all = lons_240_all(~isnan(lons_240_all));
lats_240_all = lats_240_all(~isnan(lats_240_all));

%now we can run boundary - 30 days
tic
k240_all = boundary(lons_240_all,lats_240_all,.7);
toc

save k240_all.mat k240_all -v7.3

%60 DAY
%same steps for 60 day non-convex hull
lons_480_all = lons_480_all(~isnan(lons_480_all));
lats_480_all = lats_480_all(~isnan(lats_480_all));

tic
k480_all = boundary(lons_480_all,lats_480_all,.7);
toc

save k480_all.mat k480_all -v7.3

%%Now, work on getting the density data for "z" element of plot using hh3 - https://www.mathworks.com/help/stats/hist3.html

[hist c] = hist3(data_all, 'Nbins', [1 1]*500); %CHANGE THIS IF YOU WANT SEASONS! %using [hh3 c] to return bin center
hist(hist==0) = nan;
z_all = hist;

[X, Y] = meshgrid(c{1},...
    c{2}); %can use the bin centers to create a meshgrid (required for plotting pcolorm) - want to make sure that our color grid aligns with the histogram grid

%%Manipulate z (histogram frequencies) so it is easier to understand on the
%map

z1_all = z_all/(max(max(z_all))); %makes it a proportion(scale 0-1)

%z2_all = log(z_all); %log of z (makes it easier to see)

%Finally ready to create the figure! 
%cmap = colormap(interp(256));%default # of colors in interp is 256
%cmap(1,:) = [1,1,1];
%colormap(cmap)

f1 =figure('Name','YYYY_all') %CHANGE THIS FOR YEAR
colorbar
hold on
worldmap([18 32],[-100 -80]) %don't fo as far south
pcolorm(Y,X-360,z1_all') %" ' " flips it on the acis so it displays correctly %USE Z2 FOR LOG
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
g(1) = geoshow('coral_gom3.shp','FaceColor', '#F0A9DD', 'FaceAlpha',1, 'EdgeColor', '#F0A9DD')
g(2)=plotm(lats_120_all(k120_all),lons_120_all(k120_all),':','Color','#E02F8A', 'LineWidth',2, 'DisplayName','15 day non-convex hull')
g(3)=plotm(lats_240_all(k240_all),lons_240_all(k240_all),'--','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','30 day non-convex hull')
g(4)=plotm(lats_480_all(k480_all),lons_480_all(k480_all),'-','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','60 day non-convex hull')
legend ([g(1), g(2),g(3),g(4)], {'Reef areas', '15 Day', '30 Day', '60 Day'})
legend (g, 'Location', 'eastoutside')
set(gca,'ColorScale','log') %makes it the log of the proportion - best option for visualization and interpretation!
title("YYYY - All")

exportgraphics(f1,'YYYY_all.png','Resolution',600) %CHANGE YEAR


%DECEMBER RELEASE (1:20)
%lats and lons for 15 day non-convex hull - timestep 120 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_120_dec = [bigstruct(1:20).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_120_dec = lons_120_dec(1:120,:);
lons_120_dec = lons_120_dec(:);

lats_120_dec = [bigstruct(1:20).lat];
lats_120_dec = lats_120_dec(1:120,:);
lats_120_dec = lats_120_dec(:);

%lats and lons for 30 day non-convex hull - timestep 240 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_240_dec = [bigstruct(1:20).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_240_dec = lons_240_dec(1:240,:);
lons_240_dec = lons_240_dec(:);

lats_240_dec = [bigstruct(1:20).lat];
lats_240_dec = lats_240_dec(1:240,:);
lats_240_dec = lats_240_dec(:);

%lats and lons for 60 day non-convex hulls - timestep 480 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_480_dec = [bigstruct(1:20).lon];
lons_480_dec = lons_480_dec(1:480,:);
lons_480_dec = lons_480_dec(:);

lats_480_dec = [bigstruct(1:20).lat];
lats_480_dec = lats_480_dec(1:480,:);
lats_480_dec = lats_480_dec(:);

%%Create non-covex hulls with boundary function - remove nan's and make
%coefficient 0.5 for intermediate enveloping of points (1 = most compact, 0
%= least compact) - https://www.mathworks.com/help/matlab/ref/boundary.html

%15 DAY
%need to subest first so the index location in boundary line up
lons_120_dec = lons_120_dec(~isnan(lons_120_dec));
lats_120_dec = lats_120_dec(~isnan(lats_120_dec));

%now we can run boundary - 15 day
tic
k120_dec = boundary(lons_120_dec,lats_120_dec,.7);
toc

save k120_dec.mat k120_dec -v7.3

%30 DAY
%need to subest first so the index location in boundary line up
lons_240_dec = lons_240_dec(~isnan(lons_240_dec));
lats_240_dec = lats_240_dec(~isnan(lats_240_dec));

%now we can run boundary - 30 days
tic
k240_dec = boundary(lons_240_dec,lats_240_dec,.7);
toc

save k240_dec.mat k240_dec -v7.3

%60 DAY
%same steps for 60 day non-convex hull
lons_480_dec = lons_480_dec(~isnan(lons_480_dec));
lats_480_dec = lats_480_dec(~isnan(lats_480_dec));

tic
k480_dec = boundary(lons_480_dec,lats_480_dec,.7);
toc

save k480_dec.mat k480_dec -v7.3

%%Now, work on getting the density data for "z" element of plot using hh3 - https://www.mathworks.com/help/stats/hist3.html

[hist c] = hist3(data_dec, 'Nbins', [1 1]*500); %CHANGE THIS IF YOU WANT SEASONS! %using [hh3 c] to return bin center
hist(hist==0) = nan;
z_dec = hist;

[X, Y] = meshgrid(c{1},...
    c{2}); %can use the bin centers to create a meshgrid (required for plotting pcolorm) - want to make sure that our color grid aligns with the histogram grid

%%Manipulate z (histogram frequencies) so it is easier to understand on the
%map

z1_dec = z_dec/(max(max(z_dec))); %makes it a proportion(scale 0-1)

%z2_dec = log(z_dec); %log of z (makes it easier to see)

%Finally ready to create the figure! 
%cmap = colormap(interp(256));%default # of colors in interp is 256
%cmap(1,:) = [1,1,1];
%colormap(cmap)

f2 =figure('Name','YYYY_dec') %CHANGE THIS FOR YEAR
colorbar
hold on
worldmap([18 32],[-100 -80]) %don't fo as far south
pcolorm(Y,X-360,z1_dec') %" ' " flips it on the acis so it displays correctly %USE Z2 FOR LOG
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
g(1) = geoshow('coral_gom3.shp','FaceColor', '#F0A9DD', 'FaceAlpha',1, 'EdgeColor', '#F0A9DD')
g(2)=plotm(lats_120_dec(k120_dec),lons_120_dec(k120_dec),':','Color','#E02F8A', 'LineWidth',2, 'DisplayName','15 day non-convex hull')
g(3)=plotm(lats_240_dec(k240_dec),lons_240_dec(k240_dec),'--','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','30 day non-convex hull')
g(4)=plotm(lats_480_dec(k480_dec),lons_480_dec(k480_dec),'-','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','60 day non-convex hull')
legend ([g(1), g(2),g(3),g(4)], {'Reef areas', '15 Day', '30 Day', '60 Day'})
legend (g, 'Location', 'eastoutside')
set(gca,'ColorScale','log') %makes it the log of the proportion - best option for visualization and interpretation!
title("YYYY - December")

exportgraphics(f2,'YYYY_dec.png','Resolution',600) %CHANGE YEAR


%OCTOBER RELEASE (21:40)
%lats and lons for 15 day non-convex hull - timestep 120 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_120_oct = [bigstruct(21:40).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_120_oct = lons_120_oct(1:120,:);
lons_120_oct = lons_120_oct(:);

lats_120_oct = [bigstruct(21:40).lat];
lats_120_oct = lats_120_oct(1:120,:);
lats_120_oct = lats_120_oct(:);

%lats and lons for 30 day non-convex hull - timestep 240 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_240_oct = [bigstruct(21:40).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_240_oct = lons_240_oct(1:240,:);
lons_240_oct = lons_240_oct(:);

lats_240_oct = [bigstruct(21:40).lat];
lats_240_oct = lats_240_oct(1:240,:);
lats_240_oct = lats_240_oct(:);

%lats and lons for 60 day non-convex hulls - timestep 480 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_480_oct = [bigstruct(21:40).lon];
lons_480_oct = lons_480_oct(1:480,:);
lons_480_oct = lons_480_oct(:);

lats_480_oct = [bigstruct(21:40).lat];
lats_480_oct = lats_480_oct(1:480,:);
lats_480_oct = lats_480_oct(:);

%%Create non-covex hulls with boundary function - remove nan's and make
%coefficient 0.5 for intermediate enveloping of points (1 = most compact, 0
%= least compact) - https://www.mathworks.com/help/matlab/ref/boundary.html

%15 DAY
%need to subest first so the index location in boundary line up
lons_120_oct = lons_120_oct(~isnan(lons_120_oct));
lats_120_oct = lats_120_oct(~isnan(lats_120_oct));

%now we can run boundary - 15 day
tic
k120_oct = boundary(lons_120_oct,lats_120_oct,.7);
toc

save k120_oct.mat k120_oct -v7.3

%30 DAY
%need to subest first so the index location in boundary line up
lons_240_oct = lons_240_oct(~isnan(lons_240_oct));
lats_240_oct = lats_240_oct(~isnan(lats_240_oct));

%now we can run boundary - 30 days
tic
k240_oct = boundary(lons_240_oct,lats_240_oct,.7);
toc

save k240_oct.mat k240_oct -v7.3

%60 DAY
%same steps for 60 day non-convex hull
lons_480_oct = lons_480_oct(~isnan(lons_480_oct));
lats_480_oct = lats_480_oct(~isnan(lats_480_oct));

tic
k480_oct = boundary(lons_480_oct,lats_480_oct,.7);
toc

save k480_oct.mat k480_oct -v7.3

%%Now, work on getting the density data for "z" element of plot using hh3 - https://www.mathworks.com/help/stats/hist3.html

[hist c] = hist3(data_oct, 'Nbins', [1 1]*500); %CHANGE THIS IF YOU WANT SEASONS! %using [hh3 c] to return bin center
hist(hist==0) = nan;
z_oct = hist;

[X, Y] = meshgrid(c{1},...
    c{2}); %can use the bin centers to create a meshgrid (required for plotting pcolorm) - want to make sure that our color grid aligns with the histogram grid

%%Manipulate z (histogram frequencies) so it is easier to understand on the
%map

z1_oct = z_oct/(max(max(z_oct))); %makes it a proportion(scale 0-1)

%z2_oct = log(z_oct); %log of z (makes it easier to see)

%Finally ready to create the figure! 
%cmap = colormap(interp(256));%default # of colors in interp is 256
%cmap(1,:) = [1,1,1];
%colormap(cmap)

f3 =figure('Name','YYYY_oct') %CHANGE THIS FOR YEAR
colorbar
hold on
worldmap([18 32],[-100 -80]) %don't fo as far south
pcolorm(Y,X-360,z1_oct') %" ' " flips it on the acis so it displays correctly %USE Z2 FOR LOG
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
g(1) = geoshow('coral_gom3.shp','FaceColor', '#F0A9DD', 'FaceAlpha',1, 'EdgeColor', '#F0A9DD')
g(2)=plotm(lats_120_oct(k120_oct),lons_120_oct(k120_oct),':','Color','#E02F8A', 'LineWidth',2, 'DisplayName','15 day non-convex hull')
g(3)=plotm(lats_240_oct(k240_oct),lons_240_oct(k240_oct),'--','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','30 day non-convex hull')
g(4)=plotm(lats_480_oct(k480_oct),lons_480_oct(k480_oct),'-','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','60 day non-convex hull')
legend ([g(1), g(2),g(3),g(4)], {'Reef areas', '15 Day', '30 Day', '60 Day'})
legend (g, 'Location', 'eastoutside')
set(gca,'ColorScale','log') %makes it the log of the proportion - best option for visualization and interpretation!
title("YYYY - October")

exportgraphics(f3,'YYYY_oct.png','Resolution',600) %CHANGE YEAR

%JULY RELEASE (41:60)
%lats and lons for 15 day non-convex hull - timestep 120 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_120_jul = [bigstruct(41:60).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_120_jul = lons_120_jul(1:120,:);
lons_120_jul = lons_120_jul(:);

lats_120_jul = [bigstruct(41:60).lat];
lats_120_jul = lats_120_jul(1:120,:);
lats_120_jul = lats_120_jul(:);

%lats and lons for 30 day non-convex hull - timestep 240 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_240_jul = [bigstruct(41:60).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_240_jul = lons_240_jul(1:240,:);
lons_240_jul = lons_240_jul(:);

lats_240_jul = [bigstruct(41:60).lat];
lats_240_jul = lats_240_jul(1:240,:);
lats_240_jul = lats_240_jul(:);

%lats and lons for 60 day non-convex hulls - timestep 480 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_480_jul = [bigstruct(41:60).lon];
lons_480_jul = lons_480_jul(1:480,:);
lons_480_jul = lons_480_jul(:);

lats_480_jul = [bigstruct(41:60).lat];
lats_480_jul = lats_480_jul(1:480,:);
lats_480_jul = lats_480_jul(:);

%%Create non-covex hulls with boundary function - remove nan's and make
%coefficient 0.5 for intermediate enveloping of points (1 = most compact, 0
%= least compact) - https://www.mathworks.com/help/matlab/ref/boundary.html

%15 DAY
%need to subest first so the index location in boundary line up
lons_120_jul = lons_120_jul(~isnan(lons_120_jul));
lats_120_jul = lats_120_jul(~isnan(lats_120_jul));

%now we can run boundary - 15 day
tic
k120_jul = boundary(lons_120_jul,lats_120_jul,.7);
toc

save k120_jul.mat k120_jul -v7.3

%30 DAY
%need to subest first so the index location in boundary line up
lons_240_jul = lons_240_jul(~isnan(lons_240_jul));
lats_240_jul = lats_240_jul(~isnan(lats_240_jul));

%now we can run boundary - 30 days
tic
k240_jul = boundary(lons_240_jul,lats_240_jul,.7);
toc

save k240_jul.mat k240_jul -v7.3

%60 DAY
%same steps for 60 day non-convex hull
lons_480_jul = lons_480_jul(~isnan(lons_480_jul));
lats_480_jul = lats_480_jul(~isnan(lats_480_jul));

tic
k480_jul = boundary(lons_480_jul,lats_480_jul,.7);
toc

save k480_jul.mat k480_jul -v7.3

%%Now, work on getting the density data for "z" element of plot using hh3 - https://www.mathworks.com/help/stats/hist3.html

[hist c] = hist3(data_jul, 'Nbins', [1 1]*500); %CHANGE THIS IF YOU WANT SEASONS! %using [hh3 c] to return bin center
hist(hist==0) = nan;
z_jul = hist;

[X, Y] = meshgrid(c{1},...
    c{2}); %can use the bin centers to create a meshgrid (required for plotting pcolorm) - want to make sure that our color grid aligns with the histogram grid

%%Manipulate z (histogram frequencies) so it is easier to understand on the
%map

z1_jul = z_jul/(max(max(z_jul))); %makes it a proportion(scale 0-1)

%z2_jul = log(z_jul); %log of z (makes it easier to see)

%Finally ready to create the figure! 
%cmap = colormap(interp(256));%default # of colors in interp is 256
%cmap(1,:) = [1,1,1];
%colormap(cmap)

f4 =figure('Name','YYYY_jul') %CHANGE THIS FOR YEAR
colorbar
hold on
worldmap([18 32],[-100 -80]) %don't fo as far south
pcolorm(Y,X-360,z1_jul') %" ' " flips it on the acis so it displays correctly %USE Z2 FOR LOG
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
g(1) = geoshow('coral_gom3.shp','FaceColor', '#F0A9DD', 'FaceAlpha',1, 'EdgeColor', '#F0A9DD')
g(2)=plotm(lats_120_jul(k120_jul),lons_120_jul(k120_jul),':','Color','#E02F8A', 'LineWidth',2, 'DisplayName','15 day non-convex hull')
g(3)=plotm(lats_240_jul(k240_jul),lons_240_jul(k240_jul),'--','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','30 day non-convex hull')
g(4)=plotm(lats_480_jul(k480_jul),lons_480_jul(k480_jul),'-','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','60 day non-convex hull')
legend ([g(1), g(2),g(3),g(4)], {'Reef areas', '15 Day', '30 Day', '60 Day'})
legend (g, 'Location', 'eastoutside')
set(gca,'ColorScale','log') %makes it the log of the proportion - best option for visualization and interpretation!
title("YYYY - July")

exportgraphics(f4,'YYYY_jul.png','Resolution',600) %CHANGE YEAR


%APRIL RELEASE (61:80)
%lats and lons for 15 day non-convex hull - timestep 120 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_120_apr = [bigstruct(61:80).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_120_apr = lons_120_apr(1:120,:);
lons_120_apr = lons_120_apr(:);

lats_120_apr = [bigstruct(61:80).lat];
lats_120_apr = lats_120_apr(1:120,:);
lats_120_apr = lats_120_apr(:);

%lats and lons for 30 day non-convex hull - timestep 240 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_240_apr = [bigstruct(61:80).lon]; %could also change the ":" here to extract run lines 1:x for each release season
lons_240_apr = lons_240_apr(1:240,:);
lons_240_apr = lons_240_apr(:);

lats_240_apr = [bigstruct(61:80).lat];
lats_240_apr = lats_240_apr(1:240,:);
lats_240_apr = lats_240_apr(:);

%lats and lons for 60 day non-convex hulls - timestep 480 with ouput
%frequency ever 3 hrs (10,800 sec)
lons_480_apr = [bigstruct(61:80).lon];
lons_480_apr = lons_480_apr(1:480,:);
lons_480_apr = lons_480_apr(:);

lats_480_apr = [bigstruct(61:80).lat];
lats_480_apr = lats_480_apr(1:480,:);
lats_480_apr = lats_480_apr(:);

%%Create non-covex hulls with boundary function - remove nan's and make
%coefficient 0.5 for intermediate enveloping of points (1 = most compact, 0
%= least compact) - https://www.mathworks.com/help/matlab/ref/boundary.html

%15 DAY
%need to subest first so the index location in boundary line up
lons_120_apr = lons_120_apr(~isnan(lons_120_apr));
lats_120_apr = lats_120_apr(~isnan(lats_120_apr));

%now we can run boundary - 15 day
tic
k120_apr = boundary(lons_120_apr,lats_120_apr,.7);
toc

save k120_apr.mat k120_apr -v7.3

%30 DAY
%need to subest first so the index location in boundary line up
lons_240_apr = lons_240_apr(~isnan(lons_240_apr));
lats_240_apr = lats_240_apr(~isnan(lats_240_apr));

%now we can run boundary - 30 days
tic
k240_apr = boundary(lons_240_apr,lats_240_apr,.7);
toc

save k240_apr.mat k240_apr -v7.3

%60 DAY
%same steps for 60 day non-convex hull
lons_480_apr = lons_480_apr(~isnan(lons_480_apr));
lats_480_apr = lats_480_apr(~isnan(lats_480_apr));

tic
k480_apr = boundary(lons_480_apr,lats_480_apr,.7);
toc

save k480_apr.mat k480_apr -v7.3

%%Now, work on getting the density data for "z" element of plot using hh3 - https://www.mathworks.com/help/stats/hist3.html

[hist c] = hist3(data_apr, 'Nbins', [1 1]*500); %CHANGE THIS IF YOU WANT SEASONS! %using [hh3 c] to return bin center
hist(hist==0) = nan;
z_apr = hist;

[X, Y] = meshgrid(c{1},...
    c{2}); %can use the bin centers to create a meshgrid (required for plotting pcolorm) - want to make sure that our color grid aligns with the histogram grid

%%Manipulate z (histogram frequencies) so it is easier to understand on the
%map

z1_apr = z_apr/(max(max(z_apr))); %makes it a proportion(scale 0-1)

%z2_apr = log(z_apr); %log of z (makes it easier to see)

%Finally ready to create the figure! 
%cmap = colormap(interp(256));%default # of colors in interp is 256
%cmap(1,:) = [1,1,1];
%colormap(cmap)

f5 =figure('Name','YYYY_apr') %CHANGE THIS FOR YEAR
colorbar
hold on
worldmap([18 32],[-100 -80]) %don't fo as far south
pcolorm(Y,X-360,z1_apr') %" ' " flips it on the acis so it displays correctly %USE Z2 FOR LOG
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
g(1) = geoshow('coral_gom3.shp','FaceColor', '#F0A9DD', 'FaceAlpha',1, 'EdgeColor', '#F0A9DD')
g(2)=plotm(lats_120_apr(k120_apr),lons_120_apr(k120_apr),':','Color','#E02F8A', 'LineWidth',2, 'DisplayName','15 day non-convex hull')
g(3)=plotm(lats_240_apr(k240_apr),lons_240_apr(k240_apr),'--','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','30 day non-convex hull')
g(4)=plotm(lats_480_apr(k480_apr),lons_480_apr(k480_apr),'-','Color','#E02F8A', 'LineWidth',1.5, 'DisplayName','60 day non-convex hull')
legend ([g(1), g(2),g(3),g(4)], {'Reef areas', '15 Day', '30 Day', '60 Day'})
legend (g, 'Location', 'eastoutside')
set(gca,'ColorScale','log') %makes it the log of the proportion - best option for visualization and interpretation!
title("YYYY - April")

exportgraphics(f5,'YYYY_apr.png','Resolution',600) %CHANGE YEAR



%STEP 3: Distance - Fig. 3

%%Use this code if you have netcdf output and distance were calculated by
%%the cms- THIS IS TOTAL DISTANCE TRAVELED

% % Initialize the distance matrix with zeros
% distmat_15_YYYY = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));
% distmat_30_YYYY = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));
% distmat_60_YYYY = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));
% distmat_90_YYYY = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));
% 
% 
% % Loop through each structure in bigstruct
% for r = 1:size(bigstruct, 2)
%     % Access the distance matrix for the current structure
%     currentMatrix = bigstruct(r).distance; 
% 
%     % Find the maximum for each column and store it in distmat
%     for p = 1:size(currentMatrix, 2)
%         distmat_15_YYYY(r, p) = max(currentMatrix(1:120, p)); % Get max for 15 days
%         distmat_30_YYYY(r, p) = max(currentMatrix(1:240, p)); % Get max for 30 days
%         distmat_60_YYYY(r, p) = max(currentMatrix(1:480, p)); % Get max for 60 days
%         distmat_90_YYYY(r, p) = max(currentMatrix(1:721, p)); % Get max for 90 days
%     end
% end
% 
% % Now distmat contains the maximum distances for each of the 2500 columns for each of the 80 structures
% writematrix(distmat_15_YYYY,'dist_15_YYYY')
% writematrix(distmat_30_YYYY,'dist_30_YYYY')
% writematrix(distmat_60_YYYY,'dist_60_YYYY')
% writematrix(distmat_90_YYYY,'dist_90_YYYY')

% % Initialize the distance matrix with zeros
% distmat = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));

% % Loop through each structure in bigstruct
% for r = 1:size(bigstruct, 2)
%     % Access the distance matrix for the current structure
%     currentMatrix = bigstruct(r).distance; % use r instead of p
% 
%     % Find the maximum for each column and store it in distmat
%     for p = 1:size(currentMatrix, 2)
%         distmat(r, p) = max(currentMatrix(:, p)); % Get max of each column
%     end
% end
% 
% % Now distmat contains the maximum distances for each of the 2500 columns for each of the 80 structures
% writematrix(distmat,'dist_euc_YYYY')

%Interested in EUCLIDIAN DISTANCE - direct from point A to point B, so use
%the following code to calculate

r = 80 %This is for release lines
p = 1:2500; %This is for the particles in each release

% Initialize the distance matrix with zeros
eucmat_15_YYYY = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));
eucmat_30_YYYY = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));
eucmat_60_YYYY = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));
eucmat_90_YYYY = zeros(size(bigstruct, 2), size(bigstruct(1).distance, 2));

%Loop for 15 days
for r = 1:size(bigstruct,2)
    for p = 1:size(bigstruct(r).lon,2)
        pdists = zeros(size(bigstruct(r).lon,1),1);
        for i = 1:120 %Timestep for 15 days
            pdists_15(i) = pos2dist([bigstruct(r).lon(1,p),bigstruct(r).lat(1,p)], [bigstruct(r).lon(i,p),bigstruct(r).lat(i,p)]); % If you don't save two variables it just outputs the first one (Haversine)
        end
        eucmat_15_YYYY(r,p) = max(pdists_15);     
    end
end

%Loop for 30 days
for r = 1:size(bigstruct,2)
    for p = 1:size(bigstruct(r).lon,2)
        pdists = zeros(size(bigstruct(r).lon,1),1);
        for i = 1:240 %Timestep for 30 days
            pdists_30(i) = pos2dist([bigstruct(r).lon(1,p),bigstruct(r).lat(1,p)], [bigstruct(r).lon(i,p),bigstruct(r).lat(i,p)]); % If you don't save two variables it just outputs the first one (Haversine)
        end
        eucmat_30_YYYY(r,p) = max(pdists_30);     
    end
end

%Loop for 60 days
for r = 1:size(bigstruct,2)
    for p = 1:size(bigstruct(r).lon,2)
        pdists = zeros(size(bigstruct(r).lon,1),1);
        for i = 1:480 %Timestep for 60 days
            pdists_60(i) = pos2dist([bigstruct(r).lon(1,p),bigstruct(r).lat(1,p)], [bigstruct(r).lon(i,p),bigstruct(r).lat(i,p)]); % If you don't save two variables it just outputs the first one (Haversine)
        end
        eucmat_60_YYYY(r,p) = max(pdists_60);     
    end
end

%Loop for 90 days
for r = 1:size(bigstruct,2)
    for p = 1:size(bigstruct(r).lon,2)
        pdists = zeros(size(bigstruct(r).lon,1),1);
        for i = 1:721 %Timestep for 60 days
            pdists_90(i) = pos2dist([bigstruct(r).lon(1,p),bigstruct(r).lat(1,p)], [bigstruct(r).lon(i,p),bigstruct(r).lat(i,p)]); % If you don't save two variables it just outputs the first one (Haversine)
        end
        eucmat_90_YYYY(r,p) = max(pdists_90);     
    end
end

writematrix(eucmat_15_YYYY,'euc_15_YYYY')
writematrix(eucmat_30_YYYY,'euc_30_YYYY')
writematrix(eucmat_60_YYYY,'euc_60_YYYY')
writematrix(eucmat_90_YYYY,'euc_90_YYYY')


%STEP 4
%%Proportion of particles connected to buffered reef areas

%THIS SECTION ASSUMES YOU HAVE ALREADY CALCULATED YOUR POLYGONS FOR
%CONNECTIONS - I do this in the connectivity.m file


%read in new buffered shapefile
coral_buff = shaperead('coral_buff.shp')

%Lon/lat data from cms model output stored in format: data_month = [lon_month, lat_month];

%December connectivity
in_dec_YYYY = zeros(length(coral_buff),1);
on_dec_YYYY = zeros(length(coral_buff),1);

for i = 1:length(coral_buff)
    [in, on] = inpolygon((data_dec(:,1)-360), data_dec(:,2),coral_buff(i).X, coral_buff(i).Y);
    in_dec_YYYY(i,1) = sum(in);
    on_dec_YYYY(i,1) = sum(on);
end

%October connectivity
in_oct_YYYY = zeros(length(coral_buff),1);
on_oct_YYYY = zeros(length(coral_buff),1);

for i = 1:length(coral_buff)
    [in, on] = inpolygon((data_oct(:,1)-360), data_oct(:,2),coral_buff(i).X, coral_buff(i).Y);
    in_oct_YYYY(i,1) = sum(in);
    on_oct_YYYY(i,1) = sum(on);
end

%July connectivity
in_jul_YYYY = zeros(length(coral_buff),1);
on_jul_YYYY = zeros(length(coral_buff),1);

for i = 1:length(coral_buff)
    [in, on] = inpolygon((data_jul(:,1)-360), data_jul(:,2),coral_buff(i).X, coral_buff(i).Y);
    in_jul_YYYY(i,1) = sum(in);
    on_jul_YYYY(i,1) = sum(on);
end

%April connectivity
in_apr_YYYY = zeros(length(coral_buff),1);
on_apr_YYYY = zeros(length(coral_buff),1);

for i = 1:length(coral_buff)
    [in, on] = inpolygon((data_apr(:,1)-360), data_apr(:,2),coral_buff(i).X, coral_buff(i).Y);
    in_apr_YYYY(i,1) = sum(in);
    on_apr_YYYY(i,1) = sum(on);
end

writematrix(in_apr_YYYY,'in_apr_YYYY.txt')
writematrix(on_apr_YYYY,'on_apr_YYYY.txt')

writematrix(in_jul_YYYY,'in_jul_YYYY.txt')
writematrix(on_jul_YYYY,'on_jul_YYYY.txt')

writematrix(in_oct_YYYY,'in_oct_YYYY.txt')
writematrix(on_oct_YYYY,'on_oct_YYYY.txt')

writematrix(in_dec_YYYY,'in_dec_YYYY.txt')
writematrix(on_dec_YYYY,'on_dec_YYYY.txt')

save 'YYYY_processing.mat' '-v7.3'