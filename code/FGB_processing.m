%%End to end processing script for FGB disease backtracking

%STEP 1: Seascape tidying script updated January 2024
%%This processes raw trajectory file outputs from CMS into a usable
%%structure - "bigstruct". Also created a quick map to visualize the
%%trajectory files.

%STEP 2:Generate particle density maps with 15/30/60 day non-convex hulls
%for all release periods (seasons)and annual aggregation. 

%STEP 3:Generates the matrix imported into R to develop Figure 3-
%%maximum distance traveled by particles. 

%This script was run for each year - find and replace the year to modify.

%STEP 1: Processing if NETCDF

trajlist = dir(fullfile(pwd,'/traj*'));

% Initialize the structure
bigstruct = struct();

% Loop over the range of NetCDF files
for i = 1:length(trajlist)
    % Construct the filename
    filename = sprintf('traj_file_%02d.nc', i);  % Adjust format for two digits
    
    % Check if the file exists
    if isfile(filename)
        % Read the variables from the NetCDF file
        bigstruct(i).time = ncread(filename, 'time');
        bigstruct(i).location = ncread(filename, 'location');
        bigstruct(i).lon = ncread(filename, 'lon');
        bigstruct(i).lat = ncread(filename, 'lat');
        bigstruct(i).depth = ncread(filename, 'depth');
        bigstruct(i).distance = ncread(filename, 'distance');
        bigstruct(i).exitcode = ncread(filename, 'exitcode');
        bigstruct(i).releasedate = ncread(filename, 'releasedate');
    else
        warning('File %s does not exist.', filename);
    end
end

%remove particles if depth is greater than 400m
for i = 1:length(bigstruct)
    if bigstruct(i).depth > 400
        % Get all field names
        fields = fieldnames(bigstruct);
        % Set all fields to NaN
        for j = 1:length(fields)
            bigstruct(i).(fields{j}) = NaN;
        end
    end
end


f0 =figure('Name','YYYY_traj')
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

exportgraphics(f0,'YYYY_traj.png','Resolution',600) %CHANGE YEAR


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
%coefficient 0.5 for intermediate enveloping of points (1 = most compact, 0
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
cmap = colormap(interp(256));%default # of colors in interp is 256
cmap(1,:) = [1,1,1];
colormap(cmap)

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
cmap = colormap(interp(256));%default # of colors in interp is 256
cmap(1,:) = [1,1,1];
colormap(cmap)

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
cmap = colormap(interp(256));%default # of colors in interp is 256
cmap(1,:) = [1,1,1];
colormap(cmap)

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
cmap = colormap(interp(256));%default # of colors in interp is 256
cmap(1,:) = [1,1,1];
colormap(cmap)

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

%MAKING FIGURE OF DISTANCE TRAVELED
% r = %This is for release lines
% p = 1:2500 %This is for the particles in each release

distmat = zeros(size(bigstruct,2),size(bigstruct(1).lon,2));
for r = 1:size(bigstruct,2)

    for p = 1:size(bigstruct(r).lon,2)
        pdists = zeros(size(bigstruct(r).lon,1),1);
        for i = 1:size(bigstruct(r).lon,1)
            pdists(i) = pos2dist([bigstruct(r).lon(1,p),bigstruct(r).lat(1,p)], [bigstruct(r).lon(i,p),bigstruct(r).lat(i,p)]); % If you don't save two variables it just outputs the first one (Haversine)
        end

        distmat(r,p) = max(pdists);
    end
end
writematrix(distmat,'dist_YYYY') %UPDATE WITH YEAR

save 'YYYY_processing.mat' '-v7.3'