%Script used for connectivity analysis in manuscript

% Map showing the coral reef areas in the gulf
figure() 
hold on
worldmap([18 32],[-100 -80]) 
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
g(1) = geoshow('coral_gom3.shp','FaceColor', '#F0A9DD', 'FaceAlpha',1, 'EdgeColor', '#F0A9DD')

%Make a table from shapefile
saveGeoTable = readgeotable('coral_gom3.shp');
%There are 370 polygon areas in this shapefile

%Want to 1) create larger polygons that are a buffer around reef areas and
%2) divide the region into areas (WG -western gulf/veracruz, NY - northern
%yucatan/ scorpion reef, NM - northern mesoamerican reef, C - Cuba, FL-
%Florida Reef Tract, and FGB)

%Step 1 - buffer (only need to do this step 1x)
%https://www.mathworks.com/help/map/ref/bufferm.html

gom_coral = shaperead('coral_gom3.shp','UseGeoCoords',true)


    bufwidth = 0.05;
    direction = 'outPlusInterior';
for i = 101:length(gom_coral)
    corallat = gom_coral(i).Lat;
    corallon = gom_coral(i).Lon;
    [latbuf,lonbuf] = bufferm(corallat,corallon,bufwidth,direction);
    gom_coral(i).Latbuf = latbuf;
    gom_coral(i).Lonbuf = lonbuf;
end

%Save the new buffered lat lon as a shapefile so I can plot both on a map
%Remove original lat/lon
rmfield(gom_coral,'Lat');
rmfield(gom_coral,'Lon');

%Replace buffer names with just Lat/Lon
newName = 'Lat';
oldName = 'Latbuf';
[gom_coral.(newName)] = gom_coral.(oldName);
gom_coral = rmfield(gom_coral,oldName);

newName = 'Lon';
oldName = 'Lonbuf';
[gom_coral.(newName)] = gom_coral.(oldName);
gom_coral = rmfield(gom_coral,oldName);

%Create shapefile
shapewrite(gom_coral, "coral_buff.shp")


figure() 
hold on
worldmap([18 32],[-100 -80]) 
geoshow('coastL1.shp','FaceColor', '#D3D3D3', 'FaceAlpha',1) %alpha = transparentcy (1=opaque)
geoshow('coral_buff.shp','DisplayType','polygon','FaceColor', '#F0A9DD', 'FaceAlpha',1, 'EdgeColor', '#F0A9DD')
geoshow('coral_gom3.shp','FaceColor', 'red', 'FaceAlpha',1, 'EdgeColor', '#F0A9DD')

clear;clc

%%Step 2: Proportion of particles connected to buffered reef areas

%read in new buffered shapefile
coral_buff = shaperead('coral_buff.shp')

%Lon/lat data from cms model output stored in format: data_month = [lon_month, lat_month];

%December connectivity
in_dec_2001 = zeros(length(coral_buff),1);
on_dec_2001 = zeros(length(coral_buff),1);

for i = 1:length(coral_buff)
    [in, on] = inpolygon((data_dec(:,1)-360), data_dec(:,2),coral_buff(i).X, coral_buff(i).Y);
    in_dec_2001(i,1) = sum(in);
    on_dec_2001(i,1) = sum(on);
end

%October connectivity
in_oct_2001 = zeros(length(coral_buff),1);
on_oct_2001 = zeros(length(coral_buff),1);

for i = 1:length(coral_buff)
    [in, on] = inpolygon((data_oct(:,1)-360), data_oct(:,2),coral_buff(i).X, coral_buff(i).Y);
    in_oct_2001(i,1) = sum(in);
    on_oct_2001(i,1) = sum(on);
end

%July connectivity
in_jul_2001 = zeros(length(coral_buff),1);
on_jul_2001 = zeros(length(coral_buff),1);

for i = 1:length(coral_buff)
    [in, on] = inpolygon((data_jul(:,1)-360), data_jul(:,2),coral_buff(i).X, coral_buff(i).Y);
    in_jul_2001(i,1) = sum(in);
    on_jul_2001(i,1) = sum(on);
end

%April connectivity
in_apr_2001 = zeros(length(coral_buff),1);
on_apr_2001 = zeros(length(coral_buff),1);

for i = 1:length(coral_buff)
    [in, on] = inpolygon((data_apr(:,1)-360), data_apr(:,2),coral_buff(i).X, coral_buff(i).Y);
    in_apr_2001(i,1) = sum(in);
    on_apr_2001(i,1) = sum(on);
end

writematrix(in_apr_2001,'in_apr_2001.txt')
writematrix(on_apr_2001,'on_apr_2001.txt')

writematrix(in_jul_2001,'in_jul_2001.txt')
writematrix(on_jul_2001,'on_jul_2001.txt')

writematrix(in_oct_2001,'in_oct_2001.txt')
writematrix(on_oct_2001,'on_oct_2001.txt')

writematrix(in_dec_2001,'in_dec_2001.txt')
writematrix(on_dec_2001,'on_dec_2001.txt')

%Divide 370 polygons in coral shapefile into regions

coral_buff = readgeotable('coral_buff.shp')

%FGB - X < -90 & Y > 25
%WG - X < -93 > -88 & Y  < 25 18.4452,-94.7478
%Meso - X < -88 > -86 & Y <25
%Cuba+islands - X > 86 & Y < 23.5
%FL - X > -85 & Y > 23.5

%LINES 25 & 26 FGB
FGB = (coral_buff.X_coord < -90) & (coral_buff.Y_coord > 25);
disp(find(FGB))
numel(find(FGB))

%Western Gulf
WG = (coral_buff.X_coord < -93) & (coral_buff.Y_coord < 25);
disp(find(WG))
numel(find(WG))

%Northern Yucatan
NY = (coral_buff.X_coord > -93) & (coral_buff.X_coord < -88) & (coral_buff.Y_coord < 25);
disp(find(NY))
numel(find(NY))

%Meso - CHANGE TO 88.3
Meso = (coral_buff.X_coord > -88.3) & (coral_buff.X_coord < -85) & (coral_buff.Y_coord < 21);
disp(find(Meso))
numel(find(Meso))

%Cuba
C = (coral_buff.X_coord > -85) & (coral_buff.Y_coord < 23);
disp(find(C))
numel(find(C))

%Florida
FL = (coral_buff.X_coord > -83) & (coral_buff.Y_coord > 23);
disp(find(FL))
numel(find(FL))

%sum all connection for polygons in regions
%divide by (number of particles * number of timesteps) - proportion of all
%particle locations

% Example logical arrays (replace these with your actual logical data)
C = rand(370, 1) > 0.5;      % Example logical array C
FGB = rand(370, 1) > 0.5;    % Example logical array FGB
FL = rand(370, 1) > 0.5;     % Example logical array FL
Meso = rand(370, 1) > 0.5;   % Example logical array Meso
NY = rand(370, 1) > 0.5;     % Example logical array NY
WG = rand(370, 1) > 0.5;     % Example logical array WG

% Initialize an empty cell array to store results
combinedResults = {};

% Define the names of the dataframes
names = {'C', 'FGB', 'FL', 'Meso', 'NY', 'WG'};
dataframes = {C, FGB, FL, Meso, NY, WG};

% Iterate over each dataframe
for i = 1:length(dataframes)
    logicalArray = dataframes{i};
    dataframeName = names{i};
    
    % Create a new column with the dataframe name where logicalArray is 1
    newColumn = repmat({''}, 370, 1);  % Initialize with empty strings
    newColumn(logicalArray) = {dataframeName};  % Assign name where logical is 1
    
    % Concatenate to the combined results
    combinedResults = [combinedResults; newColumn];  % Append to combinedResults
end

% Convert combined results to a table
resultTable = table((1:370)', combinedResults, 'VariableNames', {'Index', 'Source'});

% Display the result
disp(resultTable);
