% Author: Nils Tack
% Compiled with Matlab 2022b
% Script used to prepare masks from image sequences before performing PIV
% computations, and plot the resulting velocity and vorticity fields

%% Clear workspace
close all; clearvars

%% Set paths
path_main = 'E:\Nils\MATLAB scripts\veloVortPlot\Data'; % set path to 'Data' folder

% Other paths
path_PIVraw = fullfile(path_main,'PIVraw');         % raw PIV data
path_image = fullfile(path_main,'images');          % image sequence
path_PIVclean = fullfile(path_main,'PIVclean');     % clean velocity data without headers
path_masksBW = fullfile(path_main,'masksBW');       % BW masks produced from original images
path_outlines = fullfile(path_main,'outlines');     % outlines with x and y coordinates produced from BW masks (.csv)
% path_pressure = fullfile(path_main,'pressure');     % pressure data (.csv)
% path_forces = fullfile(path_main,'forces');         % force data

% Set file type for image sequence (.jpg or .tif or .png)
imFileType = '.jpg'; % !!! warning: png cause issues with colormap overlay.

%% Select images based on the increment used for PIV data - also rename images
incrPIV = 1;% same as increment originally set in DaVis to calculate velocity fields
selectImages(path_image,incrPIV,imFileType); % select images corresponding to each velocity field. Also renames files starting at 0001. File name becomes 'B".

%% OPTIONAL - Rotate images if necessary
% This step is completely optional. The orientation will not affect the
% creation of masks and the calculation of the pressure fields.
% To save time, it is best to orient all the original frames properly
% outside Matlab. The software used by Photron (PFV4) can be used to export
% rotated image sequences. WARNING! Rotating images using Microsoft’s
% rotate tool straight from the folder where all the images are found does
% NOT work.

opts.rotateSequence = 0;  % option to enable rotation of image sequence. 0 = no rotation, 1 = 90° rotation counter-clockwise; created as a failsafe to prevent accidental rotation

if opts.rotateSequence == 1
rotateImages(path_image,imFileType,180) % (folder containing images, file type, rotation angle); positive angle = CCW, negative angles = CW
end

%% Options for BW masks
opts.inverseMask    = 0; % inverse the BW masks to make white areas black and vice-versa. Required if subject is dark and background white (brightfield PIV)
opts.maskImage      = 0; % keep only area of interest (uses ginput)
opts.checkMasks     = 1; % check the quality of the masks as they get exported
opts.plotComparison = 1; % plots mask-making steps (muted automatically for mask export)
opts.exportForDaVis = 1; % export inversed masks for DaVis because Davis uses black (binary 0) to mask an area

%% Test mask and image processing options
% Read 1 test image
i1 = 50; % test frame number
D_image = dir(fullfile(path_image,['*',imFileType])); %extract images from image folder

% display title as image ID
temp = split(D_image(i1).name,'.');
filename = temp{1}; clear temp;
disp(filename)

% Import image data
I = importdata(quickfilepath(D_image(i1)));
channelID = 1; % RGB channel of the image; check camera specs for which layer contains the most information. Most sensors are more sensitive to green, but for PIV, this will depend on the wavelength of the laser
I = I(:,:,channelID); % select only one of the three layers of the image; 

% Option to inverse the original grayscale image. makeMasks3 requires subject to be white on a black background -
% use this option if brightfield PIV.
if opts.inverseMask
    figure; hold on
    I = imcomplement(I);
    imshow(I)
    title('Inverse colors')
end

% Masking unnecessary parts of the image to isolate the subkect
if opts.maskImage
    [I,x_ginput,y_ginput] = maskMasked(I); % select two oposing corners of a rectangle including the area of interest
end

% Subplots to see if masking settings look good
% Open 'makeMask3' function to change parameters if necessary
% Parameters include masking unwanted areas like tank/flume edges or
% reflections, contrast and gamma corrections, smoothing options.
I4 = makeMask3(I, opts);

%% Export all the BW masks
opts.plotComparison = 0; % enable or disable comparison figures to check the quality of the masks as they get exported
% -> this is important to ensure that the PIV data is computed properly
% without artifacts at the bounday

if opts.checkMasks % check individual masks as they get exported
    figTemp = figure;
end

fprintf('Exporting masks...');

for i=1:length(D_image)
    
    % display title as image ID
    temp = split(D_image(i).name,'.');
    filename = temp{1}; clear temp;

    % Read original image data
    I = importdata(quickfilepath(D_image(i)));
    I = I(:,:,channelID); % Convert to grayscale using only one of the RGB channels
    
    % Option to inverse the original grayscale image. makeMasks3 requires subject to be white on a black background -
    % use this option if brightfield PIV.
    if opts.inverseMask
        I = imcomplement(I);
        imshow(I)
        title('Inverse colors')
    end

        % Masking unnecessary parts of the image to isolate the subject
    if opts.maskImage
        ITemp = uint8(zeros(length(I(:,1)),length(I(1,:)))); % initiate new mask as black
        tempI = I(y_ginput(1):y_ginput(2),x_ginput(1):x_ginput(2));
        ITemp(y_ginput(1):y_ginput(2),x_ginput(1):x_ginput(2)) = tempI;
        I = ITemp;
    end
    
    % Make a binary mask
     I4 = makeMask3(I,opts);
    
    % Inverse the BW masks for DaVis
     if opts.exportForDaVis
     I4 = imcomplement(I4);
     end

    % Option to check masks during export
    if opts.checkMasks
       imshow(I4)
       title(sprintf('Frame %i',i))
       pause(0.05);
       clf
    end
        
    % Export mask
    filenameBW = sprintf('BW_%05g',i);
    imwrite(I4,fullfile(path_masksBW,[filenameBW,'.tif'])) % exports as tif     
end

close(figTemp)
fprintf('done\n');

%% Export outlines for other masking purposes (pressure computations, plotting)
% options
opts.multiOutlines = 1; % option to enable the storage of multiple outlines in the same file; otherwise uses only the biggest blob created during the BW mask export step 

% set BW mask directory
D_BW = dir(fullfile(path_masksBW,'*.tif')); 

% Make outlines from mask
f1 = 1; % test frame number
BWtemp = flipud(importdata(quickfilepath(D_BW(f1)))); % import BW mask
BWdiag = bwmorph(BWtemp,'diag'); % diagonal fill
BWtempDilate = bwmorph(BWdiag,'thicken',1); % dilate mask by 1px so the outline falls exactly on the edge of the mask
BWedge = bwboundaries(BWtempDilate,'noholes'); % Find edge of the subject

if opts.multiOutlines
    outlinesRaw = cell(length(BWedge)*2-1,1);
    a = 1;
        for ii = 1:length(BWedge)
        outlinesRaw{a,1} = BWedge{ii};
        outlinesRaw{a+1,1} = [NaN,NaN];
        a = a+2;
        end
    outline = fliplr((cell2mat(outlinesRaw)/scale)-(1/scale)); % final matrix containing the edge of multiple outlines
else
    outline = [(BWedge{1}(:,2)/scale)-(1/scale),(BWedge{1}(:,1)/scale)-(1/scale)]; % offset outline by the equivalent of 1px because bwboundaries starts the image at (x,y)=1 and outputs the outline at (x,y)=1.5. The outline should be set at the center of the poixel, thus (x,y)=0.5 when the origin of the BW mask is set to 0. Set the scale for the outline, in m. 
end

figure;
hold on
imagesc([0,(size(BWtempDilate,2)-1)/scale],[0,(size(BWtempDilate,1)-1)/scale],BWdiag); colormap('gray') % plot mask using BW image directly
plot(outline(:,1), outline(:,2), 'r','LineWidth',2);
axis equal

%% Export all the outlines
% opts.checkOutlines = 1; % option to check that the outlines look good
% path_outlines = fullfile(filepath,'outlines'); % outlines produced from BW masks (.csv)

fprintf('Exporting outlines...');

for i = 1:size(D_BW,1)
    BWtemp = flipud(importdata(quickfilepath(D_BW(i)))); % import mask
    BWdiag = bwmorph(BWtemp,'diag'); % diagonal fill
    BWtempDilate = bwmorph(BWdiag,'thicken',1); % dilate mask by 1px so the outline falls exactly on the edge of the mask
    BWedge = bwboundaries(BWtempDilate,'noholes'); %Find edge of the copepod

    if opts.multiOutlines
        outlinesRaw = cell(length(BWedge)*2-1,1);
        a = 1;
            for ii = 1:length(BWedge)
            outlinesRaw{a,1} = BWedge{ii};
            outlinesRaw{a+1,1} = [NaN,NaN];
            a = a+2;
            end
        outline = fliplr((cell2mat(outlinesRaw)/scale)-(1/scale)); % final matrix containing the edge of multiple outlines
    else
        outline = [(BWedge{1}(:,2)/scale)-(1/scale),(BWedge{1}(:,1)/scale)-(1/scale)]; % offset outline by the equivalent of 1px because bwboundaries starts the image at (x,y)=1 and outputs the outline at (x,y)=1.5. The outline should be set at the center of the poixel, thus (x,y)=0.5 when the origin of the BW mask is set to 0. Set the scale for the outline, in m. 
    end
    
%     BWoutline = BWedge{1}; % outline
%     outline = [(BWoutline(:,2)/scale)-(1/scale),(BWoutline(:,1)/scale)-(1/scale)]; % set the scale for the outline, in m
    filename = sprintf('iface_%05g', i); % Number sequentially
    writematrix(outline,fullfile(path_outlines,[filename,'.csv']))
    progressCount2(i,size(D_BW,1)); % display export progress
end

fprintf('done\n');

%% Convert raw velocity data in the proper format (without column headers) and units
D_PIVraw = dir(fullfile(path_PIVraw,'*.txt')); % set directory; accepted file type .txt

% Options
filestub = 'B';         % Default file name
opts.convertUnits = 1;  % Option to enable (1) or disable (0) unit conversion

fprintf('PIV files conversion...');
 for i=1:length(D_PIVraw)
     % Generate name of PIVclean files
     num = num2str(i,'%05d'); % create file number with 5 digits
     rf = strcat(filestub, num); % concatenate file number; filestub + file number
     
     % Import raw data
     A = importdata(quickfilepath(D_PIVraw(i))); % read data
     tempVelo = A.data; % velocity data

     % Convert units if necessary
     if opts.convertUnits == 1
     tempVelo(:,1:2) = tempVelo(:,1:2)/1000; % !!!! convert mm to m because DaVis exports x and y data in mm!!!!
     end

     % Export clean and converted data
     writematrix(tempVelo,fullfile(path_main,'PIVclean',[sprintf('B%05d',i),'.csv'])); % export clean data; accepted format .csv

     progressCount2(i,length(D_PIVraw)); % display export progress
end

fprintf('done\n');

%% Extract u and v from all the velocity files
D_velo = dir(fullfile(path_PIVclean,'*.csv')); % create the directory containing all the velocity files

% Initiate u and v matrices
VeloData = importdata(quickfilepath(D_velo(1)));
X = reshape(VeloData(:,1),length(unique(VeloData(:,1))),length(VeloData(:,1))/length(unique(VeloData(:,1))))';
Y = reshape(VeloData(:,2),length(unique(VeloData(:,1))),length(VeloData(:,1))/length(unique(VeloData(:,1))))';

uVelo = zeros(size(X,1),size(X,2),length(D_velo));
vVelo = zeros(size(X,1),size(X,2),length(D_velo));
UVspeed = zeros(size(X,1),size(X,2),length(D_velo));

fprintf('Velocity files extraction...');
    for i = 1:length(D_velo)
        tempVeloData = importdata(quickfilepath(D_velo(i)));
        uVelo(:,:,i) = reshape(tempVeloData(:,3),length(unique(tempVeloData(:,1))),length(tempVeloData(:,1))/length(unique(tempVeloData(:,1))))';
        vVelo(:,:,i) = reshape(tempVeloData(:,4),length(unique(tempVeloData(:,1))),length(tempVeloData(:,1))/length(unique(tempVeloData(:,1))))';
        UVspeed(:,:,i) = sqrt(uVelo(:,:,i).^2 + vVelo(:,:,i).^2); % compute flow velocity magnitude
    end
fprintf('done\n');

%% Set 0 masked values to NaN fo bulk flow subtraction
% uVelo(uVelo == 0) = NaN;
% vVelo(vVelo == 0) = NaN;
% UVspeed(UVspeed == 0) = NaN;
findIndex = uVelo(:,:,1); % store one velocity field to find indices for U and V = 0
uVelo([1:10 67:end],:,:) = NaN;
uVelo(:,[1:4 84:end],:) = NaN;
vVelo([1:10 67:end],:,:) = NaN;
vVelo(:,[1:4 84:end],:) = NaN;
%% Set scale from PIV data to scale outlines
scale = 7850.4; % px/m (applies to all the videos)

%% Calculate time-averaged u, v, and vorticity
% calculate mean U and V
meanU = mean(uVelo,3); % average performed across 3D matrix
meanV = mean(vVelo,3); % average performed across 3D matrix
meanSpeed = sqrt(meanU.^2 + meanV.^2);  % average performed across 3D matrix 

% Calculate mean vorticity
meanVortz=curl(X,Y,meanU,meanV);

%% Plot mean velocity field
% Plot the time-averaged u or v component to measure the mean bulk flow
opts.flowDirection = 1; % 1 = u (x direction); 2 = v (y direction); 3 = uv (actual direction)

figure; hold on
axis equal
vecScale = 0.1; % vector length scale
maxScale = 0.05; % max velocity magnitude

% Trace velocity profile across the field of view (run only if the flow speed is not known)
npoints = 200; % number of points along velocity profile to be extracted (for interpolation)

if opts.flowDirection == 1
   quiverColor(X,Y,meanU,meanV*0,maxScale,vecScale); 
   profile_velo = quickprofile(X,Y,meanU,npoints); % profile_u is the u velocity profile across X

elseif opts.flowDirection == 2
    quiverColor(X,Y,meanU*0,meanV,maxScale,vecScale);
    profile_velo = quickprofile(X,Y,meanV,npoints); % profile_v is the u velocity profile across X

elseif opts.flowDirection == 3
    quiverColor(X,Y,meanU,meanV,maxScale,vecScale); 
end


%% Calculate mean flow speed in the flume (run only if the flow speed is not known)
meanFlow = mean(profile_velo); % mean flow in the flume (oposite to the direction of swimming)
fprintf('Mean v flow speed: %f m.s-1\n',meanFlow)

%% mean V flow if flow rate is already known
meanFlow = 0; % mean V flow in m s-1

%% Subtract U from the original velociy profiles to reveal vortices (horizontal flow)
SubuVelo = uVelo-meanFlow;
fprintf('Mean flow subtracted...done\n')

%% Plot vorticity and velocity fields
% options
opts.plotVectors = 1;   % Option to overlay the velocity vectors
opts.plotMask = 0;      % Plot masks (instead of leaving the vorticity field empty)
opts.tightPlot = 1;     % Option to plot specific section of the velocity/vorticity fields

warning('off','MATLAB:polyshape:repairedBySimplify');   % disable unecessary polyshape warning
warning('off','MATLAB:polyshape:boundary3Points');      % disable unecessary polyshape warning

% Custom vorticity map
vorticity_colormap = cmocean('balance');  % define colormap
hLines = 200; % number of isolines

f = 1; % frame of interest

vortz = curl(X,Y,SubuVelo(:,:,f),vVelo(:,:,f)); % compute vorticity field for the frame of interest

% Increase the number of points in the vorticity field to ensure a smooth
% cutout of the field by the mask. The mask is made using outlines directly
% computed by the BW masks
smoothness = 0.0003; % parameter controlling the resolution of the vorticity field (smaller values = finer mesh)
[Xq,Yq] = meshgrid(min(X(:)):smoothness:max(X(:)),min(Y(:)):smoothness:max(Y(:))); % makes Y and Y data grid with higher resolution
Vq = interp2(X,Y,vortz,Xq,Yq); % interpolates vorticity (subVortz) data from original XY grid to higher resolution grid

% Fing the outlines from the BW masks
BWtemp = flipud(importdata(quickfilepath(D_BW(f)))); % import mask
BWdiag = bwmorph(BWtemp,'diag'); % diagonal fill
BWtempDilate = bwmorph(BWdiag,'thicken',1); % dilate mask by 1px so the outline falls exactly on the edge of the mask
B = bwboundaries(BWtempDilate,8);
outlinesRaw = cell(length(B)*2-1,1);
a = 1;
for ii = 1:length(B)
outlinesRaw{a,1} = B{ii};
outlinesRaw{a+1,1} = [NaN,NaN];
a = a+2;
end
        
outlines = cell2mat(outlinesRaw)/scale; % final matrix containing the outlines with holes

% Delete vorticity within the mask
[in] = inpolygon(Xq,Yq,outlines(:,2),outlines(:,1)); % finds XY data in high-res grid that falls in-on mask
in = (in-1).^2; % transforms 1 to 0 for conversion to NaN, and 0 to 1
in(in == 0) = NaN;
Vq = Vq.*in;

figure; hold on
I = importdata(quickfilepath(D_image(f))); % import corresponding image
imagesc(flipud(I),'XData',[min(X,[],'all') max(X,[],'all')],'YData',[min(Y,[],'all') max(Y,[],'all')]) % plot image
contourf(Xq,Yq,Vq,hLines,'edgecolor','none','FaceAlpha',0.75);   % Plot vorticity fields; change the opacity if needed (from 0 to 1)

pgon = polyshape(outlines(:,2),outlines(:,1));

if opts.plotMask == 1
plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','black')
end

% plot quiver
if opts.plotVectors
    % Remove vectors within the mask (often added by DaVis if the fill all
    % option is used or smoothing is enabled)
        in = inpolygon(X,Y,pgon.Vertices(:,1),pgon.Vertices(:,2));
        in = (in-1).^2; % transforms 1 to 0 for conversion to NaN, and 0 to 1
        in(in == 0) = NaN;
        uVeloMask = SubuVelo(:,:,f).*in; 
        vVeloMask = vVelo(:,:,f).*in; 
end

% Set axes length
axis equal
if opts.tightPlot == 1
    % Change parameters as needed
    axisXmin = 0.008; %was 0.040 (in m)
    axisYmin = 0.045; %was 0.005 (in m)
    axisXmax = 0.125; % was 0.05
    axisYmax = 0.12; % was 0.05

    axis([axisXmin axisXmax axisYmin axisYmax])
    axis off
else
    axisXmin = 0; %was 0.040 (in m)
    axisYmin = 0; %was 0.005 (in m)
    axisXmax = max(unique(X)); % was 0.05
    axisYmax = max(unique(Y)); % was 0.05

    xlim([axisXmin axisXmax])
    ylim([axisYmin axisYmax])
%     
    % Plot properties
    xlabel('distance (m)')
    ylabel('distance (m)')
end

    % Scale arrow for flow velocity
    originArrow = 0.005;
    qscale = 0.07;
    vecDensity = 1; % parameter to plot 1/n vector. 1 plots all the vectors
    p1 = [axisXmin+originArrow axisYmin+originArrow]; % First Point was 0.003
    p2 = [axisXmin+originArrow axisYmin+originArrow+0.1*qscale];% Second Point | scale arrow = 50 mm s-1
    dp = p2-p1; % Difference
    quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'r') 
    quiver(X(1:vecDensity:end,1:vecDensity:end),Y(1:vecDensity:end,1:vecDensity:end),qscale*uVeloMask(1:vecDensity:end,1:vecDensity:end),qscale*vVeloMask(1:vecDensity:end,1:vecDensity:end),0,'k','linewidth',0.5) % velocity vectors


% scalebar
scaleBarXmax = axisXmax-0.005; % was 0.003
scaleBarXmin = scaleBarXmax-0.02; % m scale (2cm) 
scaleBarY = axisYmin+0.005; % was 0.003
plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3)

% Colormap properties
colormap(vorticity_colormap)
hcb = colorbar;
maxVortMag = 8; % maximum vorticity magnitude value
clim([-maxVortMag maxVortMag]) % in s-1
hcb.Label.String = 'vorticity (s^{-1})';
hcb.Label.FontSize = 12;
hcb.Label.FontName = 'Arial';

formatFigure


%% Export all the frames
% Set export file path
filepath_svg = fullfile(path_main,'Figures','VortImOverlay-svg');
filepath_jpg = fullfile(path_main,'Figures','VortImOverlay-jpg');
mkdir(filepath_svg)
mkdir(filepath_jpg)


figure;
% Export loop
for ff=1:length(D_velo)

vortz = curl(X,Y,SubuVelo(:,:,ff),vVelo(:,:,ff)); % compute vorticity field for the frame of interest

% Increase the number of points in the vorticity field to ensure a smooth
% cutout of the field by the mask. The mask is made using outlines directly
% computed by the BW masks
% smoothness = 0.0003; % parameter controlling the resolution of the vorticity field (smaller values = finer mesh)
% [Xq,Yq] = meshgrid(min(X(:)):smoothness:max(X(:)),min(Y(:)):smoothness:max(Y(:))); % makes Y and Y data grid with higher resolution
Vq = interp2(X,Y,vortz,Xq,Yq); % interpolates vorticity (subVortz) data from original XY grid to higher resolution grid

% Fing the outlines from the BW masks
BWtemp = flipud(importdata(quickfilepath(D_BW(ff)))); % import mask
BWdiag = bwmorph(BWtemp,'diag'); % diagonal fill
BWtempDilate = bwmorph(BWdiag,'thicken',1); % dilate mask by 1px so the outline falls exactly on the edge of the mask
B = bwboundaries(BWtempDilate,8);
outlinesRaw = cell(length(B)*2-1,1);
a = 1;
for ii = 1:length(B)
outlinesRaw{a,1} = B{ii};
outlinesRaw{a+1,1} = [NaN,NaN];
a = a+2;
end
        
outlines = cell2mat(outlinesRaw)/scale; % final matrix containing the outlines with holes

% Delete vorticity within the mask
[in] = inpolygon(Xq,Yq,outlines(:,2),outlines(:,1)); % finds XY data in high-res grid that falls in-on mask
in = (in-1).^2; % transforms 1 to 0 for conversion to NaN, and 0 to 1
in(in == 0) = NaN;
Vq = Vq.*in;

hold on
I = importdata(quickfilepath(D_image(ff))); % import corresponding image
imagesc(flipud(I),'XData',[min(X,[],'all') max(X,[],'all')],'YData',[min(Y,[],'all') max(Y,[],'all')]) % plot image
contourf(Xq,Yq,Vq,hLines,'edgecolor','none','FaceAlpha',0.75);   % Plot vorticity fields; change the opacity if needed (from 0 to 1)

pgon = polyshape(outlines(:,2),outlines(:,1));
if opts.plotMask == 1
plot(pgon,'FaceColor','white','FaceAlpha',1,'EdgeColor','black')
end

% plot quiver
if opts.plotVectors
    % Remove vectors within the mask (often added by DaVis if the fill all
    % option is used or smoothing is enabled)

        in = inpolygon(X,Y,pgon.Vertices(:,1),pgon.Vertices(:,2));
        in = (in-1).^2; % transforms 1 to 0 for conversion to NaN, and 0 to 1
        in(in == 0) = NaN;
        uVeloMask = SubuVelo(:,:,ff).*in; 
        vVeloMask = vVelo(:,:,ff).*in; 
end

% Set axes length
axis equal
if opts.tightPlot == 1
%     % Change parameters as needed
%     axisXmin = 0.004; %was 0.040 (in m)
%     axisYmin = 0.04; %was 0.005 (in m)
%     axisXmax = 0.126; % was 0.05
%     axisYmax = 0.12; % was 0.05
% 
    axis([axisXmin axisXmax axisYmin axisYmax])
    axis off
else
%     axisXmin = 0; %was 0.040 (in m)
%     axisYmin = 0; %was 0.005 (in m)
%     axisXmax = max(unique(X)); % was 0.05
%     axisYmax = max(unique(Y)); % was 0.05

    xlim([axisXmin axisXmax])
    ylim([axisYmin axisYmax])
%     
    % Plot properties
    xlabel('distance (m)')
    ylabel('distance (m)')
end


    % Scale arrow for flow velocity
%     originArrow = 0.005;
%     qscale = 0.07;
%     vecDensity = 1; % parameter to plot 1/n vector. 1 plots all the vectors
%     p1 = [axisXmin+originArrow axisYmin+originArrow]; % First Point was 0.003
%     p2 = [axisXmin+originArrow axisYmin+originArrow+0.1*qscale];% Second Point | scale arrow = 50 mm s-1
%     dp = p2-p1; % Difference
    quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',0.7,'LineWidth', 2,'Color', 'r') 
    quiver(X(1:vecDensity:end,1:vecDensity:end),Y(1:vecDensity:end,1:vecDensity:end),qscale*uVeloMask(1:vecDensity:end,1:vecDensity:end),qscale*vVeloMask(1:vecDensity:end,1:vecDensity:end),0,'k','linewidth',0.5) % velocity vectors



% scalebar
% scaleBarXmax = axisXmax-0.005; % was 0.003
% scaleBarXmin = scaleBarXmax-0.02; % m scale (2cm) 
% scaleBarY = axisYmin+0.005; % was 0.003
plot([scaleBarXmin; scaleBarXmax], [scaleBarY; scaleBarY], '-k', 'LineWidth', 3)

% Colormap properties
colormap(vorticity_colormap)
hcb = colorbar;
% maxVortMag = 10; % maximum vorticity magnitude value
clim([-maxVortMag maxVortMag]) % in s-1
hcb.Label.String = 'vorticity (s^{-1})';
hcb.Label.FontSize = 12;
hcb.Label.FontName = 'Arial';

formatFigure

%Export paths
FilepathSVG = fullfile(filepath_svg,sprintf('vort%04d',ff));
FilepathJPG = fullfile(filepath_jpg,sprintf('vort%04d',ff));
% print('-depsc2', '-painters', FilepathSVG) %export eps
%print('-dpng', Filepath) %export png
print('-dsvg', '-vector',FilepathSVG) %exports to svg (preferred)
print('-djpeg', FilepathJPG)
%print('-dpdf', '-painters', '-bestfit', Filepath) %export pdf
%print('-dtiff', filepath) %export tif

disp(ff) 
hold off
clf
end

fprintf('Vorticity fields exported\n');
beep on; beep



