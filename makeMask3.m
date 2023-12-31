function BW5 = makeMask3(I, opts)

%  Enhance contrast and preserve edges (may work in specific cases)
  edgeThreshold = 0.4; 
  amount = 0.5; % the closer to 1, the higher the contrast
   
    % vertical correction
    % I(:,200:500) = localcontrast(I(:,200:500),edgeThreshold, amount);% excellent to enhance contrast and preserve edges
    
    % horizontal correction
       %I(250:650,:) = localcontrast(I(250:650,:),edgeThreshold, amount);
       %I(1250:1800,:) = localcontrast(I(1250:1800,:),edgeThreshold, amount);
         % I(1500:1900,:) = localcontrast(I(1500:1900,:),edgeThreshold, amount);
    
    
% Adjust contrast and gamma correction on original image
I2 = imadjust(I,[0.30 1],[],0.9); % change values accordingly. (Image, [low_in high_in], [], gamma); If gamma is less than 1, then imadjust weights the mapping toward higher (brighter) output values.

% adaptive binarization
%BW2 = imbinarize(I2,'adaptive','Sensitivity',0.45); % Sensitivity factor for adaptive thresholding; higher = more details, lower = less details
% Number in the range [0, 1]. A high sensitivity value leads to thresholding more pixels as foreground, at the risk of including some background pixels.
BW2 = imbinarize(I2, 0.4);

% Mask unwanted horizontal and/or vertical sections of the images
    % horizontal
    BW2(1:100,:) = 0;
    BW2(size(BW2,1)-200:size(BW2,1),:) = 0;

    % vertical
    %BW2(:,1:600) = 0;                       % left
    %BW2(:,size(BW2,1)-400:size(BW2,1)) = 0; % right

% Clean mask
% Extract the largest blob only
%BW3 = bwareafilt(BW2, 1); %  keeps the n largest objects; change based on the number of masks needed in one image
%-> alternatively, select only blobs of a certain size.
props = regionprops(BW2, 'Area'); % Find areas of each blob and put into a structure.
allAreas = [props.Area];
BW3 = bwareafilt(BW2, [100,max(allAreas)]); % select all ther blob ranging from the largest to a specific, smaller threshold

% Fill holes
% BW3 = imfill(BW3, 'holes');
BW3 = ~bwareaopen(~BW3, 10); % awesome to only fill small holes. Perfect for complicated shapes and systems like multiple legs/appendages

% If anterior section of the fish is too dark
% se = strel('disk',40);
% BW3(:,1:800) = imclose(BW3(:,1:800),se);

% Smooth edges
% Dilate the Image to preserve small details (like the tip of the tail)
se90 = strel('line',10,90);
se0 = strel('line',10,0);

% Thicken the mask
BW4 = bwmorph(BW3,'thicken',5); % was 10

% Blur the image
windowWidth = 15; % was 25
BWblurr = conv2(double(BW4), ones(windowWidth)/windowWidth^2, 'same');

% Threshold again.
BWblurr = BWblurr > 0.5;

% Shrink mask
BW5 = bwmorph(BWblurr,'shrink',2.5); % 1/2 of original dilation since blur was cut to 1/2 already

% Plot combination of original, blurred, and final mask images
if opts.plotComparison
f = figure;
f.Position(1:2) = [0 0];
subplot(2,2,1)
imshow(I)
title('Original image')
subplot(2,2,2)
imshow(I2)
title('Gamma correction')
subplot(2,2,3)
imshow(BW2)
title('Adaptive binarization + removing edges')
subplot (2,2,4)
imshow(BW5);
title('Smoothed mask')
end

end