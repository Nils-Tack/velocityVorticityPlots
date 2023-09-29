%% Plot pressure fields over actual image
% A: image data
% XX: meshgrid X data from velocity files
% YY: meshgrid Y data from velocity files
% subVortz: meshgrid vorticity data
% vel: velocity data
% vortMag: magnitude of scalebar for vorticity fields
% smoothness: smoothness parameter for masks (depends on data but 0.0002 to
% 0.0005 works well). <0.0001 it time consuming
% style: 1 = Classic, 2 = Classic + subtraction of low values, 3 = Low
% values black, 4 = low values black + subtraction of low values
% transparency: transparency paramneter of pressure field; for option 1&3,
% range = 100-150 with max = 255; for option 2&4, range =255-500.
% maskOn = uses outline to mask pressure field
% outline: import outline coordinates

function [Vq] = VorticityOverlayOnTif2(A,XX,YY,subVortz,vortMag,smoothness,transparency,opts,outline)
% opts.maskOn = maskOn;
% opts.style =   style;

[Xq,Yq] = meshgrid(min(XX(:)):smoothness:max(XX(:)),min(YY(:)):smoothness:max(YY(:))); % makes Y and Y data grid with higher resolution
Vq = interp2(XX,YY,subVortz,Xq,Yq); % interpolates vorticity (subVortz) data from original XY grid to higher resolution grid

if opts.maskOn
 [in] = inpolygon(Xq,Yq,outline(:,1),outline(:,2)); % finds XY data in high-res grid that falls in-on mask
 in = (in-1).^2; % transforms 1 to 0 for conversion to NaN, and 0 to 1
 in(in == 0) = NaN;
 Vq = Vq.*in; % all the data remains the same but values in masks are replaced by NaN
 end

if opts.style == 1||opts.style == 2
% red_blue=customColormapOriginal(99, [0 0 1], [1 0 0], [1 1 1]);% custom colormap for vorticity (blue, red, white)
red_blue=cmocean('balance');
else
red_blue=customColormapOriginal(99, [0 0 1], [1 0 0], [0 0 0]);% custom colormap for vorticity (blue, red, black)
end

imshow(A,'XData',[min(XX,[],'all') max(XX,[],'all')],'YData',[min(YY,[],'all') max(YY,[],'all')])
% imshow(A,'XData',[0 0.046179],'YData',[0 0.033191])
set(gca,'YDir','normal')
hold on;

% overlays vorticity field
hLines = 75;
[~, hContour]=contourf(Xq,Yq,Vq,hLines,'edgecolor','none');

colormap(red_blue)
hcb = colorbar;
caxis([-vortMag vortMag])
hcb.Label.String = 'vorticity (s^{-1})';
hcb.Label.FontSize = 12;
hcb.Label.FontName = 'Helvetica';

if opts.style == 1 || opts.style == 3
% transparency loop
drawnow;  % this is important, to ensure that FacePrims is ready in the next line!
hFills = hContour.FacePrims;  % array of TriangleStrip objects
[hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
for idx = 1 : numel(hFills)
   hFills(idx).ColorData(4) = transparency;   % default=255
end

elseif opts.style == 2 || opts.style == 4
drawnow;
% colorbar;
hFills = hContour.FacePrims;
[hFills.ColorType] = deal('truecoloralpha'); 
AlphaGradient=abs(linspace(-transparency,transparency,hLines+1));%The actual values that I want to use for the alpha values
% AlphaGradient = test;
for idx = 1:numel(hFills)
hFills(idx).ColorData(4) = AlphaGradient(idx);
end



% % transparency loop
% drawnow;  % this is important, to ensure that FacePrims is ready in the next line!
% hFills = hContour.FacePrims;  % array of TriangleStrip objects
% [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
% for idx = 1 : numel(hFills)
%    hFills(idx).ColorData(4) = transparency;   % default=255
% end
end