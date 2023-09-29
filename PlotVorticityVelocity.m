function PlotVorticityVelocity(mask,X,Y,mag,vort,subU,subV,outline,frame)
opts.arrow = 0;
opts.vectors = 1;
opts.maskOn = mask;

axisMagnitude = mag;

% Plot the vorticity fields with a custom colormap
% purp_orng=customColormapOriginal(99, [0 0 1], [1 0 0], [1 1 1]);% uses customColormap function

% if opts.style == 1||opts.style == 2
% red_blue=customColormapOriginal(99, [0 0 1], [1 0 0], [1 1 1]);% custom colormap for vorticity (blue, red, white)
red_blue=cmocean('balance');
% else
% red_blue=customColormapOriginal(99, [0 0 1], [1 0 0], [0 0 0]);% custom colormap for vorticity (blue, red, black)
% end
 
figure; hold on
% dataVortzz = cell2mat(dataVortz(ff-1));
contourf(X,Y,vort,75,'edgecolor','none') %was 75


caxis([-axisMagnitude axisMagnitude])
colormap(red_blue)%custom colormap

hcb = colorbar;
hcb.Label.String = 'vorticity (s^{-1})';
hcb.Label.FontSize = 12;
hcb.Label.FontName = 'Helvetica';
axis equal 
axis off
title(sprintf('Frame: %d\n',frame))

qscale=0.02; % scale for vector arrows

% Scale arrow for flow velocity
if opts.arrow
p1 = [axisXmin+0.003 axisYmin+0.003]; % First Point
p2 = [axisXmin+0.003+1*qscale axisYmin+0.003];% Second Point (m s-1)
dp = p2-p1;% Difference
quiver(p1(1),p1(2),dp(1),dp(2),0,'MaxHeadSize',0.3,'LineWidth', 2,'Color', 'r') 
end

set(gcf, 'Position', get(0, 'Screensize'));
set(gca,'Color','b') % set color to black

if opts.vectors % plots velocity vectors on top of vorticity fields
    quiver(X,Y,qscale*subU,qscale*subV,0,'k')
end

if opts.maskOn
%     patch(outline(:,1)/1000,outline(:,2)/1000,'k')
end
end