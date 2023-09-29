function [profile_u,profile_r,profile_x,profile_y,xp,yp] = quickprofile(X,Y,U,npoints,xp,yp,ylab)
% Create a quick velocity profile along a line in a velocity field.
% If points are not defined, get them from a plot using ginput.

%% Options
plotOn = 1;

%%
if nargin < 5
    [xp,yp] = ginput(2);
    hold on
    plot(xp(1),yp(1),'go')
    plot(xp(2),yp(2),'ro')
    plot(xp,yp,'k--')
end

profile_x = linspace(xp(1),xp(2),npoints)'; % x coordinates of the interpolated points along the line
profile_y = linspace(yp(1),yp(2),npoints)'; % y coordinates of the interpolated points along the line
profile_r = sqrt((profile_x - profile_x(1)).^2 +...
                 (profile_y - profile_y(1)).^2); % enables plotting alng one line
profile_u = interp2(X,Y,U,profile_x,profile_y); % interpolate along the line

%%
if plotOn == 1
    figure
    plot(profile_r,profile_u,'.')
    xlabel('length of profile (m)')
    if nargin < 7
        ylabel('velocity (m s^{-1})')
    else
        ylabel(ylab)
    end
end
