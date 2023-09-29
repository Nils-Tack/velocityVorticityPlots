function formatFigure
% figure properties
set(gca,'linewidth',1) % axes linewidth
set(gca, 'Layer', 'top') % move the axes to the top layer
set(groot, 'defaultAxesXColor', [0,0,0], ... % set axes color to black instead of charcoal grey
               'defaultAxesYColor', [0,0,0], ...
               'defaultAxesZColor', [0,0,0]);