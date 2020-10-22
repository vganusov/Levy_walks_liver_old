function varargout = labelPlotMSD(obj, ha)
%%LABELPLOTMSD A convenience method to set the axes labels.
%
% obj.labelPlotMSD(ha) sets the axis label of the axes with
% the specified handle ha. It is meant for axes containing the
% plot of the mean-square-displacement.
%
% hl = obj.plotMSD(...) returns the handle to the generated
% labels.

if nargin < 2
    ha = gca;
end

hl = NaN(2, 1);
hl(1) = xlabel(ha, ['Log_{10} Delay (' obj.time_units ')']);
hl(2) = ylabel(ha, ['Log_{10} MSD (' obj.space_units '^2)']);

xl = xlim(ha);
%xlim(ha, [0 xl(2)]);
xlim(ha, [0.5, 4]);
yl = ylim(ha);
%ylim(ha, [0 yl(2)]);
ylim(ha, [0, 4.5]);
box(ha, 'off')

if nargout > 0
    varargout{1} = hl;
end
end