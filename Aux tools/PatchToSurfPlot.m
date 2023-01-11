function [X,Y,Z,Shear,FibDev] = PatchToSurfPlot(P,Grid,PlotVar)
% This function takes the P array which is set up to plot the individual
% draping cells as patches using the 'patch' function and converts it into
% grid vectors X,Y,Z. These are plotted using the 'surf' function.

X = reshape(mean(P(:,:,1),2),Grid-1);
Y = reshape(mean(P(:,:,2),2),Grid-1);
Z = reshape(mean(P(:,:,3),2),Grid-1);
Shear = reshape(mean(P(:,:,4),2),Grid-1);
FibDev = reshape(mean(P(:,:,5),2),Grid-1);

if strcmpi(PlotVar,'Shear')
    surf(X,Y,Z,Shear,'edgecolor','none')
elseif strcmpi(PlotVar,'FibDev')
    surf(X,Y,Z,FibDev,'edgecolor','none')
end
colormap jet
axis equal
colorbar

xlabel('x')
ylabel('y')
zlabel('z')

end