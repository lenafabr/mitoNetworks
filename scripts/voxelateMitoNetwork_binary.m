%% load in mito network
load('~/grants/NIH_mitoaging/mitonetwork.mat')

%% interpolate network paths to be at roughly equivalent separations

dxpath = 0.1; % desired separation
NT.reinterpolateEdgePaths(dxpath);

%% get network graph for sharing
[NTgraph,netcoords] = makeGraphEdgePath(NT);
%%
options = struct();
options.edgeplotopt = {'MarkerSize',20};
NT.plotNetwork(options)


%% set up cubical voxels
% voxel size (in um)
vxsize = 0.1625; 
%vxsize = 0.1083; 
%vxsize = 0.065; 
%vxsize = 0.0433;

buffer = 1; % length of approximate buffer zone on each side (in um)

minvals = min(netcoords);
maxvals = max(netcoords);
rangevals = maxvals-minvals;

viewboxmin = minvals - buffer-mitorad;
% number of voxels
nvxls = ceil((rangevals + 2*buffer)/vxsize);
nvxl = max(nvxls)

% get voxel point coordinates

vxvals = viewboxmin(1) + vxsize*((1:nvxl)+1/2);
vyvals = viewboxmin(2) + vxsize*((1:nvxl)+1/2);
vzvals = viewboxmin(3) + vxsize*((1:nvxl)+1/2);

[VX,VY,VZ] = meshgrid(vxvals,vyvals,vzvals);

Vcoords = [VX(:) VY(:) VZ(:)];

%% for each coordinate, decide if its close enough to network points
mitoimgbinary = zeros(nvxl,nvxl,nvxl);
for xc = 1:nvxl
    if (mod(xc,10)==0); disp([xc nvxl]); end
    for yc = 1:nvxl
        for zc = 1:nvxl
            diffs = [vxvals(xc) vyvals(yc) vzvals(zc)] - netcoords;
            dists = sum(diffs.^2,2);
            if (min(dists)<mitorad^2)
                mitoimgbinary(xc,yc,zc) = 1;
            end
        end
    end
end

%%

%% plot a slice through the image
imshow(mitoimgbinary(:,:,70),[])
set(gcf,'Position',[100 100 500 500])
set(gca,'Position',[0.05 0.05 0.9 0.9])
%% plot a maximum intensity projection
maxprojimg = max(mitoimgbinary,[],3);

imshow(maxprojimg,[])
set(gcf,'Position',[100 100 500 500])
set(gca,'Position',[0.05 0.05 0.9 0.9])

%% save results
filename = sprintf('../results/mitonetwork_binary_vx%g.mat',vxsize*1000)
save(filename,'mitoimgbinary','NTgraph','netcoords','vxsize','mitorad','nquad','dxpath')

%% check saved results
load('../results/mitonetwork_binary_vx43.3.mat')

%%
maxprojimg = max(mitoimgbinary,[],3);

imshow(maxprojimg,[])
%%
plot(NTgraph,'XData',netcoords(:,1),'YData',netcoords(:,2),'ZData',netcoords(:,3))