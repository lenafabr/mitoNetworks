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
options.plotedgepath = true;
NT.plotNetwork(options)

%% get quadrature points for each sphere
mitorad = 0.17; % mitochondrial radius;
Nr = 40; % truncation order for radial gaussian quadrature

[rquadpts,rquadw]=lgwt(Nr,0,mitorad);

%lebedev quadrature points
lebedevfile = './lebedev/lebedev_031.txt';
lebdata = dlmread(lebedevfile);
phiquadpts = lebdata(:,1)*pi/180;
thquadpts = lebdata(:,2)*pi/180;
% weights for lebedev integration, so that integral of 1 over the surface
% gives 4 pi
lebquadw = lebdata(:,3)*4*pi;

% overall coords for all integration points
Rmat = repmat(rquadpts,1,size(thquadpts,1));
thmat = repmat(thquadpts',size(rquadpts,1),1);
phimat = repmat(phiquadpts',size(rquadpts,1),1);
sphquadcoords = [Rmat(:),thmat(:),phimat(:)];

% overall weights for all integration points
quadweights = (rquadpts.^2.*rquadw)*lebquadw'; 
quadweights = quadweights(:);

nquad = length(quadweights)

%% get matrix of pairwise distances between all nodes along edgepaths
npath = size(netcoords,1);
Adist = zeros(npath,npath);
for pc1 = 1:npath
    diffs = netcoords - netcoords(pc1,:);
    dists = sqrt(sum(diffs.^2,2));
    Adist(pc1,:) = dists;
end

%% Place quadrature points within spheres all along the network
% with global search for overlap

% matrices of r and angles
Rmat = repmat(rquadpts,1,size(thquadpts,1));
thmat = repmat(thquadpts',size(rquadpts,1),1);
phimat = repmat(phiquadpts',size(rquadpts,1),1);
sphquadcoords = [Rmat(:),thmat(:),phimat(:)];

allquadpts = []; allquadw = [];

for pc = 1:npath
    if (mod(pc,10)==0); disp([pc npath]); end
    cent = netcoords(pc,:);
    
    pts = sphere2cart(sphquadcoords,cent);
    curweights = quadweights;
    
    checknodes = find(Adist(pc,1:pc-1)<2*mitorad);
    for pc2 = checknodes % check for overlap with all prior spheres
        %if (Adist(pc,pc2)>2*mitorad); continue; end % no overlap possible
        diffs = pts - netcoords(pc2,:);
        dists = sum(diffs.^2,2);
        ind = dists>mitorad.^2; % keep only non-overlapping points
        pts = pts(ind,:);
        curweights = curweights(ind,:);
    end
    allquadpts = [allquadpts; pts];
    allquadw = [allquadw; curweights];
end

%% slow plotting of all quadrature points
options = struct();
options.edgeplotopt = {'MarkerSize',20};
NT.plotNetwork(options)
%
hold all
plot3(allquadpts(:,1),allquadpts(:,2),allquadpts(:,3),'r.')
hold off

%% set up cubical voxels
% voxel size (in um)
%vxsize = 0.1625; 
%vxsize = 0.1083; 
%vxsize = 0.065; 
vxsize = 0.0433;

buffer = 1; % length of approximate buffer zone on each side (in um)

minvals = min(allquadpts);
maxvals = max(allquadpts);
rangevals = maxvals-minvals;

viewboxmin = minvals - buffer;
% number of voxels
nvxls = ceil((rangevals + 2*buffer)/vxsize);
nvxl = max(nvxls)

%% map each quadrature point to a particular voxel index
vxind = floor((allquadpts - viewboxmin)/vxsize)+1;

% sum of quadrature weights for all points that fall in a pixel
mitoimg = accumarray(vxind,allquadw,[nvxl,nvxl,nvxl]);
%% plot a slice through the image
imshow(mitoimg(:,:,50),[])
set(gcf,'Position',[100 100 500 500])
set(gca,'Position',[0.05 0.05 0.9 0.9])
%% plot a maximum intensity projection
maxprojimg = max(mitoimg,[],3);

imshow(maxprojimg,[])
set(gcf,'Position',[100 100 500 500])
set(gca,'Position',[0.05 0.05 0.9 0.9])

%% save results
filename = sprintf('../results/mitonetwork_vx%g.mat',vxsize*1000)
save(filename,'mitoimg','NTgraph','netcoords','vxsize','mitorad','nquad','dxpath')

% -------------------------------------------------
%% check saved results
load('../results/mitonetwork_vx43.3.mat')

%%
maxprojimg = max(mitoimg,[],3);

imshow(maxprojimg,[])
%%
plot(NTgraph,'XData',netcoords(:,1),'YData',netcoords(:,2),'ZData',netcoords(:,3))