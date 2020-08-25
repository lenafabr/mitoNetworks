function allmitoimgs = voxelateMitoNetwork(NT,vxsizes,options)
% starting with a network, object convert it into voxels of a particular
% brightness
% vxsizes = list of desired cubical voxel sizes, in um
% returns a cell list, each entry is a 3D image for one of the voxel sizes

% default options
opt = struct();
opt.dxpath = 0.1; % point separation for tracing out edges
opt.mitorad = 0.17; % assumed mitochondrial radius (in um)
opt.Nr = 30; % truncation order for radial gaussian quadrature
% file containing lebedev quadrature
opt.lebedevfile = './lebedev/lebedev_031.txt';

opt.buffer = 1; % size of border buffer in the image, in um

opt.dodisplay = 1; % how much displaying to do

if (exist('options','var'))
    opt = copyStruct(options,opt);
end

%% interpolate network paths to be at roughly equivalent separations
NT.reinterpolateEdgePaths(opt.dxpath);

%% get network graph for sharing
[NTgraph,netcoords] = makeGraphEdgePath(NT);
%%
if (opt.dodisplay>1)
    dispoptions = struct();
    dispoptions.edgeplotopt = {'MarkerSize',20};
    dispoptions.plotedgepath = true;
    NT.plotNetwork(dispoptions)
end
%% get quadrature points for each sphere

[rquadpts,rquadw]=lgwt(opt.Nr,0,opt.mitorad);

%lebedev quadrature points
lebdata = dlmread(opt.lebedevfile);
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
    
    checknodes = find(Adist(pc,1:pc-1)<2*opt.mitorad);
    for pc2 = checknodes % check for overlap with all prior spheres
        %if (Adist(pc,pc2)>2*mitorad); continue; end % no overlap possible
        diffs = pts - netcoords(pc2,:);
        dists = sum(diffs.^2,2);
        ind = dists>opt.mitorad.^2; % keep only non-overlapping points
        pts = pts(ind,:);
        curweights = curweights(ind,:);
    end
    allquadpts = [allquadpts; pts];
    allquadw = [allquadw; curweights];
end

%% slow plotting of all quadrature points
if (opt.dodisplay>2)
    dispoptions = struct();
    dispoptions.edgeplotopt = {'MarkerSize',20};
    NT.plotNetwork(dispoptions)
    %
    hold all
    plot3(allquadpts(:,1),allquadpts(:,2),allquadpts(:,3),'r.')
    hold off
end
%% set up cubical voxels
for vxc = 1:length(vxsizes)
    % voxelate for a given voxel size
    vxsize = vxsizes(vxc);
    
    minvals = min(allquadpts);
    maxvals = max(allquadpts);
    rangevals = maxvals-minvals;
    
    viewboxmin = minvals - opt.buffer;
    % number of voxels
    nvxls = ceil((rangevals + 2*opt.buffer)/vxsize);
    nvxl = max(nvxls)
    
    %% map each quadrature point to a particular voxel index
    vxind = floor((allquadpts - viewboxmin)/vxsize)+1;
    
    % sum of quadrature weights for all points that fall in a pixel
    mitoimg = accumarray(vxind,allquadw,[nvxl,nvxl,nvxl]);
    
    %% plot a slice through the image
    %imshow(mitoimg(:,:,50),[])
    %set(gcf,'Position',[100 100 500 500])
    %set(gca,'Position',[0.05 0.05 0.9 0.9])
    %% plot a maximum intensity projection
    if (opt.dodisplay>0)
        figure(vxc)
        maxprojimg = max(mitoimg,[],3);
        
        imshow(maxprojimg,[])
        set(gcf,'Position',[100 100 500 500])
        set(gca,'Position',[0.05 0.05 0.9 0.9])
        
        title(sprintf('%s, voxels %g um',NT.Name,vxsize),'Interpreter','none')
    end
    
    allmitoimgs{vxc} = mitoimg;
    
end
