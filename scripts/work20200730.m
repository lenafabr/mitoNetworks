addpath('../../networktools')

%% load in mito network
load('~/grants/NIH_mitoaging/mitonetwork.mat')

%% interpolate network paths to be at roughly equivalent separations

dxwant = 0.1; % desired separation
NT.reinterpolateEdgePaths(dxwant);

%% get network graph for sharing
[NTgraph,netcoords] = makeGraphEdgePath(NT);
%%
options = struct();
options.edgeplotopt = {'MarkerSize',20};
NT.plotNetwork(options)

%% plot spheres around edgepath
path = NT.edgepath{1};
mitorad = 0.170;

[x0 y0 z0] = sphere(50);

hold all
for ec = 1:NT.nedge
    [ec NT.nedge]
    path = NT.edgepath{ec};
    for pc = [1:size(path,1)]
        x = x0*mitorad+path(pc,1); y = y0*mitorad + path(pc,2); z = z0*mitorad+path(pc,3);
        h = surfl(x, y, z);
        shading flat
        set(h, 'FaceAlpha', 0.3)
    end
end
hold off

%% define a trial curve

npt = 10;
xpt = (1:npt)'; 
ypt = 3*sin(2*pi/10*(1:npt))';
pts = [xpt ypt zeros(npt,1)]*0.2;

plot3(pts(:,1),pts(:,2),pts(:,3),'.-')
axis equal

%% 
mitorad = 0.17; % mitochondrial radius;
Nr = 10; % truncation order for radial gaussian quadrature

[rquadpts,rquadw]=lgwt(Nr,0,mitorad)

%lebedev quadrature points
lebedevfile = './lebedev/lebedev_011.txt';
lebdata = dlmread(lebedevfile);
phiquadpts = lebdata(:,1)*pi/180;
thquadpts = lebdata(:,2)*pi/180;
% weights for lebedev integration, so that integral of 1 over the surface
% gives 4 pi
lebquadw = lebdata(:,3)*4*pi;

%% test quadrature by computing sphere volume

% quantity we are integrating * weights
integweights = (rquadpts.^2.*rquadw)*lebquadw'; 

quadvol = sum(integweights(:))
Vtrue = 4/3*pi*mitorad.^3;

% fraction error
fracerr = (quadvol-Vtrue)/Vtrue

%% test quadrature by finding volume of spherical cap
% matrices of r and angles
Rmat = repmat(rquadpts,1,size(thquadpts,1));
thmat = repmat(thquadpts',size(rquadpts,1),1);
phimat = repmat(phiquadpts',size(rquadpts,1),1);
sphcoords = [Rmat(:),thmat(:),phimat(:)];
cartcoords = sphere2cart(sphcoords,[0 0 0]);

% only integrate above a certain cutoff
caph = 0.05;
ind = find(cartcoords(:,1)>mitorad-caph);

% quantity we are integrating * weights
integweights = (rquadpts.^2.*rquadw)*lebquadw'; 

quadvol = sum(integweights(ind))

Vtrue = pi/3*caph^2*(3*mitorad-caph)

% fraction error
fracerr = (quadvol-Vtrue)/Vtrue

% ---------
%% Test case: linear network
NT = NetworkObj();
NT.nodepos = [0 0 0; 1 0 0; 2 0 0; 3 0 0];
NT.edgenodes = [1 2; 2 3; 3 4];
NT.setupNetwork()
NT.interpolateEdgePaths(7);
NT.setCumEdgeLen();

options = struct();
options.edgeplotopt = {'MarkerSize',20};
NT.plotNetwork(options)


%% get quadrature points for each sphere
mitorad = 0.17; % mitochondrial radius;
Nr = 50; % truncation order for radial gaussian quadrature

[rquadpts,rquadw]=lgwt(Nr,0,mitorad);

%lebedev quadrature points
lebedevfile = './lebedev/lebedev_083.txt';
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

length(quadweights)


%% Place quadrature points within spheres all along the network, only considering neighbor overlap

% matrices of r and angles
Rmat = repmat(rquadpts,1,size(thquadpts,1));
thmat = repmat(thquadpts',size(rquadpts,1),1);
phimat = repmat(phiquadpts',size(rquadpts,1),1);
sphquadcoords = [Rmat(:),thmat(:),phimat(:)];

donenode = false(NT.nnode,1);
allquadpts = []; allquadw = [];

for ec = 1:NT.nedge
    [ec NT.nedge]
    n1 = NT.edgenodes(ec,1);
    n2 = NT.edgenodes(ec,2);
    path = NT.edgepath{ec};
    
    if (~donenode(n1)) % do the first node of the edge
        pts = sphere2cart(sphquadcoords,path(1,:));
        curweights = quadweights;
        
        % remove all points intersecting with neighbor sphere on other
        % edges
        for ec2 = NT.nodeedges(n1,1:NT.degrees(n1))
            if (NT.edgenodes(ec2,2) == n1) % incoming edge
                neighb = NT.edgepath{ec2}(end-1,:);
                diffs = pts - neighb;
                dists = sum(diffs.^2,2);
                ind = dists>mitorad.^2; % keep only non-overlapping points
                pts = pts(ind,:);
                curweights = curweights(ind,:);
            end
        end
        allquadpts = [allquadpts; pts];
        allquadw = [allquadw; curweights];
        
        donenode(n1) = true;
    end
    
    % do remaining nodes on edge
    for pc = 2:size(path,1)
        pts = sphere2cart(sphquadcoords,path(pc,:));
        
        
        if (pc==size(path,1))
            if(donenode(n2))
                continue % last node already taken care of
            else
                donenode(n2) = true;
            end
        end
        
        % remove points overlapping with previous sphere
        diffs = pts - path(pc-1,:);
        dists = sum(diffs.^2,2);
        ind = dists>mitorad.^2;        
        allquadpts = [allquadpts; pts(ind,:)];
        allquadw = [allquadw; quadweights(ind)];
        
    end
end

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

%%
options = struct();
options.edgeplotopt = {'MarkerSize',20};
NT.plotNetwork(options)
%
hold all
plot3(allquadpts(:,1),allquadpts(:,2),allquadpts(:,3),'r.')
hold off

%% estimate total volume
Vtot = sum(allquadw)

L = norm(NT.nodepos(end,:)-NT.nodepos(1,:))

Vtrue = pi*mitorad^2*L + 4/3*pi*mitorad^3

fracerror = (Vtot-Vtrue)/Vtrue

%% set up cubical voxels
%vxsize = 0.1083; % voxel size (in um)
%vxsize = 0.1625; % voxel size (in um)
vxsize = 0.0433

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

% check array accumulation
check = [23 20 25];
ind = find(vxind(:,1)==check(1) & vxind(:,2)==check(2) & vxind(:,3) == check(3)); 
sum(allquadw(ind))
mitoimg(check(1),check(2),check(3))

%% plot a slice through the image
imshow(mitoimg(:,:,56),[])

%% plot a maximum intensity projection
maxprojimg = max(mitoimg,[],3);

imshow(maxprojimg,[])

%% 