addpath('../../contractingER/scripts')

%% load in mito network
load('~/UCSD/grants/NIH_mitoaging/mitonetwork.mat')

%%
options = struct();
options.edgepath = NT.edgepath;
options.edgeplotopt = {'MarkerSize',10};
plotNetwork(NT.nodepos,NT.edgenodes,options)

%% plot spheres around edgepath
path = NT.edgepath{1};
cylrad = 0.170;

[x0 y0 z0] = sphere(100);

hold all
for ec = 1:NT.nedge
    [ec NT.nedge]
    path = NT.edgepath{ec};
    for pc = [1:size(path,1)]
        x = x0*cylrad+path(pc,1); y = y0*cylrad + path(pc,2); z = z0*cylrad+path(pc,3);
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
lebedevfile = './lebedev/lebedev_009.txt';
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
%% convert quadratures around each bead to cartesian space
