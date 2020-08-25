% directory name where files containing network structure are stored
dirname = '/data/proj/mitochondrialNetworks/Viana2020MendeleyDataset/';

% name for a specific network
% you will need these files:
% name_skeleton.txt (containing skeleton coords)
% name.coo (containing node coords)
% name.gnet (containing connectivity)
name = '020614_1_4_005'; 

%% load skeleton
data = readtable([dirname sprintf('%s_skeleton.txt',name)]);

skeldata = [data.Points_0 data.Points_1 data.Points_2]

%% plot skeleton
plot3(skeldata(:,1),skeldata(:,2),skeldata(:,3),'r.')
axis equal

%% get network node and connectivity information, save as network object
nodedata = dlmread([dirname sprintf('%s.coo',name)]);
edgedata = dlmread([dirname sprintf('%s.gnet',name)]);
edgenodes = edgedata(2:end,1:2)+1;

NT = NetworkObj();
NT.nodepos = nodedata;
NT.edgenodes = edgenodes
NT.setupNetwork();
NT.Name = name;

%% map skeleton to network edges
mapSkeleton2Network(NT,skeldata)

%% plot skeleton with edge paths
NT.plotNetwork(struct('plotedgepath',true))

hold all
plot3(skeldata(:,1),skeldata(:,2),skeldata(:,3),'r.')
axis equal
hold off

%%
% voxelate the network for several different voxel sizes
allmitoimgs = voxelateMitoNetwork(NT,[0.043,0.065,0.1083,0.1625],options)