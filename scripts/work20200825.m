dirname = '/data/proj/mitochondrialNetworks/Viana2020MendeleyDataset/';
name = '020614_1_9_155'; % specific image name

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

%% plot network and skeleton
NT.plotNetwork()
hold all
plot3(skeldata(:,1),skeldata(:,2),skeldata(:,3),'r.')
axis equal

%% map skeleton to network edges
mapSkeleton2Network(NT,skeldata)

%% plot skeleton with edge paths
NT.plotNetwork(struct('plotedgepath',true))

hold all
plot3(skeldata(:,1),skeldata(:,2),skeldata(:,3),'r.')
axis equal

%%
% try voxelating
allmitoimgs = voxelateMitoNetwork(NT,[0.043,0.065],options)

% ---------------------------------
%% --------- old scripts to parse what aidan had already processed ---------
%% pull up a network structure
NT = NetworkObj('~/proj/networkDiffusion/mitochondria/networkFiles_WT_withLoops_2019Mar8/RFP-100914_244_99_899.net');
%%
%NT = NetworkObj('~/proj/networkDiffusion/mitochondria/fullskeletons/network_wt2.out');
%%
plotNetwork(NT.nodepos,NT.edgenodes)
axis equal


%% load skeleton
skeldata= dlmread('~/proj/networkDiffusion/mitochondria/WT_skeletons_LCC/RFP-100914_244_95_887.ske');

%% plot network and skeleton

%plotNetwork(NT.nodepos,NT.edgenodes)
%axis equal
%hold all
plot3(skeldata(:,1),skeldata(:,2),skeldata(:,3),'r.')
axis equal
hold off

%% map skeleton to edge paths
mapSkeleton2Network(NT,skeldata)

%% plot with edge paths
plotNetwork(NT.nodepos,NT.edgenodes)
axis equal
hold all
%plot3(skeldata(:,1),skeldata(:,2),skeldata(:,3),'r.')
hold off

cmat = lines(NT.nedge);
hold all
for ec = 1:NT.nedge
    edgepath = NT.edgepath{ec};
    plot3(edgepath(:,1),edgepath(:,2),edgepath(:,3),'.-','Color',cmat(ec,:))
end
hold off