function mapSkeleton2Network(NT,skeldata,options)
% given a network object
% and a list of datapoints that go along the non-straight backbone of the
% structure (including nodes, points are out of order)
% trace out paths for each of the individual edges
% and save them in NT.edgepath

% default parameters
opt = struct();
% relative weighting of distance and angle in deciding which point comes
% next
opt.weightdist = 10;
opt.weightang = 1;
opt.rhocutoff = 0.3; % minimal cosine(angle) for consecutive nodes
opt.distcutoffscl= 8; % max distance allowed, as a factor of shortest step between consecutive datapoints

% display edgepaths
opt.dodisplay = 0;

if (nargin>2)
    opt =copyStruct(options,opt,1);
end

%% go through skeleton, find consecutive stretches of points that lead between 2 nodes
if (opt.dodisplay)
    NT.plotNetwork()
    axis equal
    hold all
    plot3(skeldata(:,1),skeldata(:,2),skeldata(:,3),'r.')
end

% do assume that the network connectivity is correctly known
% trace according to a weighted utility function involving node distance
% and angle

% pick out data points that match nodes
whichnode = zeros(size(skeldata,1),1);
for nc = 1:NT.nnode
    dists = sqrt(sum((skeldata-NT.nodepos(nc,:)).^2,2));
    [mindist,minind] = min(dists);
    if mindist>0.1 
        error('failed to find data point matching this node')
    end    
    whichnode(minind) = nc;
    nodeinds(nc) = minind;
end

%% 
% which data points have been assigned already
done = false(size(skeldata,1),1);
%done(whichnode>0) = true;
% how many undone edges for each node
degreeleft= NT.degrees;

alldiffs = diff(skeldata,1);
allstepsizes = sqrt(sum(alldiffs.^2,2));
distcutoff = min(allstepsizes)*opt.distcutoffscl;
allmin = min(allstepsizes);

cmat = jet(NT.nedge);

nedgedone = 0;
doneedge = zeros(NT.nedge,1);

for nc = 1:NT.nnode
    nodepos = NT.nodepos(nc,:);
    edgepath = nodepos;
    while degreeleft(nc)>0
        % find the nearest unprocessed datapoint to this node                        
        undoneind = find(~done);
        
        dists = sqrt(sum((skeldata(undoneind,:)-nodepos).^2,2));
        dists(dists<1e-3) = inf; % ignore 0 distance
        
        % do not allow direct node-node connection      
        dists(whichnode(undoneind)>0) = inf;
       
        
        [mindist,minind] = min(dists);
        minind = undoneind(minind);
        
        if (mindist>distcutoff)
            error('failed to find sufficiently close node!')
        end
        edgepath = [nodepos; skeldata(minind,:)];        
        
        ct = 2;
        curdone = false(size(skeldata,1),1);
        curdone(minind) = true;    
        curdone(nodeinds(nc)) = true;
        while 1     
            laststep = edgepath(end,:)-edgepath(end-1,:);
            
            % find the nearest next node among those in the right direction
            diffs = skeldata(undoneind,:) - edgepath(end,:);
            dists = sqrt(sum(diffs.^2,2));
            rhos = sum(diffs.*laststep,2)./dists./norm(laststep);  
            
            possible = find(rhos>opt.rhocutoff & dists < distcutoff);
            if (isempty(possible)) % relax criterion
                error('no possible next point')
            %    possible = find(rhos>rhocutoff2);
            end                     
            
            utility = -opt.weightang*rhos(possible) + opt.weightdist*(dists(possible)-allmin).^2;                       
            
           % [mindist,minind] = min(dists(possible));
           [minutil,minind] = min(utility);
            nextind = undoneind(possible(minind));
            
            ct = ct+1;
            
            edgepath(end+1,:) = skeldata(nextind,:);
            curdone(nextind) = true;            
                        
            n2 = whichnode(nextind);
           
            if (n2>0)
                % reached the end of the edge   
                
                if (opt.dodisplay)
                    plot3(edgepath(:,1),edgepath(:,2),edgepath(:,3),'o','LineWidth',2,'Color',cmat(nedgedone+1,:))
                end
                
                foundedge = false; % successfully mapped an edge
                for cc = 1:NT.degrees(nc)
                    ec = NT.nodeedges(nc,cc);
                    if (doneedge(ec)); continue; end
                    if (NT.edgenodes(ec,1)==nc & NT.edgenodes(ec,2)==n2)
                        NT.edgepath{ec} = edgepath;
                        foundedge = true;                         
                        break;
                    elseif (NT.edgenodes(ec,2)==nc & NT.edgenodes(ec,1)==n2)
                        NT.edgepath{ec} = flipud(edgepath);
                        foundedge = true;
                        break;
                    end
                end                
                
                if (foundedge) % done with this edge
                    
%                     if (doneedge(ec))
%                         % this edge has already been done!
%                         disp(sprintf('Adding extra edge to network. %d %d',nc,n2))
%                         NT.edgenodes(end+1,:) = [nc n2];
%                         newec = size(NT.edgenodes,1);
%                         NT.edgepath{newec} = edgepath;                        
%                         doneedge(newec) = 1;
%                     else
%                         doneedge(ec) = 1;
%                     end
                    doneedge(ec) = 1;
                    
                    done(curdone) = true;
                    done(whichnode>0) = false; % nodes never marked as done
                    degreeleft(nc) = degreeleft(nc)-1;
                    degreeleft(n2) = degreeleft(n2)-1;
                    done(nodeinds(degreeleft==0)) = true; % fully done nodes
                    nedgedone = nedgedone+1;
                    break
                else
                    error('failed to map edge')
                end                
            end
        end
    end
end

NT.setupNetwork()

hold off
end