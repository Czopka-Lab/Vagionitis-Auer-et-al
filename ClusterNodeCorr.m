clear 
% 
% the distribution of intercluster distance is 18.2±8 µm (measured values from Figure 1G)
% the distribution of internode distance is 30.9±13.2 µm (measured values from Figure 1H)


numSim = 1000;                  % the simulation is looped for 1000 times
NodeClusterWidth = 2;           % the node and cluster width is assumed to be 2um 

% clusters and nodes are counted as colocalizing when the difference of the
% cluster and node middle is less than 2um (in other words there is a part
% of the node and cluster that is overlapping which is the same criterion that was used
% for the data analysis)





%% prospective analysis: does the position of a cluster correspond with the location of a node

numCluster = 32;  % as we could analyse 32 clusters that remained and axons got myelinated 
numNodes = 32;    

for ll = 1:numSim
    
    % simulate axons with cluster and myelinated axons with nodes
    
    width = NodeClusterWidth/2;  % initalizes the width, width starts at 1 (half the cluster width) to get the middle of the cluster 
    clusterPosition = 0;   % initalizes the cluster position;   
    for ii = 2:numCluster+1      
       intClusterdist(ii) = normrnd(18.2,8.0);  % randomly generates intercluster distances from the given distribution
       clusterPosition(ii) = clusterPosition(ii-1) + intClusterdist(ii)+width; %width is added to add the cluster width as the distribution only gives the intercluster distance
       width=NodeClusterWidth; % node width
    end

    widthnode = NodeClusterWidth; % initalizes the width, width starts at 1 (half the node width) to get the middle of the node 
    nodePosition = 0;      % initalizes the node position;   
    for ii = 2:numNodes+1       
      intNodedist(ii) = normrnd(30.9,13.2);
      nodePosition(ii) = nodePosition(ii-1) + intNodedist(ii)+widthnode; %add nodewidth as the distribution only gives the internode distance
      widthnode=NodeClusterWidth; 
    end
    
    % checks if the simulated axon with nodes is at least as long as
    % the simulated axon with clusters - if not it stops the simulation
    if max(nodePosition)-max(clusterPosition)<0
        disp('increase the number of nodes to be generated for prospective analysis');
        break
    end
    
    % analyses if the the position of a cluster corresponds with the location of a node
    
    CoLo(ll) = 0; % initlizes the number of colocalized cluster and nodes
    
    for ii = 2:length(clusterPosition)              % starts at 2 as the first position is always set to 0 and will not be analysed
        for jj = 2:length(nodePosition)
            if abs(clusterPosition(ii) - nodePosition(jj)) < NodeClusterWidth     %checks if the cluster position correlates with the cluster/node position, if the difference is smaller than the node width, and a width of 2um is assumed the two overlap *(same criteria was used for the analysisof the real data)
                CoLo(ll)= CoLo(ll) + 1;             % if the above criterion is met, it counts it as a colocalization of cluster and node
            end
        end
    end
        


end
disp('Prospective Analysis')
m = round(mean(CoLo),2);
disp(sprintf('From %d analysed clusters, the mean number of clusters that colocalize with nodes is %.f',numCluster, m));
ratio = (mean(CoLo)/numCluster)*100 ;% gives percentage of clusters and nodes that colocalized
SD=(std(CoLo)/numCluster)*100;
disp(sprintf('The percentage of clusters that colocalize with nodes is %.f +- %.f ', ratio, SD));
maxCoLo= (max(CoLo)/numCluster)*100;
disp(sprintf('The maximum percentage of colocalization observed by chance is %.f',maxCoLo))

%% plots a histogram for the relative number of colocalized clusters-nodes
figure
subplot(2,2,1)
PercentCoLo =((CoLo/numCluster)*100);
h1 = histogram(PercentCoLo,'FaceColor',[0.7 0.7 0.7]);
h1.BinWidth = 3;
hold on
title('distribution of % of cluster that became a node at the same position')
ylabel('probability/number of observations')
xlabel('% of cluster that became a node at the same position')
legend('simulated data')
subplot(2,2,2)
[f,xi] = ksdensity(PercentCoLo);
xline(mean(PercentCoLo),'k','LineWidth',2)
hold on
plot(xi,f,'k','LineWidth',3)
hold on
xline(37.5,'r','LineWidth',5) % plots the measured value from the real data
xline(62.5,'Color',[1 0.4 0.4],'LineWidth',5) % plots the measured value from the real data - including heminodes
title('distribution of % of cluster that became a node at the same position')
ylabel('probability')
xlabel('% of cluster that became at the same position')
legend('probability density','mean value - simulated','measured value - nodes','measured value - nodes and heminodes')

%% retrospective anlysis: does the position of a node correspond with the location of a cluster

numNodes = 33;    % in the retrospective anlysis 33 nodes were analysed
numCluster = 70;  % we have to generate more clusters than nodes as their distance is
% shorter than node distance, to model a similar long axon stretch for both axons

for ll = 1:numSim
    
    % simulate axons with cluster and myelinated axons with nodes
    width = NodeClusterWidth; % width starts at half the cluster width to get the middle of the cluster 
    clusterPosition = 0;
    intClusterdist = 0;
    for ii = 2:numCluster+1      
       intClusterdist(ii) = normrnd(18.2,8.0);
       clusterPosition(ii) = clusterPosition(ii-1) + intClusterdist(ii)+width; %x is added to add the cluster width as the distribution only gives the intercluster distance
       width=NodeClusterWidth; 
    end

    widthnode = NodeClusterWidth; % width starts at half the node width to get the middle of the node 
    nodePosition = 0;
    intNodesdist  = 0;
    for ii = 2:numNodes+1      
      intNodedist(ii) = normrnd(30.9,13.2);
      nodePosition(ii) = nodePosition(ii-1) + intNodedist(ii)+widthnode; %x is added to add the cluster width as the distribution only gives the intercluster distance
      widthnode=NodeClusterWidth; 
    end
    % checks if the simulated axons with clusters is at least as long as
    % the simulated axon with nodes - if not it stops the simulation
    if max(clusterPosition)-max(nodePosition)<0
        disp('increase the number of clusters to be generated for retrospective analysis');
        break
    end
        
    
    % analyses if the the position of a cluster corresponds with the location of a node
    CoLoRetro(ll) = 0; % initlizes the number of colocalized cluster and nodes
    for ii = 2:length(nodePosition)
        for jj = 2:length(clusterPosition)
            if abs(nodePosition(ii)-clusterPosition(jj)) < NodeClusterWidth     %checks if the node position correlates with the cluster position, if the difference is smaller than the node width, and a width of 2um is assumed the two overlap *(same criteria was used for the analysisof the real data)
                CoLoRetro(ll)= CoLoRetro(ll) + 1;
            end
        end
    end


end

% print the results in text form

disp('Retrospective Analysis')
m = round(mean(CoLoRetro),2);
disp(sprintf('From %d analysed nodes, the mean number of nodes that colocalize with clusters is %.f',numNodes, m));
ratio = (mean(CoLoRetro)/numNodes)*100 ;% gives percentage of clusters and nodes that colocalized
SD=(std(CoLoRetro)/numNodes)*100;
disp(sprintf('The percentage of nodes that colocalize with clusters is %.f +- %.f ', ratio, SD));
maxCoLoRetro= (max(CoLoRetro)/numNodes)*100;
disp(sprintf('The maximum percentage of colocalization observed by chance is %.f',maxCoLoRetro))



%% generate histogram for retrospective analysis
subplot(2,2,3)
% calculate percentage of colocalization
PercentRetro = ((CoLoRetro/numNodes)*100);
h2 = histogram(PercentRetro,'FaceColor',[0.4 0.4 0.4]);
h2.BinWidth = 3;
title('distribution of % of nodes that had a cluster at the same position')
ylabel('probability/number of observations')
xlabel('% of nodes that had a cluster at the same position')
legend('simulated data')
subplot(2,2,4)
[f,xi] = ksdensity(PercentRetro);
xline(((m/numNodes)*100),'k','LineWidth',2)
hold on
plot(xi,f,'k','LineWidth',3)
hold on
xline(55,'r','LineWidth',5) % plots the measured value from the real data - only nodes
xline(66.7,'Color',[1 0.4 0.4],'LineWidth',5) % plots the measured value from the real data - nodes and heminodes
title('distribution of % of nodes that had a cluster at the same position')
ylabel('probability')
xlabel('% of nodes that had a cluster at the same position')
legend('probability density','mean value - simulated','measured value - nodes','measured value - nodes and heminodes')


max(CoLo)
max(CoLoRetro)
