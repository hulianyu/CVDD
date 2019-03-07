% [P,C]=RS(data, amount)
% Function perform RLS clustering algorithm
% Input data:
%    F - file name with data to be clasterized, file must be located
%       in the current folder. File must be in ASCII *.txt format. Data inside -
%       matrix of size [N d], where N - number of data points, d -
%       dimensions
%      - name of data matrix, already loaded to the Matlab workspace
%
%    amount - number of centroids to be selected randomly from the
%       data, default is 15
%    T - optional, the number of interations to be performed, default
%       is 1000
%    K - optional, the number of iterations of k-means algorithm,
%       default is 2 
%    outC - optional, filename to be used to save result centroids,
%       if empty, centroids.txt used
%    outP - optional, filename to be used to save result clusters,
%       if empty, clusters.txt used   
% 
%    Output:
%    P - new centroids
%    C - new clusters
%
%
% Yevgeniya Kobrina, 2010
% yevgeniya.kobrina@uef.fi

% edit: Lianyu Hu, 2019


function [P,C]=RS(S, amount)
tic;
% S = load([F '.txt']);
% outC = 'centroids.txt';
% outP = 'partitions.txt';

% T=1000;
T = 1000;
% K=2;
% K=20;
K = 20;
v=0;
% amount=15;
C = selectRepresentatives(S, amount);

[P, Do]= optimalPartition(S, C);
prev_MSE = solutionMSE(S, C, P);
if(v==1)
    h=figure();
end
for i=1:T
    % New Soulution %
    newC = C;
    newP = P;
    D = Do;
    [newReps, j, index] = randomSwap(S, newC, newP);
    % Tuning new solution %
    [newP, D] = localRepartition(S, newReps, newP, D, j, index);
       
    [newP, newC]=k_means(S,amount,newP,K);
    
    % Evaluate solution %
    curr_MSE = solutionMSE(S, newC, newP);
    if (curr_MSE < prev_MSE)
        prev_MSE = curr_MSE;
        C = newC;
        P = newP;
        Do=D;
        accepted=1;
    end
    % Visualize solution %
    if(v~=0 && accepted==1)
        Visualize(h,S,P,C, i);
    end
end
    tElapsed = toc;
    ShowResults(tElapsed, prev_MSE);
%     WriteOutput(outC, C, outP, P);
end


%%==================================================
% internal functions
function newC = selectRepresentatives(S, amount)
    % Creates random codebook
    newC=S(ceil(rand(amount,1)*size(S,1)),:);
end

function [P, D] = optimalPartition(S, C)
    % Does a partition of data S
    dis=CalculateDistancesGeneral(S, C);
    clear S C;
    [D, P] = min(dis,[],1);
    clear dis;
end

function [C, j, i] = randomSwap(S, C, P)
    % select new centroid
    s = size(S);
    i=ceil(rand(1) * s(1));
    % select centroid to be deleted 
    s = size(C);
    j=ceil(rand(1) * s(1));
    ms = P==j;
    if (ms(i)==1 || sum(ms)==0)
       % recursive call until the appropriate cluster found 
%        [newC, j, i] = randomSwap(S, C, P); 
       [~, j, i] = randomSwap(S, C, P); 
    else    
        C(j, :) = S(i,:);
    end    
end

function [P, D] = localRepartition(S, C, P, D, j, index)
   [P, D]=RepartitionOldCluster(S, P, C, D, j);
   [P, D]=RepartitionNewCluster(S, P, C, D, index, j);
end

function [P, D]=RepartitionOldCluster(S, P, C, D, j)
    % retrive indexes which belong to cluster j
    relIndexes = P(:)==j;
    % get the points of the cluster j
    releasedCluster = S(relIndexes,:);
    % calculate distance from each point of cluster j to each centroid   
     dis=CalculateDistancesGeneral(releasedCluster, C);  
    % retrive indexes of C with minimal distance to points of
    % cluster j
    [Drel, NC] = min(dis,[],1);
    clear dis releasedCluster;
    % update data in main cluster array
    P(relIndexes) = NC;
    D(relIndexes) = Drel;
end

function [P, D]=RepartitionNewCluster(S, P, C, D, index, j)
    %reallocate cluster
    % get the cluster id which should be divided to 2 new
    clusterId = P(index);
    % retrive indexes which belong to cluster clusterId
    relIndexes = P(:)==clusterId;
    centroindsCount = size(C,1); 
    % get the points of the cluster clusterId
    releasedCluster = S(relIndexes,:);
    % calculate distance from each point of cluster to clusterId centroid
    % infinity required to emphasize only required C
    dis = inf(centroindsCount,size(releasedCluster,1));
    % calculate distance from each point of cluster to clusterId centroid
    dis(clusterId,:)=VectorDistance(releasedCluster, C(clusterId,:));
    % calculate distance from each point of cluster to j centroid
    dis(j,:)=VectorDistance(releasedCluster, C(j,:));
    % retrive indexes of C with minimal distance to points
    [Drel, NC] = min(dis,[],1);
    P(relIndexes) = NC;
    D(relIndexes) = Drel;
end

% function dis=VectorDistance(v1, v2)
% dis=0;
%     for i=1:size(v1,2)
%         dis = dis+ (v1(:,2) - v2(1,2)).^2;
%     end
% end
function dis=VectorDistance(v1, v2)
dis=0;
    for i=1:size(v1,2)
        dis = dis+ (v1(:,i) - v2(1,i)).^2;
    end
end

function dis=CalculateDistancesGeneral(S, C)
    % Calculates distances from each data point to each centroid
    dis = zeros(size(C,1),size(S,1));
        for i=1:size(C,1)
            dis(i,:)=VectorDistance(S,C(i,:));
        end    
    clear S;
end
function value = solutionMSE(S, C, P)
    % sum of squared distances of the data object to their cluster
    % representatives
    dis=zeros(size(S,1),1);
    for i=1:size(S,2)
        dis = dis+(S(:,i) - C(P(:),i)).^2;
    end
    value = sum(dis);
    clear dis;
end
function Visualize(h, S,P,C,i)
    figure(h);
    drawData(S,P); 
    title(sprintf('Solution selected on %d iteration',i));
    hold on;
    drawCentroids(C); 
    hold off;
end

function drawData(S, P)
    marker=2;
    scatter(S(:,1),S(:,2),marker,P); 
end

function drawCentroids(C)
   scatter(C(:,1),C(:,2),'ko','filled'); 
end
function ShowResults(t, MSE)
    disp('Time:');
    disp(t);
    disp('MSE:');
    disp(MSE);

end
% function WriteOutput (outC, C, outP, P)
%     save(outC, 'C', '-ASCII', '-tabs');
%     save(outP, 'P', '-ASCII', '-tabs');
% end
