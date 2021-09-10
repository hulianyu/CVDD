% Cluster Validity index based on Density-involved Distance
%   MAIN REFERENCE:
%       - L. Hu and C. Zhong, “An internal validity index based on density-involved distance,” IEEE Access, vol. 7, pp. 40038–40051, 2019, doi: 10.1109/ACCESS.2019.2906949.

function cvddindex = CVDDIndex(X, piX)
d = pdist2(X,X,'minkowski',2); %Euclidean distance of X
try
    DD = Density_involved_distance(d, 7); %Density-involved distance of X
    cvddindex = CVDD(piX, d, DD);
catch
    % an error will be catch when the Rel function fails
    % because of 0 distances between some points (total overlap)
    % in this case, the minimum value of the index will be returned
    cvddindex = 0;
end

function CVDD = CVDD(piX,d,DD)
% -------------------------------------------------------------------------
%Aim: The matlab code of "An internal validity index based on density-involved distance"
%Algorithm 1: CVDD
% -------------------------------------------------------------------------
%Input:
%pi: the partition of X
%d: the Euclidean distance between objects in X
% -------------------------------------------------------------------------
%Output:
%results: the CVDD index of pi
% -------------------------------------------------------------------------
% Written by Lianyu Hu
% Department of Computer Science, Ningbo University
% February 2019

%% initialization
NC = length(unique(piX));
sc_list = zeros(NC,1); %separation
com_list = zeros(NC,1);%compactness
%%
for i = 1: NC
    a = find(piX == i);
    b = piX ~= i;
    n = length(a);
    if isempty(a)~=1
        %% compute the separation sep[i]
        sc_list(i,1) = min(min(DD(a,b)));
        %% compute the compactness com[i]
        try
            Ci = fast_PathbasedDist(d(a,a));
            % std2 requires Image Processing Toolbox
            com_list(i,1) = (std2(Ci)/n)*mean(Ci(:));
            %com_list(i,1) = (std(Ci(:))/n)*mean(Ci(:));
        catch
            com_list(i,1) = max(com_list);
        end
    else
        sc_list(i,1)=0;
        com_list(i,1) = max(com_list);
    end
end
%% compute the validity index CVDD
sep = sum(sc_list);
com = sum(com_list);
CVDD = sep/com;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extra functionalities
% Obtained from: https://github.com/hulianyu/CVDD (Original index implementation)
function  DD = Density_involved_distance(d, K)
% -------------------------------------------------------------------------
%Aim: The matlab code of "An internal validity index based on density-involved distance"
%compute the Density-involved distance DD
% -------------------------------------------------------------------------
%Input:
%d: the Euclidean distance between objects in X
%K: the number of neighborhoods
% -------------------------------------------------------------------------
%Output:
%results: the DD of d
% -------------------------------------------------------------------------
% Written by Lianyu Hu
% Department of Computer Science, Ningbo University
% August 2018

%% compute density-involved distance
N = length(d);
% K = 7;
[KNNG]=KNearestNeighborGraph(d,K);
Den = zeros(N,1);
for j = 1:N
    Den(j,1) = sum(d(j,KNNG{j,1}))/K;
end
fDen = Den/max(Den); %absolute-density distance factor
Rel = repmat(Den,1,N)./repmat(Den',N,1);
tmp1 = Rel + Rel';
fRel = 1-exp(-abs(tmp1-2)); %relative-density distance factor
nD = repmat(Den,1,N) + repmat(Den',N,1);
relD = nD.*fRel;
drD = d + relD; %directly density-reachable distance
conD =  fast_PathbasedDist(drD); %connectivity distance
tmp2 = sqrt(repmat(fDen,1,N).*repmat(fDen',N,1));
DD = conD.*tmp2; %density-involved distance

function [KNNG] = KNearestNeighborGraph(dist,K)
% KNearestNeighborGraph - Computes the K Nearest Neighbor graph
%
% Computes the K Nearest Neighbor Graph  (KNNG) of a set of points
% KNNG is a DIRECTED graph
%
% CALL:
% [KNNG]=KNearestNeighborGraph(dist,K)
%
% INPUT:
% dist: NxN Euclidean distance matrix between each pair of points
%               dist(i,k)=norm(data(i,:)-data(k,:), (dist(k,k)=0)
%
% OUPUT:
% KNNG: N cells KNNG{i}=[a b...f] set of index in data rows of the KNNG
% neighbors of i
%
% Author   : Michael Aupetit
%            Qatar Computing Research Institute (QCRI)
%            Hamad Bin Khalifa University (HBKU)
%            Doha, Qatar
%            maupetit@qf.org.qa/michael.aupetit@gmail.com
% Last Rev : Feb 08 2016 by Michael Aupetit

[N,N]=size(dist);
dist(eye(N)==1)=inf;
[val, indSort]=sort(dist);
KNNG=cell(N,1);

for i=1:N
    KNNG{i}=indSort(1:K,i)';
end

function PathbasedW = fast_PathbasedDist(W)
%% To generate pathbased distance matrix, i.e. minmax distance matrix
% Input: W -- Dissimilarity matrix
% Output: PathBased -- PathBased dissimilarity matrix
% Caiming Zhong, 2014/05/08
% Lianyu Hu, 2018/09/08
% function "minspantree" is in MATLAB

%Pairs = MST(W);
MST0 = minspantree(graph(W));
Pairs = MST0.Edges.EndNodes;
N = size(W,1);
PathbasedW = zeros(N);
Degrees = zeros(N, 1);

Connections = zeros(N);
for i = 1: N - 1
    Degrees(Pairs(i,1)) = Degrees(Pairs(i,1)) + 1;
    Degrees(Pairs(i,2)) = Degrees(Pairs(i,2)) + 1;
    Connections(Pairs(i,1), Degrees(Pairs(i,1))) = Pairs(i,2);
    Connections(Pairs(i,2), Degrees(Pairs(i,2))) = Pairs(i,1);
end

for i = 1: N
    if Degrees(i) == 1  % first node with one edge
        break;
    end
end

CloseT = zeros(N,1);
OpenT = zeros(N,1);
cursor_C = 0;
cursor_T = 0;

cursor_C = cursor_C + 1;
CloseT(cursor_C) = i;
Visited = zeros(N, 1);
Visited(i) = 1;

while cursor_C < N
    % nodes connected to the current node in close table
    Nodes1 = Connections( CloseT(cursor_C), 1: Degrees(CloseT(cursor_C)));
    Nodes = Nodes1(Visited(Nodes1)==0);

    for i = 1: length(Nodes)
        % Nodes to Close table
        for j = 1: cursor_C
            PathbasedW(CloseT(j), Nodes(i)) = max(PathbasedW(CloseT(j), CloseT(cursor_C)), W(CloseT(cursor_C),  Nodes(i)));
            PathbasedW(Nodes(i), CloseT(j) )= PathbasedW(CloseT(j), Nodes(i)) ;
        end

        % Nodes to Open table
        for j = 1: cursor_T
            PathbasedW(OpenT(j), Nodes(i)) = max(PathbasedW(OpenT(j), CloseT(cursor_C)), W(CloseT(cursor_C),  Nodes(i)));
            PathbasedW(Nodes(i), OpenT(j) )= PathbasedW(OpenT(j), Nodes(i)) ;
        end
    end

    % Nodes each other
    for i = 1: length(Nodes) - 1
        for j = i  + 1: length(Nodes)
            PathbasedW(Nodes(i), Nodes(j)) = max(W(Nodes(i), CloseT(cursor_C)), W(Nodes(j), CloseT(cursor_C)));
            PathbasedW(Nodes(j), Nodes(i)) = PathbasedW(Nodes(i), Nodes(j));
        end
    end
    for i = 1: length(Nodes)
        % Nodes are added into OpenT
        cursor_T = cursor_T + 1;
        OpenT(cursor_T) = Nodes(i);
    end
    cursor_C = cursor_C + 1;
    CloseT(cursor_C) = OpenT(cursor_T);
    Visited(OpenT(cursor_T)) = 1;
    cursor_T = cursor_T - 1;
end
