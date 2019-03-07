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

function  DD = Density_involved_distance(d, K)
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
end

