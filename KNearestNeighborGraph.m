function [KNNG]=KNearestNeighborGraph(dist,K)

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