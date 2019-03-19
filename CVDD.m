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

% function CVDD = CVDD_DDDE(pi,X, d) 
function CVDD = CVDD(pi,d,DD) 
%     %% compute the density-involved distance, K=7
%     try
%         DD = DDDE_involved_distance(X, d, K);
%         DD_erro = 0;
%     catch
%         DD_erro = 1;
%     end
%     if DD_erro~=1
       %% initialization
        NC = length(unique(pi));
        sc_list = zeros(NC,1); %separation
        com_list = zeros(NC,1);%compactness
       %% 
        for i = 1: NC
            a = find(pi == i);
            b = pi ~= i;
            n = length(a);
            if isempty(a)~=1
              %% compute the separation sep[i]
                sc_list(i,1) = min(min(DD(a,b)));
              %% compute the compactness com[i]
                try
                    Ci = fast_PathbasedDist(d(a,a));
                    com_list(i,1) = (std2(Ci)/n)*mean(Ci(:));
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
%     else
%         CVDD = 0;
%     end
end

