% -------------------------------------------------------------------------
%Aim: The matlab code of "Generalizing centroid index to different clustering models"
%mapping by partition similarity 
%https://www.youtube.com/watch?v=WTyfjjAbUPg
% -------------------------------------------------------------------------
%Input:
%C: the partition of X
%C_Label: the ground truth of X
% -------------------------------------------------------------------------
%Output:
%results: the Centroid Index of C
% -------------------------------------------------------------------------
% Written by Lianyu Hu
% Department of Computer Science, Ningbo University 
% February 2019

function ci = CI(C, C_Label)
    
    L1 = unique(C);
    L1(:,2) = 0;
    L2 = unique(C_Label);
    L2(:,2) = 0;
    for i = 1: length(L1)
        Idx = C==L1(i);
        id = mode(C_Label(Idx));%mapping:C -> C_Label
        L2(id,2) = L2(id,2)+1;
    end
    ci_L2 = length(find(L2(:,2)==0));
    for i = 1: length(L2)
        Idx = C_Label==L2(i);
        id = mode(C(Idx));%mapping:C_Label -> C
        L1(id,2) = L1(id,2)+1;
    end
    ci_L1 = length(find(L1(:,2)==0));
    
    ci = max(ci_L2,ci_L1);
end