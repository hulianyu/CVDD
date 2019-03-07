
function [P, C]=k_means(S,amount,P,iter)
    for i=1:iter
        C=OptimalRepresentatives(S,P,amount);
        D=CalculateDistancesGeneral(S,C);
        P=OptimalPartition(D);   
    end
% clear D;
end

%%===========================================
function C=OptimalRepresentatives(S, P, amount)
% Update centroids using means of partitions
% Calculate new centroid value
C=zeros(amount, size(S,2));
    for i=1:amount
        members = (P == i);
        if sum(members) > 0
            if any(members)
                C(i,:)=sum(S(members,:),1)/sum(members);
            else
                warning('random centre created');
                C(i,:) = S(ceil(rand(1)*size(S,1)),:);
            end
        end
    end
end

%%=========================================
function P = OptimalPartition(D)
    [d, P] = min(D,[],1);
end

function dis=VectorDistance(v1, v2)
dis=0;
    for i=1:size(v1,2)
        dis = dis+ (v1(:,i) - v2(1,i)).^2;
    end
end

function D=CalculateDistancesGeneral(S, C)
% Calculates distances from each data point to centroids
D = zeros(size(C,1),size(S,1));
    for i=1:size(C,1)
        D(i,:)=VectorDistance(S, C(i,:));
    end  
end
