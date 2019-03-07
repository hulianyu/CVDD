function CVIs = otherCVIs(X, C, dis)

%     addpath([cd '/minmax']);
%     mx = minmax(C, sim, k);
    k = 8;
    cvnn =  CVNN(C,dis,k);
    wb = WB(X,C);
    sil = Silhouette(X, C);
    ch = CalinskiHarabasz(X, C);
    db = DaviesBouldin(X, C);
    dunn = Dunns(length(unique(C)),dis,C);
    s_dbw = S_Dbw(X, C);
    ai = I(X, C);
    CVIs =  [ cvnn, wb, sil, ch, db, dunn, s_dbw, ai] ;
end

%% 
function Sil = Silhouette(X, C)
    s = evalclusters(X,C,'Silhouette');
    Sil = s.CriterionValues;
end

function db = DaviesBouldin(X, C)
    s = evalclusters(X,C,'DaviesBouldin');
    db = s.CriterionValues;
end

function ch = CalinskiHarabasz(X, C)
    s = evalclusters(X,C,'CalinskiHarabasz');
    ch = s.CriterionValues;
end

function s_dbw = S_Dbw(X, C)
    K = length(unique(C));
    
    stdev = 0;
    for i = 1: K
        C_i = find(C == i);
        sigma_C_i = (std(X(C_i, :), 1, 1)).^2;
        s1 = sigma_C_i * sigma_C_i';
        s1 = sqrt(s1);  
        stdev = stdev + s1;
    end
    stdev = sqrt(stdev) / K;
    
    scat = 0;
    den = 0;
    for i = 1: K
        % scat
        C_i = find(C == i);
        sigma_C_i = (std(X(C_i, :), 1, 1)).^2;
        s1 = sigma_C_i * sigma_C_i';
        s1 = sqrt(s1);
        
        sigma_X = (std(X, 1, 1)).^2;
        s2 = sigma_X * sigma_X';
        s2 = sqrt(s2);   
        
        scat = scat + s1/s2;
        
        % Den
        for j = 1: K
            if i == i 
                continue;
            end
            C_j = find(C == j);
            C_ij = union(C_i, C_j);
            
            d = X(C_ij ,:) - repmat(mean(X(C_ij,:),1), length(C_ij), 1);
            d = sqrt(sum(d.* d, 2));
            d1 = sum(d < stdev);
            
            d = X(C_i,:) - repmat(mean(X(C_i,:),1), length(C_i), 1);
            d = sqrt(sum(d.* d, 2));
            d2 = sum(d < stdev);
            
            d = X(C_j,:) - repmat(mean(X(C_j,:),1), length(C_j), 1);
            d = sqrt(sum(d.* d, 2));
            d3 = sum(d < stdev);     
            
            den = den + d1/max(d2, d3);
        end
        
    end
    scat = scat / K;   
    den = den /(K*(K-1));
    s_dbw = scat + den;
end   

function ai = I(X, C)
    p = 2;
    K = length(unique(C));
    
    Diff = X - repmat(mean(X), length(C), 1);
    d1 = sum(sqrt(sum(Diff .* Diff, 2)));
    d2 = 0;
    Cen = zeros(K, size(X, 2));
    for i = 1: K
        C_i = find(C==i);
        X_i = X(C_i,:);
        Diff = X_i - repmat(mean(X_i,1), length(C_i), 1);
        d2 = d2 + sum(sqrt(sum(Diff .* Diff, 2)));     
        Cen(i, :) = mean(X_i,1);
    end
    d3 = max(pdist(Cen));
    ai = (1 / K * d1 / d2 * d3)^p;
end  

% function wb = WB(X, C)
%     m = length(unique(C));
%     C_x = mean(X);
%     W = 0;
%     B = 0;
%     for i = 1: m
%         C_i = find(C==i);
%         n = length(C_i);
%         X_i = X(C_i,:);
%         Diff_w = X_i - repmat(mean(X_i,1), length(C_i), 1);
%         W = W + sum(sqrt(sum(Diff_w .* Diff_w, 2)));     
%         Diff_b = mean(X_i,1)-C_x;
%         B = B + n*(sqrt(sum(Diff_b .* Diff_b)));    
%     end
%     wb = (m*W)/B;
% end
function wb = WB(X, C)
    m = length(unique(C));
    C_x = mean(X);
    W = 0;
    B = 0;
    for i = 1: m
        C_i = find(C==i);
        n = length(C_i);
        X_i = X(C_i,:);
        Diff_w = X_i - repmat(mean(X_i,1), length(C_i), 1);
        W = W + sum(sum(Diff_w .* Diff_w, 2));     
        Diff_b = mean(X_i,1)-C_x;
        B = B + n*(sum(Diff_b .* Diff_b));    
    end
    wb = (m*W)/B;
end


function cvnn = CVNN( C,dist,K )
    % K = 8;
    NumC = length(unique(C));
    [KNNG]=KNearestNeighborGraph(dist,K);
    cluster_weight = zeros(NumC,1);
    compactness = 0;
    for i=1:NumC
        a = find(C == i);
        b = find(C ~= i);
        num_a = length(a);
        for n=1:num_a
            if ~isempty(intersect(KNNG{a(n),1},b))
                q = length(intersect(KNNG{a(n),1},b));
                cluster_weight(i) = cluster_weight(i) + q/K;
            end     
        end
        cluster_weight(i) = cluster_weight(i)/num_a;
        compactness = compactness + (2/(num_a*(num_a-1)))*sum(sum(dist(a,a)));
    end
    separation = max(cluster_weight);
    cvnn = separation + compactness;
end

function dunn=Dunns(clusters_number,distM,ind)   
%%%Dunn's index for clustering compactness and separation measurement
% dunns(clusters_number,distM,ind)
% clusters_number = Number of clusters 
% distM = Dissimilarity matrix
% ind   = Indexes for each data point aka cluster to which each data point
% belongs
i=clusters_number;
denominator=[];

for i2=1:i
    indi=find(ind==i2);
    indj=find(ind~=i2);
    x=indi;
    y=indj;
    temp=distM(x,y);
    denominator=[denominator;temp(:)];
end

num=min(min(denominator)); 
neg_obs=zeros(size(distM,1),size(distM,2));

for ix=1:i
    indxs=find(ind==ix);
    neg_obs(indxs,indxs)=1;
end

dem=neg_obs.*distM;
dem=max(max(dem));

dunn=num/dem;
end





