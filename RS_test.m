filename = char('s1','s2','s3','s4','dim32','dim64','dim128','a1','a2','unbalance');
% filename = char('pathbased','spiral','aggregation','atom','lsun','zelnik1','rings','zelnik6','triangle1','longsquare',...
%     'iris','ionosphere','wine','glass','wdbc','movement_libras','vertebral_column','yeast','leukemia1','Seeds');
%% Initialization
% add Pasi's Random Swap (RS) algorithm
addpath([cd '/RS']);


Pi_s = cell(10,1);
CVI_s = cell(10,1);
OP = zeros(10,1);
OP_s = zeros(10,1); % the OP's cluster number
for I = 1:1% select fth datasets
    %% local Initialization for each dataset
    X = load(['Datasets_all30\', strtrim(filename(I,:)), '.txt']); %load a dataset
    N = length(X(:,1)); %Number of objects
%     GT = load(['basic_benchmark\', strtrim(filename(I,:)), '_label.txt']);%Ground truth
    NC = length(unique(GT)); %Number of clusters
    d = pdist2(X,X,'minkowski',2); %Euclidean distance of X
    nc = floor(sqrt(N));
%     nc = 15;
    DD = Density_involved_distance(d, 8);%Density-involved distance of X
    Partitions = zeros(N,length(nc-1)); %Partitions pi
    C_erro = [];
    
    parfor i=1:nc-1
        disp('*********');
        disp(I);
        disp(i);
%         if i+1==NC
%             pi = GT;
%         else
%             [pi,C]=RS(X,i+1); %partitions produced by Pasi's Random Swap
%             pi = pi';
%         end
        [pi,C]=RS(X,i+1); %partitions produced by Pasi's Random Swap
        pi = pi';
%         pi = kmeans(X,i+1);
        error = hist(pi,unique(pi));
        Partitions(:,i) = pi;
        disp(i);
        if length(unique(pi))~= i+1 || min(error)<8
            C_erro = [C_erro;i];
        end
    end
    Pi_s{I,1} = Partitions;
    CAs = zeros(nc-1,10);
    parfor i=1:nc-1
        disp('##########');
        disp(I);
        disp(i);
        if ismember(i,C_erro)~=1
            CAs(i,:) = [CA(Partitions(:,i),GT),CVDD(Partitions(:,i), d, DD),otherCVIs(X, Partitions(:,i), d)];
        end
    end
    CAs(:,11) = 1:nc-1;
    CAs = CAs(CAs(:,1)~=0,:);
    CVI_s{I,1} = CAs;
    
    %% select OP
    % CVIs =  [ CA, CVDD, cvnn, wb, sil, ch, db, dunn, s_dbw, ai] ;
    [~,OP(I,1)] = max(CAs(:,1)); 
    [~,OP(I,2)] = max(CAs(:,2)); % the OPth partition from Pi{f,1} is selected
    [~,OP(I,3)] = min(CAs(:,3)); 
    [~,OP(I,4)] = min(CAs(:,4)); 
    [~,OP(I,5)] = max(CAs(:,5)); 
    [~,OP(I,6)] = max(CAs(:,6)); 
    [~,OP(I,7)] = min(CAs(:,7)); 
    [~,OP(I,8)] = max(CAs(:,8)); 
    [~,OP(I,9)] = min(CAs(:,9)); 
    [~,OP(I,10)] = max(CAs(:,10)); 
    OP_s(I,1:10) = CAs(OP(I,1:10),11)+1;
end
