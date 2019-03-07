% -------------------------------------------------------------------------
%Aim:
%The matlab code of "An internal validity index based on density-involved distance"
%Algorithm 2: CVDD-OP
%
%Output:OP
%Each row in CVI,Pi corresponds a dataset
%Each column in Pi{f,1} corresponds a partition, each row in CVI{f,1} corresponds a partition
% -------------------------------------------------------------------------
% Written by Lianyu Hu
% Department of Computer Science, Ningbo University 
% August 2018
%
% function "parfor" is in MATLAB which can speed up the loops [Parallel for loop]

%% load dataset
filename = char('pathbased','spiral','aggregation','atom','lsun','zelnik1','rings','zelnik6','triangle1','longsquare',...
    'iris','ionosphere','wine','glass','wdbc','movement_libras','vertebral_column','yeast','leukemia1','Seeds');

%% Initialization
% add Jianbo Shi's Normalized Cut
addpath([cd '/Ncut']);
% %% original setting range in our paper, M = 20*2000 = 40000.
% omax = 20;
% smax = 2000;
% a small setting range you can try
omax = 2;
smax = 200;
%% global Initialization
M = omax*smax;
Para = zeros(M,2);
Pi = cell(20,1);
CVI = cell(20,1);
OP = zeros(20,1);
OP_CA = zeros(20,1); % the OP's CA
Para(:,1) = reshape(repmat(1:omax, smax, 1), omax*smax, 1);
Para(:,2) = repmat((1:smax)',omax,1);
parameters = Para;

%% test for each dataset (we have 20 non-spherical clusters and classification datasets in our experiments)
for f = 1:2% select fth datasets
    %% local Initialization for each dataset
%     I = set(f);
    I = f;
    X = load(['Datasets_all30\', strtrim(filename(I,:)), '.txt']); %load a dataset
    N = length(X(:,1)); %Number of objects
    GT = load(['Datasets_all30\', strtrim(filename(I,:)), '_label.txt']);%Ground truth
    NC = length(unique(GT)); %Number of clusters
    d = pdist2(X,X,'minkowski',2); %Euclidean distance of X
    DD = Density_involved_distance(d, 8);%Density-involved distance of X
    Partitions = zeros(N,length(parameters)); %Partitions pi
    Max_d = max(d(:));
    C_erro = [];
    order = parameters(:,1); %parameter of W
    s = parameters(:,2); %parameter of W
    p = 0.0005; 
%     p = 0.001;
    %% compute the M partations by varied W  
    for i=1:M
%     parfor i=1:M
        filter = d.^(order(i,1))/(p*(s(i,1))*Max_d);
        W = exp(-(filter));
        try
            [pi,~] = NcutClustering(W, NC); %partitions produced by Jianbo Shi's Normalized Cut
            error = hist(pi,unique(pi));
            Partitions(:,i) = pi;
            disp(i);
            if length(unique(pi))~= NC || min(error)<8
                C_erro = [C_erro;i];
            end
        catch
            C_erro = [C_erro;i];
        end
    end
    Pi{I,1} = Partitions;
    %% compute the index
    CAs = zeros(M,10);
    % Note: you can use 'for' or 'parfor', the computing results are the same. 
    for i = 1:M  
%     parfor i = 1:M 
        if ismember(i,C_erro)~=1
            disp(I);
            disp(i);
            CAs(i,:) = [CA(Partitions(:,i),GT),CVDD(Partitions(:,i), d, DD),otherCVIs(X, Partitions(:,i), d)];
        end
    end
    
    CAs(:,11) = 1:M;
    CAs = CAs(CAs(:,1)~=0,:);
    CVI{I,1} = CAs;
    %% select OP
    %% CVIs =  [ CA, CVDD, cvnn, wb, sil, ch, db, dunn, s_dbw, ai] ;
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

    OP_CA(I,1:10) = CAs(OP(I,1:10),1);
%     save
end
clear C_erro d f filename GT I M Max_d N NC omax order p Para parameters Partitions s set smax
