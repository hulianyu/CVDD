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
filename = char('pathbased','spiral','jain','flame','aggregation','compound','r15','s1',...
    'iris','ionosphere','wine','segmentation','glass','wdbc','ecoli',...
    'movement_libras','vertebral_column','yeast','leukemia1','newthyroid');
set = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];

%% Initialization
% add Jianbo Shi's Normalized Cut
addpath([cd '/Ncut']);
%% original setting range in our paper, M = 20*2000 = 40000.
% omax = 20;
% smax = 2000;
%% a small setting range you can try
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

%% test for each dataset (we have 20 datasets in our experiments)
for f = 1:2 % select fth datasets
    %% local Initialization for each dataset
    I = set(f);
    X = load(['data_all20\', strtrim(filename(I,:)), '.txt']); %load a dataset
    N = length(X(:,1)); %Number of objects
    GT = load(['data_all20\', strtrim(filename(I,:)), '_label.txt']);%Ground truth
    NC = length(unique(GT)); %Number of clusters
    d = pdist2(X,X,'minkowski',2); %Euclidean distance of X
    Partitions = zeros(N,length(parameters)); %Partitions pi
    Max_d = max(d(:));
    C_erro = [];
    order = parameters(:,1); %parameter of W
    s = parameters(:,2); %parameter of W
    p = 0.0005; 
    %% compute the M partations by varied W  
    parfor i=1:M
        filter = d.^(order(i,1))/(p*(s(i,1))*Max_d);
        W = exp(-(filter));
        try
            [pi,~] = NcutClustering(W, NC); %partitions produced by Jianbo Shi's Normalized Cut
            error = hist(pi,unique(pi));
            Partitions(:,i) = pi;
            disp(i);
            if length(unique(pi))~= NC || min(error)<7
                C_erro = [C_erro;i];
            end
        catch
            C_erro = [C_erro;i];
        end
    end
    Pi{f,1} = Partitions;
    %% compute the index
    CAs = zeros(M,2);
    % Note: you can use 'for' or 'parfor', the computing results are the same. 
    %for i = 1:M  
    parfor i = 1:M 
        if ismember(i,C_erro)~=1
            disp(f);
            disp(i);
            CAs(i,:) = [CA(Partitions(:,i),GT),CVDD(Partitions(:,i),d)];
        end
    end
    CVI{f,1} = CAs;
    %% select OP
    [~,OP(f,1)] = max(CAs(:,2)); % the OPth partition from Pi{f,1} is selected
    OP_CA(f,1) = CVI{f,1}(OP(f,1),1);
end
clear C_erro CAs d f filename GT I M Max_d N NC omax order p Para parameters Partitions s set smax X


