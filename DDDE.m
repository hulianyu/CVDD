function Density = DDDE(X,K)

% Dimensionally distributed density estimation (DDDE) 
% - Computes the density of a set of points X
% CALL:
% Density = DDDE(X,K)
% 
% INPUT:
% X: a dataset
% K: number of neighbors
%              
% 
% Last Rev : Feb 07 2019 by Lianyu Hu
%            hly4ml@gmail.com

if mod(K,2)==1
    K = K+1;
end
Dim = size(X,2);
N = size(X,1);
Dens = zeros(N,1);
Density = zeros(N,1);
for dim = 1:Dim
    Z = X(:,dim);
    [Y,object] = sort(Z);
    S_n = sum(Y(1:K/2));
    S_p = sum(Y(K/2+2:K+1));

    Dens(1,1) = (sum(Y(2:K+1))-K*Y(1))/K;
    Density(object(1),1) = Density(object(1),1) + Dens(1,1);
    for i = 2:K/2
        Dens(i,1) = ((i-1)*Y(i)-sum(Y(1:i-1)) + sum(Y(i+1:K+1))-(K-i+1)*Y(i))/K;
        Density(object(i),1) = Density(object(i),1) + Dens(i,1);
    end
    for i = N-K/2+1:N
        Dens(i,1) = (sum(Y(i+1:N))-(N-i)*Y(i) + (K-N+i)*Y(i)-sum(Y(N-K:i-1)))/K;
        Density(object(i),1) = Density(object(i),1) + Dens(i,1);
    end
    for i = K/2+1:N-K/2
        Dens(i,1) = (S_p - S_n)/K;
        Density(object(i),1) = Density(object(i),1)+ Dens(i,1);
        if i~= N-K/2
        S_n = S_n + Y(i) - Y(i-K/2);
        S_p = S_p - Y(i+1) + Y(i+K/2+1);
        end
    end
%     Density = Density/max(Density);
end