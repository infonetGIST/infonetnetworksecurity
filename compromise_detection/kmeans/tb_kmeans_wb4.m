clc;clear;
%%
rng default; % For reproducibility
X = [randn(100,1)*0.75+10*ones(100,1);
    randn(100,1)*0.5-ones(100,1)
    randn(100,1)*0.25+ones(100,1)];
range = 1:length(X);
idx1 = 1:100;
idx2 = 101:200;

figure;
plot(range, X, '.');
title 'Randomly Generated Data';

opts = statset('Display','final','MaxIter',100);

tmp_sumd = zeros(9,1);
for itr = 2:10
    [idx,C,sumd] = kmeans(X,itr,'Distance','sqEuclidean',...
        'Replicates',10,'Options',opts);
    tmp_sumd(itr-1) = sum(sumd);
end

figure;
plot(2:10, tmp_sumd,'r-.', 'MarkerSize',12);
%%