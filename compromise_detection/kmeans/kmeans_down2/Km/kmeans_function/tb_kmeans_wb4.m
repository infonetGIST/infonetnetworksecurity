clc;clear;
%%
rng default; % For reproducibility
X = [randn(100,1)*0.75+10*ones(100,1);
    randn(100,1)*0.5-ones(100,1)
    randn(100,1)*0.25+ones(100,1)];
range = 1:length(X);
idx1 = 1:100;
idx2 = 101:200;

%% K-means options
feature_vector     = X;                                 % Input
number_of_clusters = 3;                                 % Number of Clusters
Kmeans_iteration   = 40;                                % K-means Iteration
%% Test K-means
[cluster_centers, data]  = km_fun(feature_vector, number_of_clusters, Kmeans_iteration); % K-means clusterig
est_idx1 = find(data(:,number_of_clusters+1)==1);
est_idx2 = find(data(:,number_of_clusters+1)==2);
est_idx3 = find(data(:,number_of_clusters+1)==3);
%% Plot   
figure;
plot(est_idx1,X(est_idx1,1),'r.','MarkerSize',12)
hold on
plot(est_idx2,X(est_idx2,1),'b.','MarkerSize',12)
plot(est_idx3,X(est_idx3,1),'g.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
title 'Cluster Assignments and Centroids'
hold off