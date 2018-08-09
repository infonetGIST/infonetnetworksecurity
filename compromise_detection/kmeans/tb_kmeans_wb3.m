clc;clear;
%%
rng default; % For reproducibility
X = [randn(100,1)*0.75+10*ones(100,1);
    randn(100,1)*0.5-ones(100,1)
    randn(100,1)*0.25+ones(100,1)];
range = 1:length(X);
idx1 = 1:100;
idx2 = 101:200;

% figure;
% plot(idx1, X(idx1),'r.','MarkerSize',12);
% hold on;
% plot(idx2, X(idx2),'b.','MarkerSize',12);

figure;
plot(range, X, '.');
title 'Randomly Generated Data';

opts = statset('Display','final','MaxIter',100);
[idx,C,sumd,D] = kmeans(X,3,'Distance','sqEuclidean',...
    'Replicates',10,'Options',opts);
est_idx1 = find(idx==1);
est_idx2 = find(idx==2);
est_idx3 = find(idx==3);
% est_idx4 = find(idx==4);

figure;
plot(est_idx1,X(est_idx1,1),'r.','MarkerSize',12)
hold on
plot(est_idx2,X(est_idx2,1),'b.','MarkerSize',12)
plot(est_idx3,X(est_idx3,1),'g.','MarkerSize',12)
% plot(est_idx4,X(est_idx4,1),'c.','MarkerSize',12)
% plot(C(:,1),C(:,2),'kx',...
%      'MarkerSize',15,'LineWidth',3)
% legend('Cluster 1','Cluster 2','Centroids',...
%        'Location','NW')
title 'Cluster Assignments and Centroids'
hold off


WSS= sum(sumd);

%%