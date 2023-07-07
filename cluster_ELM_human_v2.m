clc; 
clear; 
close all;

%% Loading ELM-human power-law data
load('table1_p11.mat')

%% Creating domain and co-domain for clustering purpsoes
X = [data.k_c_real(:); data.k_c_elm(:)]; % x axis is critical curvature
Y = [data.v_c_real(:); data.v_c_elm(:)]; % y axis is critical speed

N = length(data.k_c_real(:)); % number of data points, i.e., number of subjects in the experiment

%% Plotting the ELM-human power-law data
figure(1);
sz = 300;
scatter(X(1:N),Y(1:N),sz,'r.');hold on,
scatter(X(N+1:end),Y(N+1:end),sz,'b.');
xlabel('Curvature (1/m)')
ylabel('Speed (m/s)')

%% k-means clustering to detect k=2 clusters in the ELM-human power-law data

clustering_case = input('Press 1 for kmedoids, 2 for kmeans\n');
    
rng default;
K = 2;
data = [X, Y];

switch clustering_case
    case 1
        [clusterIndices, centroids] = kmedoids(data, K, 'Options', statset('Display','iter'));
    case 2
        [clusterIndices, centroids] = kmeans(data, K, 'Display', 'iter');
end

er_num = 0;
er_num_op = 0;
for idx = 1:N*2
    if (clusterIndices(idx) == 1) && idx <= N
        er_num = er_num + 1;
    elseif (clusterIndices(idx) == 2) && idx > N
        er_num = er_num + 1;
    end
    if (clusterIndices(idx) == 2) && idx <= N
        er_num_op = er_num_op + 1;
    elseif (clusterIndices(idx) == 1) && idx > N
        er_num_op = er_num_op + 1;
    end
end

rng default;
figure(2);
x1 = min(X):0.001:max(X);
x2 = min(Y):0.001:max(Y);
[x1G,x2G] = meshgrid(x1,x2);
XGrid = [x1G(:),x2G(:)]; % Defines a fine grid on the plot

switch clustering_case
    case 1
        idx2Region = kmedoids(XGrid,2,'Options',statset('Display','iter','MaxIter',1));
        fig2_title = 'kmedoids';
    case 2
        idx2Region = kmeans(XGrid,2,'MaxIter',1,'Start',centroids);
        fig2_title = 'kmeans';
end

gscatter(XGrid(:,1),XGrid(:,2),idx2Region,...
    [252/255, 225/255, 244/255;0.3010 0.7450 0.9330;0.3010 0.7450 0.9330],'..');
hold on;
scatter(X(1:N),Y(1:N),sz,'r.');hold on,
scatter(X(N+1:end),Y(N+1:end),sz,'b.');
xlabel('Curvature (1/m)')
ylabel('Speed (m/s)')
hold on
plot(centroids(:,1),centroids(:,2),'kx',...
     'MarkerSize',7,'LineWidth',2);
legend('Cluster 1','Cluster 2','Real true','ELM true','Centroids','Location','SouthEast');
title(fig2_title)
hold off;

%% Guassian mixture with 2 components

rng default;
options = statset('MaxIter',1000);

Sigma = {'diagonal','full'}; % Options for covariance matrix type
nSigma = numel(Sigma);

SharedCovariance = {true,false}; % Indicator for identical or nonidentical covariance matrices
SCtext = {'true','false'};
nSC = numel(SharedCovariance);

d = 500; % Grid length
x1 = linspace(min(X)-2, max(X)+2, d);
x2 = linspace(min(Y)-2, max(Y)+2, d);
[x1grid,x2grid] = meshgrid(x1,x2);
X0 = [x1grid(:) x2grid(:)];

figure(3);
threshold = sqrt(chi2inv(0.99,2));
count = 1;
for i = 1:nSigma
    for j = 1:nSC
        gmfit = fitgmdist(data,K,'CovarianceType',Sigma{i}, ...
            'SharedCovariance',SharedCovariance{j},'Options',options); % Fitted GMM
        clusterX = cluster(gmfit,data); % Cluster index 
        mahalDist = mahal(gmfit,X0); % Distance from each grid point to each GMM component
        % Draw ellipsoids over each GMM component and show clustering result.
        subplot(2,2,count);
        h1 = gscatter(X,Y,[ones(10,1);ones(10,1)*2]);
        hold on
            for m = 1:K
                idx = mahalDist(:,m)<=threshold;
                Color = h1(m).Color*0.75 - 0.5*(h1(m).Color - 1);
                h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
                uistack(h2,'bottom');
            end    
        %plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        title(sprintf('Sigma is %s\nSharedCovariance = %s',Sigma{i},SCtext{j}),'FontSize',8)
        legend(h1,{'Real','ELM'})
        % hold on
        % scatter(X(1:N),Y(1:N),100,'r*');hold on,
        % scatter(X(N+1:end),Y(N+1:end),100,'b*');
        hold off
        count = count + 1;
    end
end
