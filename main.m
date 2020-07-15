% main file to analyze SCC dataset
% load data
clear;
clc;
echo on;

%load the cell-cell similarity matrix from consensus clustering
M = csvread('SCC_cell-cell.csv');
%load the scaled dataset with selected 3000 genes
data = dlmread('/Volumes/GoogleDrive/Shared drives/Intermediate state/code_new/SCC/SCC3.txt','',1,1)';
%load the labels of cells in given in the original paper
true_label = csvread('SCC3_celltype.csv',0,1);
%remove the labels of unused cells
cells_removed = csvread('SCC3_cell_removed.txt',1,0);
true_label(cells_removed)=[];
%load selected gene names
fileID = fopen('SCC3_gene_annotations.txt');
y = textscan(fileID,'%s','HeaderLines',1,'Delimiter','\n');
fclose(fileID);
gene_name = strings(size(y{:}));
gene_name(:) = y{:};
index = find(true_label~=1 & true_label~=2);

%% eigenvalue
% The number of clusters is estimated by analyzing the largest gap of consecutive eigenvalues of symmetric normalized graph Laplacian
[eigenvalues] = plot_eigen_gap(M);
%% using symnmf on M
addpath('NNDSVD');
addpath('symnmf2');
No_cluster = 4; %input the number of cluster
flag = 1;
[M0,~] = nndsvd(M,No_cluster,flag);
INIT.Hinit = M0;
INIT.tol = 10^(-6);
%used symmetric NMF to get the cell-cluster matrix H
[H,~,~]=symnmf_newton(M,No_cluster,INIT);
%C is the hard clustering result of cells
[~,C] = max(H,[],2);
unique(C)
%% get CPI of each cell
cpi = zeros(size(data,1),1);
H_normailed =  H./sum(H,2);
for i = 1:length(cpi)
    for j = 1:No_cluster
        if H_normailed(i,j)>0
            cpi(i) = cpi(i) - H_normailed(i,j)*log(H_normailed(i,j));
        end
    end
    cpi(i) = cpi(i)/log(No_cluster);
end
%% Low dimension visualization
%   - (1.1) true & cluster label from NMF 
%	- (1.2) CPI
[~,X] = pca(M); %PCA  
%X = tsne(M); %tSNE
%set the folder to save figures
folder ='SCC_plot';
low_dim_plot(X,C,cpi,true_label,folder,'PCA');
%% cluster-cluster graph 
%plot cluster-cluster relation
[cluster_relation,cr] = data_cluster_plot(H,folder,'cluster-cluster');

%% visulation of cells (no regularized)
[cluster_location,cell_location] = cell_visu(data,cr,H,C,true_label,cpi,folder);

%% check the TCs and infer lineage
cut = 0.34; %threshold of CPI to find TC
index_TC = find(cpi>cut);
C_TC = C;
C_TC(index_TC) = No_cluster + 1;
% plot # in and out cells of a state/total tc
state_TC(H,No_cluster,index_TC,C)
% plot state-state TC%
state_state_TC(H,No_cluster,index_TC)
% infer transition trajectories
cluster_order1=[2,1,3,4];
cluster_order2 = [2,1,4];

%% PRE
tic
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true);%,'MaxIterations',100,'OptimalityTolerance',1e-4);
ff = makefcn(cluster_location,H_normailed,0.1/size(H_normailed,1));
x = fminunc(ff,cell_location,options); %locations of cells
toc
%plot: color of clusters based on the order of cluster now
cell_visu_PRE(cr,x,C,true_label,cpi,cluster_order1,folder,'PRE')
%% PRE (weight)
tic
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'MaxIterations',100,'OptimalityTolerance',1e-4);
[~,r] = pca(data);
ff = makefcn_weight(cluster_location,H_normailed,0.01/size(H_normailed,1),pdist(r)./max(pdist(r)));
x_weight = fminunc(ff,cell_location,options); %locations of cells
toc
cell_visu_PRE(cr,x_weight,C,true_label,cpi,cluster_order1,folder,'PRE_w')

%% order cells along possible lineages
[ordered_cell1,ordered_cell2] = order_cell(H_normailed,C,cpi,cluster_order1,cluster_order2,0.96);

%% get gene W from each lineage
lam = 10; %weight of regularization
%get W from the first transtion trajectory
[W1,W1_n,C_o] = Gene_cluster_W(data,H,M,ordered_cell1,cluster_order1,lam);
%compare the new cluster from NMF and the original one
Cal_NMI(C(sort(ordered_cell1)), C_o)


%get W from the second transtion trajectory
[W2,W2_n,C_o] = Gene_cluster_W(data,H,M,ordered_cell2,cluster_order2,lam);
%compare the new cluster from NMF and the original one
Cal_NMI(C(sort(ordered_cell2)), C_o)

%% visulation of genes
%find marker genes of each cluster
[marker1,W1_marker] = marker_gene(W1_n);
% plot violinplot of marker genes in each cluster
violinplot_gene(data(sort(ordered_cell1),:),C(sort(ordered_cell1)),data,C,W1_marker,gene_name,folder)

%%
[marker2,W2_marker] = marker_gene(W2_n);
% plot violinplot of marker genes in each cluster
violinplot_gene(data(sort(ordered_cell2),:),C(sort(ordered_cell2)),data,C,W2_marker,gene_name,folder)
%% get the transition genes of each transitions
%compute pseudotime along each lineage
[pseudotime1] = local_pseudotime(x,ordered_cell1);
cut_corr = 0.64;
%detect transition genes of each transition
[transition1] = transition_gene(C,data,cluster_order1,ordered_cell1,pseudotime1,cut_corr);

[pseudotime2] = local_pseudotime(x,ordered_cell2);
[transition2] = transition_gene(C,data,cluster_order2,ordered_cell2,pseudotime2,cut_corr);

%% Heatmap of genes
%plot heatmap of top 4 marker genes of each cluster and top 3*2 transition
%genes of each transition
[rw_gene1,rw_gene_index1,gene_cluster1] = cell_gene_forrowave7(cluster_order1,transition1,marker1,gene_name,4,3);
[fitdata_ptime1] = Rowave_ptime_given_genes_new(data',gene_name,rw_gene1,ordered_cell1,gene_cluster1,folder);

%W only has 3 columns here, ordered cluster is [2 1 3] now
[rw_gene2,rw_gene_index2,gene_cluster2] = cell_gene_forrowave7([2 1 3],transition2,marker2,gene_name,4,3);
[fitdata_ptime2] = Rowave_ptime_given_genes_new(data',gene_name,rw_gene2,ordered_cell2,gene_cluster2,folder);

%% plot transition genes along each lineage
for i = 1:length(rw_gene_index1)
    [pseudotime1] = gene_pseudotime(x,ordered_cell1,cluster_order1,rw_gene_index1(i),rw_gene1{i},data,C,6,folder,'pseudotime1');
end
for i = 1:length(rw_gene_index2)
    [pseudotime2] = gene_pseudotime(x,ordered_cell2,cluster_order2,rw_gene_index2(i),rw_gene2{i},data,C,6,folder,'pseudotime2');
end

%% visualiztion of cell with global pseudotime
branching_cluster =1; % the cluster when two transition trajecttories branch
[ordered_cell,pseudotime] = cell_visu_global_p(x_weight,cr,C,pseudotime1,pseudotime2,ordered_cell1,ordered_cell2,branching_cluster,folder);

%% function
function fcn = makefcn(a,p,lam)

fcn = @visu2;
function [f,g] = visu2(x)
% Calculate objective f
f = 0;
for j = 1:size(p,2)
    f = f + sum(p(:,j).*diag((x-ones(size(x,1),1).*a(j,:))*(x-ones(size(x,1),1).*a(j,:))'));
end
f = f - lam*sum(pdist(x,'squaredeuclidean'));

if nargout > 1 % gradient required
    g = zeros(size(x,1),2);
for j = 1:size(p,2)
    g = g + 2*p(:,j).*(x-ones(size(x,1),1).*a(j,:));
end
g = g - 2*lam*(size(x,1)*x-ones(size(x,1),1).*sum(x));

end
end
end






function fcn = makefcn_weight(a,p,lam,r)

fcn = @visu2_weight;
function [f,g] = visu2_weight(x)
% Calculate objective f
f = 0;
for j = 1:size(p,2)
    f = f + sum(p(:,j).*diag((x-ones(size(x,1),1).*a(j,:))*(x-ones(size(x,1),1).*a(j,:))'));
end
f = f - lam*sum(r*pdist(x,'squaredeuclidean')');

if nargout > 1 % gradient required
    g = zeros(size(x,1),2);

for j = 1:size(p,2)
    g = g + 2*p(:,j).*(x-ones(size(x,1),1).*a(j,:));
end
g = g - 2*lam*(sum(squareform(r),2).*x-squareform(r)*x);

end
end
end
