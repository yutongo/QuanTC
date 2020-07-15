function [newdata,gene_short_name,procell_label] = preprocessing(data,gene_name,cell_label,min_exp,gene_selected)
% filter cells don't express >95% genes
cell_not_use = find(sum(data<1e-16,2)/size(data,2)>min_exp);
data(cell_not_use,:) = [];
cell_label(cell_not_use) =[];
procell_label = cell_label;
% find gene express more than 10% cells and var>5e-3 (gene with varied expressions)
index_gene_use1 =0;
for i = 1:size(data,2)
    a = data(:,i);
    if length(find(a<1e-15))/size(data,2)>0.9 || var(a)<5e-3
        continue
    else 
        index_gene_use1 = [index_gene_use1,i];
    end
end
index_gene_use1(1) =[];
% fit 3 mixture gaussian models
%remove mean 0 with >90% component
warning off;
index_gene_use2 = 0;
g_prob =0; g_mean =0; 
for i = 1:size(index_gene_use1,2)
    a = fitgmdist(data(:,index_gene_use1(i)),3,'RegularizationValue',0.01);
    if sum(a.ComponentProportion(find(a.mu<=1e-9)))<0.9
        [~,index2]=sort(a.ComponentProportion,'descend');
        g_prob = [g_prob,abs(diff(a.ComponentProportion(index2(1:2))))];
        g_mean = [g_mean,abs(diff(a.mu(index2(1:2))))];
        index_gene_use2 = [index_gene_use2,index_gene_use1(i)];
    end
end
index_gene_use2(1) =[];
g_prob(1) =[];
g_mean(1) =[];
% rank the informative genes:
%rank 1 denote max 2 prob diff from small to big; rank 2 denote mean diff from big to small
[~,rank1] = sort(g_prob,'ascend');
[~,rank2] = sort(g_mean,'descend');
rank = zeros(1,length(index_gene_use2));
for i = 1:length(index_gene_use2)
    rank(i) = find(rank1==i)+find(rank2==i);
end
% smaller rank means better
[~,index_gene_use2_rank] = sort(rank,'ascend');
index_gene_use_ranked = index_gene_use2(index_gene_use2_rank);

%output
newdata = data(:,index_gene_use_ranked(1:gene_selected));
gene_short_name = gene_name(index_gene_use_ranked(1:gene_selected));

clear data
clear gene_name
clear cell_label
end

