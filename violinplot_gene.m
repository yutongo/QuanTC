function [] = violinplot_gene(data,C,data_new,C_new,W_new,gene_name,folder)
No_cluster = size(W_new,2);
cluster_name = cell(1,No_cluster);
for i = 1:No_cluster
    cluster_name{i} = ['C' num2str(i)];
end
for i = 1:No_cluster
    gene_id = find(W_new(:,i)~=0);
    gene_use=[];
    for j = 1:length(gene_id)
        if length(find(W_new(gene_id(j),:)~=0))==1
            gene_use=[gene_use,gene_id(j)];
        end
    end
    data_violinplot_tr_gene(data,No_cluster,gene_use',gene_name,C,folder,strcat(cluster_name{i},' genes'))
    %data_violinplot_tr_gene(data_new,gene_id,gene_name,C_new,folder,strcat(cluster_name{i},' genes_cut'))
end
end

