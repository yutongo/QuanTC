function [rw_gene_index1,gene_cluster1 ] = heatmap(data,marker1,transition1,cluster_order1,ordered_cell1,gene_name,a,b)

[rw_gene1,rw_gene_index1,gene_cluster1] = cell_gene_forrowave7(cluster_order1,transition1,marker1,gene_name,a,b);
[fitdata_ptime1] = Rowave_ptime_given_genes_new(data',gene_name,rw_gene1,ordered_cell1,gene_cluster1);

end
