function  data_violinplot_tr_gene(data,No_cluster,tr_gene_id,gene_name,label,folder,method)
sw = find(label==1,1);
d = data(:,tr_gene_id);
if sw ~=1
    label([1 sw],:) = label([sw 1],:);
    d([1 sw],:) = d([sw 1],:);
end
figure;
for i = 1:No_cluster
    cluster_name{i} = ['C' num2str(i)];
end

    for i = 1:size(tr_gene_id,1)
        subplot(ceil(size(tr_gene_id,1)/3),3,i)
        myvionplot(d(:,i),label,length(unique(label)),cluster_name)
        title(gene_name(tr_gene_id(i)) )
    end

print([ folder '\' method],'-dpdf','-r300')
end