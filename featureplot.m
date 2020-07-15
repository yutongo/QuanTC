function [ ] = featureplot(result,prodata,progene_name,gene)

cell_location= result{6};



ss=30;


[~,~,gene_index] = intersect(upper(gene),upper(progene_name),'stable');
paper_gene = progene_name(gene_index);
for i = 1:length(gene_index)
    figure
    scatter(cell_location(:,1),cell_location(:,2),ss,prodata(:,gene_index),'filled','o','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);

    cmap = flipud(colormap(hot));
        cmap(1,:) = [0.9 0.9 0.9];
        colormap(cmap)
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    title(paper_gene{i})
    box on;
end



end
