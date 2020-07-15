function [ ] = tran_gene_plot(result,data,gene_index,cluster_order1,ordered_cell1,gene_name)


C = result{2};

cell_location = result{6};


mm=6;
zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
%based on lineage and cell location to get pseudotime
%computer the Euclidean dist of neighboring cells
dist2 = zeros(1,length(ordered_cell1)-1);
for i = 1:length(dist2)
    dist2(i) = norm(cell_location(ordered_cell1(i+1),:)-cell_location(ordered_cell1(i),:));
end
pseudotime1 = [0,cumsum(dist2)];
data_used = data(ordered_cell1,gene_index);

x = 0:1/(length(ordered_cell1)-1):1; % evenly distributed
%x = pseudotime1;
y = data_used';
z = C(ordered_cell1);
idx_b1 = [];
for i = 1:length(cluster_order1)
    idx_b1 = [idx_b1;find(z==cluster_order1(i))];
end
y_b1 = y(idx_b1);

p_b1 = polyfit(x(sort(idx_b1)),y_b1,mm);
f_b1 = polyval(p_b1,x(sort(idx_b1)));

% plot gene expression along ptime (label)
for ik = 1:length(cluster_order1)
    if length(cluster_order1)~=3 || ik~=3
    scatter(x(find(z==cluster_order1(ik))),y(find(z==cluster_order1(ik))),20,mycolor(ik,:),'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
    else
        scatter(x(find(z==cluster_order1(ik))),y(find(z==cluster_order1(ik))),20,mycolor(4,:),'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
    end
end


hold on;
plot(x(sort(idx_b1)),f_b1,'LineWidth',3,'Color','k','MarkerSize',10);
%xticks([])
L = get(gca,'YLim');
set(gca,'YTick',round(linspace(L(1),L(2),2),1))

yticks([ceil(min(y)),floor(max(y))])

xticks([]);
hold off
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
    fig.Units = 'Inches';
    fig.Position = [0 0 3 2.3];
    fig.OuterPosition=fig.InnerPosition;


end
