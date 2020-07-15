function [ordered_cell,pseudotime_ordered] = global_pseudotime(data,C,branching_cluster,pseudotime1,pseudotime2,ordered_cell1,ordered_cell2,cluster_order1,cluster_order2,gene_index,gene_name,mm,folder1,folder)
zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
cell_used = union(ordered_cell1,ordered_cell2);
p1=pseudotime1;
p2=pseudotime2;

%get global pseudotime
%get branching cell (last same cell before branching)
cell_repeated = intersect(ordered_cell1,ordered_cell2,'stable');
branching_cells = cell_repeated(C(cell_repeated)==branching_cluster);
branching_cell_index = branching_cells(end);
end_cell_index1 = ordered_cell1(end);
end_cell_index2 = ordered_cell2(end);

pseudotime = zeros(1,length(cell_used));
if pseudotime2(end)>= pseudotime1(end)
    pseudotime(1:find(ordered_cell1==branching_cell_index)) = pseudotime2(1:find(ordered_cell1==branching_cell_index));
else
    pseudotime(1:find(ordered_cell2==branching_cell_index)) = pseudotime1(1:find(ordered_cell2==branching_cell_index));
end
% if both lineage converge to the same cluster, map the shorter
% pseudotime to longer ones
if C(end_cell_index1)==C(end_cell_index2)
    if pseudotime2(end)>= pseudotime1(end)
        a=find(ordered_cell1==branching_cell_index);
        pseudotime1(a:end) = (pseudotime1(a:end)-pseudotime1(a))*(pseudotime2(end)-pseudotime2(a))/(pseudotime1(end)-pseudotime1(a))+pseudotime1(a);
    else
        a=find(ordered_cell2==branching_cell_index);
        pseudotime2(a:end) = (pseudotime2(a:end)-pseudotime2(a))*(pseudotime1(end)-pseudotime1(a))/(pseudotime2(end)-pseudotime2(a))+pseudotime2(a);

    end
end



for i = 1:length(cell_used)
    if ismember(cell_used(i),ordered_cell1) && ismember(cell_used(i),ordered_cell2)
        pseudotime(i) = min(pseudotime1(ordered_cell1==cell_used(i)),pseudotime2(ordered_cell2==cell_used(i)));
    elseif ~ismember(cell_used(i),ordered_cell1) && ismember(cell_used(i),ordered_cell2)
         pseudotime(i) = pseudotime2(ordered_cell2==cell_used(i));
    elseif ismember(cell_used(i),ordered_cell1) && ~ismember(cell_used(i),ordered_cell2)
         pseudotime(i) = pseudotime1(ordered_cell1==cell_used(i));
    end
end
%reorder the global pseudotime
[pseudotime_ordered,I] = sort(pseudotime);
ordered_cell = cell_used(I);
data_used = data(ordered_cell,gene_index);
x = pseudotime_ordered;
y = data_used';
z = C(ordered_cell); 
%find index of ordered_cell1 in ordered_cell
index1 = zeros(length(ordered_cell1),1);
for i = 1:length(ordered_cell1)
    index1(i) = find(ordered_cell==ordered_cell1(i));
end
%find index of ordered_cell1 in ordered_cell
index2 = zeros(length(ordered_cell2),1);
for i = 1:length(ordered_cell2)
    index2(i) = find(ordered_cell==ordered_cell2(i));
end

z1 = zeros(length(ordered_cell),1);
z1(index1) = C(ordered_cell(index1)); 
idx_b1 = [];
for i = 1:length(cluster_order1)
    idx_b1 = [idx_b1;find(z1==cluster_order1(i))];
end
y_b1 = y(idx_b1);
p_b1 = polyfit(x(sort(idx_b1)),y_b1,mm);
f_b1 = polyval(p_b1,x(sort(idx_b1)));

z2 = zeros(length(ordered_cell),1);
z2(index2) = C(ordered_cell(index2)); 
idx_b2 = [];
for i = 1:length(cluster_order2)
    idx_b2 = [idx_b2;find(z2==cluster_order2(i))];        
end
y_b2 = y(idx_b2);
p_b2 = polyfit(x(sort(idx_b2)),y_b2,mm);
f_b2 = polyval(p_b2,x(sort(idx_b2)));

% plot gene expression along ptime (label)

for ik = 1:4
    scatter(x(find(z==ik)),y(find(z==ik)),20,mycolor((ik),:),'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
end


hold on;
plot(x(sort(idx_b1)),f_b1,'-','LineWidth',2,'Color','k','MarkerSize',10);
hold on;
plot(x(sort(idx_b2)),f_b2,':','LineWidth',2,'Color','k','MarkerSize',10);
xticks([0,1]);
hold off
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 3];
%
    fig.Units = 'Inches';
    fig.Position = [0 0 3 2.3];
    fig.OuterPosition=fig.InnerPosition;
print([folder  '\'  folder1  '\' gene_name],'-dpdf','-r600','-fillpage');
end

