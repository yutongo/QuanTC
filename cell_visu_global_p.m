function [ordered_cell,pseudotime] = cell_visu_global_p(x,cr,C,pseudotime1,pseudotime2,ordered_cell1,ordered_cell2,branching_cluster,folder)
%get global pseudotime
%get branching cell (last same cell before branching)
cell_repeated = intersect(ordered_cell1,ordered_cell2,'stable');
branching_cells = cell_repeated(C(cell_repeated)==branching_cluster);
branching_cell_index = branching_cells(end);
end_cell_index1 = ordered_cell1(end);
end_cell_index2 = ordered_cell2(end);
cell_used = union(ordered_cell1,ordered_cell2);
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



G = digraph(cr,'OmitSelfLoops');
TR = shortestpathtree(G,2);
h = plot(G);
h.LineStyle = 'none';
h.NodeLabel =[];
h.MarkerSize=0.01;
highlight(h,TR,'EdgeColor',[1 1 1],'LineStyle','-','LineWidth',1)
hold on

zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
No_cluster = length(unique(C));
ss=30;
index1 = sort(ordered_cell);
colormap(cool)
scatter(x(ordered_cell,1),x(ordered_cell,2),ss,pseudotime_ordered,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
hold off

set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'position',[0 0 1 1],'units','normalized')
box on
fig = gcf;
fig.PaperPositionMode = 'auto';
a=5; b=3.8;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, a, b], 'PaperUnits', 'Inches', 'PaperSize', [a, b])
print([folder '\' 'pseudotime'],'-dpdf','-r300','-fillpage');
end

