function [ordered_cell,pseudotime_ordered] = global_pseudotime(result,path_index,path,ordered_cell)
% Get global pseudotime of the top2 trajectories
% result{1} = H;
C = result{2};
%cpi = result{3};
% result{4} = cluster_location;
% result{5} = cell_location;
cell_location= result{6};
% result{7} = C_TC;

[pseudotime1] = local_pseudotime(cell_location,ordered_cell{path_index(1)});
[pseudotime2] = local_pseudotime(cell_location,ordered_cell{path_index(2)});
cluster_order1 = path{path_index(1)};
cluster_order2 = path{path_index(2)};

zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
cell_used = union(ordered_cell{path_index(1)},ordered_cell{path_index(2)});
p1=pseudotime1;
p2=pseudotime2;

aaa=find(path{path_index(1)}(1:length(path{path_index(2)}))-path{path_index(2)}==0);
branching_cluster = aaa(end);
%get global pseudotime
%get branching cell (last same cell before branching)
cell_repeated = intersect(ordered_cell{path_index(1)},ordered_cell{path_index(2)},'stable');
branching_cells = cell_repeated(C(cell_repeated)==branching_cluster);
branching_cell_index = branching_cells(end);
end_cell_index1 = ordered_cell{path_index(1)}(end);
end_cell_index2 = ordered_cell{path_index(2)}(end);

ordered_cell1 = ordered_cell{path_index(1)};
ordered_cell2 = ordered_cell{path_index(2)};

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



figure
colormap(cool)
scatter(cell_location(ordered_cell,1),cell_location(ordered_cell,2),30,pseudotime_ordered,'filled','o','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);

set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Pseudotime')
box on


%print([folder  '\'  folder1  '\' gene_name],'-dpdf','-r600','-fillpage');
end

