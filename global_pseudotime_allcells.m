function [pseudotime_allcells] = global_pseudotime_allcells(result,path_index,path,ordered_cell)
% Get global pseudotime of all the cells based on one trajectory
H = result{1};
C = result{2};
%cpi = result{3};
% result{4} = cluster_location;
% result{5} = cell_location;
cell_location= result{6};
% result{7} = C_TC;

H_norm = H./sum(H,2);


[pseudotime1] = local_pseudotime(cell_location,ordered_cell{path_index});
cluster_order1 = path{path_index};


pseudotime_allcells = zeros(1,length(C));

pseudotime_allcells(ordered_cell{path_index})= pseudotime1;
% get the pseudotime of cells not in the chosen trajectory based on the 
% cluster centers

remaining_cell_index = setdiff(1:length(C),ordered_cell{path_index});

% get the cluster center index
[~,I] = maxk(H_norm(ordered_cell{path_index},:),1,1);


for i = 1:length(remaining_cell_index)
    
    pseudotime_allcells(remaining_cell_index(i)) = dot(H_norm(remaining_cell_index(i),:),pseudotime1(I));
end
    
%figure
%colormap(cool)
%scatter(cell_location(ordered_cell{path_index},1),cell_location(ordered_cell{path_index},2),30,pseudotime_allcells(ordered_cell{path_index}),'filled','o','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);


figure
colormap(cool)
scatter(cell_location(:,1),cell_location(:,2),30,pseudotime_allcells,'filled','o','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);

set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Pseudotime allcells')
box on



end

