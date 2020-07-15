function [pseudotime1] = local_pseudotime(cell_location,ordered_cell1)
%based on lineage and cell location to get pseudotime
%computer the Euclidean dist of neighboring cells
dist2 = zeros(1,length(ordered_cell1)-1);
for i = 1:length(dist2)
    dist2(i) = norm(cell_location(ordered_cell1(i+1),:)-cell_location(ordered_cell1(i),:));
end
pseudotime1 = [0,cumsum(dist2)];
end

