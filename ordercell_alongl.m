function [ordered_cell] = ordercell_alongl(H_normailed,H_t1,H_t2,I_top2,cpi,cluster_order,k)
%starting cluster
a=0.63;
cell_used =find(H_t1+H_t2>a);
I = intersect(intersect(find(I_top2(:,1)==cluster_order(1)),find(I_top2(:,2)==cluster_order(2))),cell_used);
[~,index1] = sort(cpi(I),'ascend');
a = setdiff(find(H_normailed(:,cluster_order(1))>=k),I);
[~,ia] = sort(cpi(a),'ascend');
ordered_cell = [a(ia);I(index1)];
% exclude the cells beyong the lineage
for i = 2:length(cluster_order)-1
    I1 = find(I_top2(:,1)==cluster_order(i));
    I0 = intersect(intersect(I1,find(I_top2(:,2)==cluster_order(i-1))),cell_used);%interaction previous cluster
    I2 = intersect(intersect(I1,find(I_top2(:,2)==cluster_order(i+1))),cell_used);%interaction later cluster
    [~,I0_order] = sort(cpi(I0),'descend');
    [~,I2_order] = sort(cpi(I2),'ascend');
    ordered_cell = [ordered_cell;I0(I0_order);find(H_normailed(:,cluster_order(i))>=k);I2(I2_order)];
end
%last cluster
I1 = find(I_top2(:,1)==cluster_order(end));
I0 = intersect(intersect(I1,find(I_top2(:,2)==cluster_order(end-1))),cell_used);
[~,I0_order] = sort(cpi(I0),'descend');
%ordered_cell = [ordered_cell;I0(I0_order);find(H(:,cluster_order(end))==1)];
a = setdiff(find(H_normailed(:,cluster_order(end))>=k),I0);
[~,ia] = sort(cpi(a),'descend');
ordered_cell = [ordered_cell;I0(I0_order);a(ia)];
end

