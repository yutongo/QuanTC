function [] = CPI_o_cell(cluster_order1,ordered_cell1,C,cpi)
% plot cells ordered along one transition trajectory wirh cpi value
for i = 1:length(cluster_order1)
    cell_o_c{i} = find(C(ordered_cell1)==cluster_order1(i));
end

zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(8,:) = [0.5 1 0];

if length(cluster_order1)==3
    mycolor(3,:)=[];
end

for j = 1:length(ordered_cell1)
        line([j,j],[0,cpi(ordered_cell1(j))],'Color',mycolor(cluster_order1==C(ordered_cell1(j)),:),'LineWidth',1)
    hold on
end
xlim([0 length(ordered_cell1)])

box off
%axis off
end

