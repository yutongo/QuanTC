function [ordered_cell1,ordered_cell2] = order_cell(H,C,cpi,cluster_order1,cluster_order2,a)
[I_top2,H_t1,H_t2] = toptwo(H);
%get the order of cells along each transition trajcetory
[ordered_cell1] = ordercell_alongl(H,H_t1,H_t2,I_top2,cpi,cluster_order1,a);
CPI_o_cell(cluster_order1,ordered_cell1,C,cpi)
[ordered_cell2] = ordercell_alongl(H,H_t1,H_t2,I_top2,cpi,cluster_order2,a);
figure;
CPI_o_cell(cluster_order2,ordered_cell2,C,cpi)
end

