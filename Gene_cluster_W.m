function [W_new1,W_new1_n,C_o] = Gene_cluster_W(data,H,M,ordered_cell1,cluster_order1,lam)
%first lineage
data_n = data(sort(ordered_cell1),:)./sum(data(sort(ordered_cell1),:),1); %normalize each gene wrt cells
H_n = H(sort(ordered_cell1),unique(cluster_order1))./sum(H(sort(ordered_cell1),unique(cluster_order1)),2); %normalize cluster sum
W_0 = (data_n'*M(sort(ordered_cell1),sort(ordered_cell1))*H_n)/(H_n'*H_n);
W_0(W_0<0)=0;
W_0=W_0';
%NMF
[H_new1,W_new1,~] = nmf_knownH(data_n,length(cluster_order1),H_n,lam,H_n,W_0);
W_new1 = W_new1';
W_new1_n = W_new1./sum(W_new1,2);
[~,C_o] = max(H_new1,[],2);
end

