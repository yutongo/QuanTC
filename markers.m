function [marker1,transition1 ] = markers(result,prodata,M,cluster_order1,ordered_cell1,lam)

H = result{1};
C = result{2};
% result{3} = cpi;
% result{4} = cluster_location;
% result{5} = cell_location;
x = result{6};
% result{7} = C_TC;


%first lineage
data_n = prodata(sort(ordered_cell1),:)./sum(prodata(sort(ordered_cell1),:),1); %normalize each gene wrt cells
H_n = H(sort(ordered_cell1),:)./sum(H(sort(ordered_cell1),:),2); %normalize cluster sum
W_0 = (data_n'*M(sort(ordered_cell1),sort(ordered_cell1))*H_n)/(H_n'*H_n);
W_0(W_0<0)=0;
W_0=W_0';
[H_new1,W_new1,~] = nmf_knownH(data_n,length(cluster_order1),H_n,lam,H_n,W_0);%nmf_knownH(data_n,3,H_n,lam,'h0', H_n,'w0',W_0,'alg','als');
W_new1 = W_new1';
W_new1_n = W_new1./sum(W_new1,2);%max(W_new,[],1); %sum(W_new);
[~,C_o] = max(H_new1,[],2);
Cal_NMI(C(sort(ordered_cell1)), C_o)
display(['Consistency of clustering: ' num2str(Cal_NMI(C(sort(ordered_cell1)), C_o))])
[marker1,~] = marker_gene(W_new1_n);
[pseudotime1] = local_pseudotime(x,ordered_cell1);
cut_corr = 0.64;%0.5;
[transition1] = transition_gene(C,prodata,cluster_order1,ordered_cell1,pseudotime1,cut_corr);
end