function [tg1_index,corr_gene,var_gene] = transition_gene(C,data,cluster_order1,ordered_cell1,pseudotime1,cut_corr)
% check all the corr choose the highest/lowest 
% transition genes
tg1_index = cell(length(cluster_order1)-1,2);
%cut_corr=0.64;
%[~,C_w] = max(W_new1,[],2);
corr_gene = zeros(size(data,2),length(cluster_order1)-1);
var_gene = zeros(size(data,2),1);
for i = 1:length(cluster_order1)-1
    cell_used = find(C(ordered_cell1)==cluster_order1(i) | C(ordered_cell1)==cluster_order1(i+1));
    if i==1 || i==length(cluster_order1)-1 % not two parts
        cell_notused = find(C(ordered_cell1)~=cluster_order1(i) & C(ordered_cell1)~=cluster_order1(i+1));
        for j = 1:size(data,2)
            var_gene(j) = var(data(ordered_cell1(cell_notused),j));
        end
    else
        a = find(C(ordered_cell1)==cluster_order1(i),1)-1; %previous
        if length(cluster_order1)<4
            for j = 1:size(data,2)
                var_gene(j) = var(data(ordered_cell1(1:a),j));
            end
        else
            b = find(C(ordered_cell1)==cluster_order1(i+2),1);
            for j = 1:size(data,2)
                var_gene(j) = var(data(ordered_cell1(1:a),j))+var(data(ordered_cell1(b,end),j));
            end
        end
    end
    for j = 1:size(data,2)
        corr_gene(j,i) = corr(pseudotime1(cell_used)',data(ordered_cell1(cell_used),j),'Type','Spearman');
    end
end
[~,Co_max] = max(abs(corr_gene),[],2);
for i = 1:length(cluster_order1)-1
    %p11 = intersect(find((corr_gene(:,i)>cut_corr) & (C_w'==cluster_order1(i+1))), find(Co_max==i)');
    p1_score = corr_gene(:,i);
    p11 = intersect(find((p1_score>cut_corr)), find(Co_max==i)');
    %decreasing
    %p12 = intersect(find(corr_gene(:,i)<-cut_corr & C_w'==cluster_order1(i)), find(Co_max==i)');
    p12 = intersect(find(p1_score<-cut_corr ), find(Co_max==i)');
    [~,I1]=sort(p1_score(p11),'descend');
    [~,I2]=sort(p1_score(p12),'ascend');
    tg1_index{i,1} = p11(I1);
    tg1_index{i,2} = p12(I2);
end

end

