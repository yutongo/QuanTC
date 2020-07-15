function [rw_gene,rw_gene_index,gene_cluster] = cell_gene_forrowave7(cluster_order,transition1,marker,gene_name,a,b)
% a: # of each marker
% b: # of each transition (acsending/descending)
rw_gene_index = [];gene_cluster=[];
for i = 1:length(cluster_order)
    if length(marker{cluster_order(i)})>=a
        rw_gene_index = [rw_gene_index,marker{cluster_order(i)}(1:a)'];
        gene_cluster = [gene_cluster,i*ones(1,a)];
    else
        rw_gene_index = [rw_gene_index,marker{cluster_order(i)}'];
        gene_cluster = [gene_cluster,i*ones(1,length(marker{cluster_order(i)}))];
    end
    if i<length(cluster_order)
        for j = 2:-1:1
            if length(transition1{i,j})>=b
                rw_gene_index = [rw_gene_index,transition1{i,j}(1:b)'];
                gene_cluster = [gene_cluster,(i+length(cluster_order))*ones(1,b)];
            else
                rw_gene_index=[rw_gene_index,transition1{i,j}'];
                gene_cluster = [gene_cluster,(i+length(cluster_order))*ones(1,length(transition1{i,j}))];
            end
        end
    end
end


% %check if have duplicates
[~, ind] = unique(rw_gene_index);
duplicate_ind = setdiff(1:size(rw_gene_index, 2), ind);
duplicate_value = rw_gene_index(duplicate_ind);



if isempty(duplicate_value) %no duplicates
    rw_gene = gene_name(rw_gene_index);
else
    % keep trans genes pick new de genes
    %duplicate_ind show the second index (in clusters)
    index_cluster = [];u_i=[];
    for l = 1:length(duplicate_value)
        index_dul_all = find(rw_gene_index==duplicate_value(l));
        cluster_dul = gene_cluster(find(rw_gene_index==duplicate_value(l)));
        %cluster index of the dul genes
        index_cluster = [index_cluster,cluster_dul(cluster_dul>length(cluster_order))]; 
        %find the index of dup genes in all gene list
        u_i = [u_i,index_dul_all(cluster_dul==index_cluster(l))];
    end
%     index_cluster = gene_cluster(duplicate_ind);
%     u_i = unique(index_cluster);
%     k = 0;
    for i = 1:length(u_i)
        if ismember(u_i(i),transition1{index_cluster(i)-length(cluster_order),1})
            change_gene = transition1{index_cluster(i)-length(cluster_order),1}(b+1:end);
        else
            change_gene = transition1{index_cluster(i)-length(cluster_order),2}(b+1:end);
        end
        ii = change_gene(find(~ismember(change_gene,rw_gene_index),1));
        if ~isempty(ii)
            rw_gene_index(u_i(i)) = ii;
        else
            rw_gene_index(u_i(i)) = Inf;
            gene_cluster(u_i(i)) = Inf;
        end
    end
    rw_gene_index(rw_gene_index==Inf)=[];
    gene_cluster(gene_cluster==Inf) = [];
%     for i = 1:length(u_i)
%         for j = 1:length(find(index_cluster==u_i(i)))
%             ii = DE_gene{u_i(i)}(a+j+1);
%             if ismember(ii,rw_gene_index)
%                 ii = DE_gene{u_i(i)}(a+j+1);
%             end
%             rw_gene_index(duplicate_ind(k)) = ii(end);
%             k = k+1;
%         end
%     end
    rw_gene = gene_name(rw_gene_index);
end

end

