function [marker,W_marker] = marker_gene(W)
% input: gene*cluster matrix
 W_n = W./sum(W,2);
 W_marker = zeros(size(W_n));
 marker = cell(size(W,2),1);
 [~,C_g] = max(W_n,[],2); %cluster of genes
 for i = 1 : size(W,2)
     index = find(C_g==i);
     [~,I] = sort(W_n(index,i),'descend');
     remove_index=[];
     if ~isempty(I)
         for j=1:length(I)
             a = W_n(index(I(j)),:);
             b = maxk(a,2);
             if b(1)-b(2)<0.03
                 remove_index=[remove_index,j];
             end
         end
     end
     I(remove_index)=[];
     marker{i} = index(I);
     if length(I)>=9
        W_marker(marker{i}(1:9),i) = W_n(marker{i}(1:9),i);
     else
         W_marker(marker{i},i) = W_n(marker{i},i);
     end
end

