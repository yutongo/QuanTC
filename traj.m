function [path,ordered_cell] = traj(result,start_cluster,No_cluster)
% compute the potential trajectory
H = result{1};
cluster_location = result{4};
C_TC = result{7};
cpi = result{3};

index_TC = find(C_TC==No_cluster+1);

[I_top2,H_t1,H_t2] = toptwo(H./sum(H,2));

%TC% between clusters
number_transiton = zeros(No_cluster);
for i = 1:No_cluster
    for j = 1:No_cluster
        kk =0;
        for k = 1:size(H,1)
            % only check with TC
            if isequal(I_top2(k,:),[i,j]) && ismember(k,index_TC) && H_t1(k)+H_t2(k)>0.63
                kk= kk+1;
            end
        end
        number_transiton(i,j) = kk;
    end
end
number_transiton2 = number_transiton+number_transiton.';
number_transiton_p = number_transiton2./length(C_TC); %tc% among all cells

% plot TC% between clusters on 2d
%h = plot(cluster_location(:,1),cluster_location(:,2),'.','MarkerSize',30);
h=scatter(cluster_location(:,1),cluster_location(:,2),100,'b','filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
hold on
dx = 0.01;
for i = 1:No_cluster
    
    text(cluster_location(i,1)+dx,cluster_location(i,2)+dx,['C' num2str(i)])
    for j = i+1:No_cluster
    text((cluster_location(i,1)+cluster_location(j,1))/2,(cluster_location(i,2)+cluster_location(j,2))/2,[sprintf('%.2f',(number_transiton_p(i,j))*100) '%'],'FontSize',11)
    plot([cluster_location(i,1),cluster_location(j,1)],[cluster_location(i,2),cluster_location(j,2)],'k')
    hold on
    end
end
hold off
axis off
box off


% given the start cluster, choose the least TC% to be the end cluster
a = number_transiton_p(start_cluster,:);
a(a == 0 ) = NaN;
[~,end_cluster] = min(a);


% find all the path from start cluster to end cluster and give cells
% involved 
path = pathbetweennodes(number_transiton_p, start_cluster, end_cluster);
ordered_cell = {};
if No_cluster == 4;
    for i = 1:size(path,1)
        ordered_cell{i} = ordercell_alongl(H./sum(H,2),H_t1,H_t2,I_top2,cpi,path{i},0.96);
    end
else
    for i = 1:size(path,1)
        ordered_cell{i} = ordercell_alongl(H./sum(H,2),H_t1,H_t2,I_top2,cpi,path{i},0.99955);
    end
end

for i = 1:size(path,1)
    display(['trajectory: ' num2str(path{i})  ', percentage of cells involved: ' num2str(length(ordered_cell{i})/length(C_TC))])
end

end