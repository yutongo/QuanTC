function [ ] = low_dim_plot(X,C,cpi,true_label,folder,method)
% plot cluster on 2-dimensional space
% input
% - X: low_dim data
% - C: clustering result
% - H: H from sym nmf
% - method: dimension reduction method, 'tsne' or 'pca'
%
% Output
%   - (1.1) true & cluster label from NMF 
%	- (1.2) CPI


zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
ss = 60;
No_cluster = length(unique(C));
Marker_type = ['o','p','s','d','+'];





subplot(1,2,1)


for i = 1:No_cluster
    for j = 1:length(unique(true_label))
        index1 = find(C==i & true_label==j);
        if Marker_type(1)=='p'
            scatter(X(index1,1),X(index1,2),ss*2,mycolor(i,:),'filled',Marker_type(j),'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
            hold on
        else
            scatter(X(index1,1),X(index1,2),ss,mycolor(i,:),'filled',Marker_type(j),'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
            hold on
        end
    end
end


box off;
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('cluster from NMF')
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
hold off




subplot(1,2,2)

for i = 1:length(unique(true_label))
    if Marker_type(1)=='p'
        scatter(X(find(true_label==i),1),X(find(true_label==i),2),2*ss,cpi(find(true_label==i)),'filled',Marker_type(i));
        hold on
    else
        scatter(X(find(true_label==i),1),X(find(true_label==i),2),ss,cpi(find(true_label==i)),'filled',Marker_type(i));
        hold on
    end
end
%scatter(X(:,1),X(:,2),ss,diff,'filled');
colorbar
title('CPI')
set(gca,'xtick',[]);
set(gca,'ytick',[]);




print([folder '\' method],'-dpdf','-r300','-fillpage'); %'-dpdf',
end

