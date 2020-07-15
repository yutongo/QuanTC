function [] = cell_visu_PRE(result, true_label, label_legend)
% result{1} = H;
C = result{2};
cpi = result{3};
% result{4} = cluster_location;
% result{5} = cell_location;
cell_location= result{6};
% result{7} = C_TC;



%PRE visualization


%get the locaion of cluster centers
No_cluster = length(unique(C));


zzz = get(gca,'colororder');
mycolor = zeros(13,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
mycolor(12,:) =     [0, 0.75, 0.75];
mycolor(13,:) = [0, 0.5, 0];




ss=30;
 



%plot clustering
for i = 1:No_cluster
        index1 = find(C==i);
        scatter(cell_location(index1,1),cell_location(index1,2),ss,mycolor(i,:),'filled','o','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
        hold on
end
hold off
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('Clustering')
cluster_name = cell(1,No_cluster);
for i = 1:No_cluster
    cluster_name{i} = ['C' num2str(i)];
end
legend(cluster_name,'FontSize',12,'Location','best')
box on



% plot true_label
figure
for i = 1:length(unique(true_label))
        index1 = find(true_label==i);
        scatter(cell_location(index1,1),cell_location(index1,2),ss,mycolor(i,:),'filled','o','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
        hold on
end
hold off
set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('True labels')
legend(label_legend,'FontSize',12,'Location','best')
box on



% plot CPI
figure

scatter(cell_location(:,1),cell_location(:,2),ss,cpi,'filled','o','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);

set(gca,'xtick',[]);
set(gca,'ytick',[]);
title('CPI labels')
box on






%plot CPI
% 
% ss=30;
%  for i = 1:No_cluster
%     for j = 1:length(unique(true_label))
%         index1 = find(C==i & true_label==j);
%             
%             if Marker_type(j)=='p'
%                 scatter(cell_location(index1,1),cell_location(index1,2),ss*2,cpi(index1),'filled',Marker_type(j),'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
%                 hold on
%             else
%                 scatter(cell_location(index1,1),cell_location(index1,2),ss,cpi(index1),'filled',Marker_type(j),'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
%                 hold on
%             end
%     end
%  end
%  colorbar
% xlim([-1.1 1.1])
% ylim([-1.15 1.1])
% hold off
% 
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% 
% box on
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% 
% a=5; b=3.8;
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, a, b], 'PaperUnits', 'Inches', 'PaperSize', [a, b])
%     fig.OuterPosition=fig.InnerPosition;
% 
% print([folder '\cell_visu_CPI'],'-dpdf','-r300','-fillpage');
end

