function [cluster_location,cell_location] = cell_visu(data,cr,H,C,true_label,cpi,folder)
%Input data cr(cluster-cluster relation)
%Output visualization of cells in 2d


%plot potential transiton trajectory
G = digraph(cr,'OmitSelfLoops');
TR = shortestpathtree(G,2);
h = plot(G);
h.LineStyle = 'none';
h.NodeLabel =[];
h.MarkerSize=0.01;
highlight(h,TR,'EdgeColor',[1 1 1],'LineStyle','-','LineWidth',1)
hold on


%get the locaion of cluster centers
No_cluster = length(unique(C));
cluster_location = zeros(No_cluster,2);
for i = 1:No_cluster
    cluster_location(i,1) = h.XData(i);
    cluster_location(i,2) = h.YData(i);
end

zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];


cell_location = zeros(size(data,1),2);
for i = 1:size(data,1)
    cell_location(i,1) = H(i,:)*cluster_location(:,1);
    cell_location(i,2) = H(i,:)*cluster_location(:,2);
end



ss=50;

Marker_type = ['o','p','s','d','+'];




%plot clustering
for i = 1:No_cluster
    for j = 1:length(unique(true_label))
        index1 = find(C==i & true_label==j);
        if Marker_type(1)=='p'
            scatter(cell_location(index1,1),cell_location(index1,2),ss*2,mycolor(i,:),'filled',Marker_type(j),'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
            hold on
        else
            scatter(cell_location(index1,1),cell_location(index1,2),ss,mycolor(i,:),'filled',Marker_type(j),'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
            hold on
        end
    end
end


 
 
 

hold off
set(gca,'xtick',[]);
set(gca,'ytick',[]);
xlim([-1.1 1.1])
ylim([-1.15 1.1])


%legend boxoff 
set(gca,'position',[0 0 1 1],'units','normalized')
box on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

    fig.Units = 'Inches';
    fig.Position = [0 0 3 2.3];
    fig.OuterPosition=fig.InnerPosition;
%title('color by truelabel')
print([folder '\cell_visu_cluster'],'-dpdf','-r300','-fillpage');


%plot CPI
figure
ss=30;
 for i = 1:No_cluster
    for j = 1:length(unique(true_label))
        index1 = find(C==i & true_label==j);
            
            if Marker_type(j)=='p'
                scatter(cell_location(index1,1),cell_location(index1,2),ss*2,cpi(index1),'filled',Marker_type(j),'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
                hold on
            else
                scatter(cell_location(index1,1),cell_location(index1,2),ss,cpi(index1),'filled',Marker_type(j),'MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
                hold on
            end
    end
 end
 colorbar
xlim([-1.1 1.1])
ylim([-1.15 1.1])
hold off

set(gca,'xtick',[]);
set(gca,'ytick',[]);

box on
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

    fig.Units = 'Inches';
    fig.Position = [0 0 3 2.3];
    fig.OuterPosition=fig.InnerPosition;

print([folder '\cell_visu_CPI'],'-dpdf','-r300','-fillpage');

end

