function  [cluster_relation,cr] = data_cluster_plot(H,folder,method)
%get likehoods of cell belonging to clusters
H_new = H./sum(H,2);
cluster_relation = H_new'*H_new;

cr = cluster_relation;
%ignore small value
cr(cr==min(min(cr)))=0;
if size(H,2)>3
    
cr(cr==min(min(cr(cr~=0))))=0;
end
figure;
G = graph(cr,'OmitSelfLoops');
TR = shortestpathtree(G,2);
h = plot(G);
highlight(h,TR,'EdgeColor','r')


print([folder '\' method],'-dpdf','-r300','-fillpage');
end

