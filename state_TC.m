function [] = state_TC(H,No_cluster,index_TC,C)
% plot # in and out cells of a state/total tc
[I_top2,~,~] = toptwo(H);
p_TC = zeros(No_cluster);
 for i = 1:No_cluster
     a = union(find(I_top2(index_TC,1)==i),find(I_top2(index_TC,2)==i)); 
     %TC around cluster i
     for j = 1:No_cluster
        p_TC(i,j) = length(find(C(index_TC(a))==j))/length(index_TC);
     end
 end
 figure;
b = bar(p_TC,'stacked');
zzz = get(gca,'colororder');
mycolor = zeros(11,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
for i =1:No_cluster
    b(i).FaceColor = mycolor(i,:);
    b(i).FaceAlpha = 0.8;
    b(i).EdgeAlpha = 0;
end
yticks([0 0.3 0.6 0.9])
box off
end

 