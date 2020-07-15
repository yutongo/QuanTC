function [ ] = pTC(result,No_cluster)
H = result{1};
C = result{2};
C_TC = result{7};

[I_top2,~,~] = toptwo(H);
index_TC =  find(C_TC == No_cluster+1);




% # in and out cells of a state/total tc
p_TC = zeros(No_cluster);
 for i = 1:No_cluster
     a = union(find(I_top2(index_TC,1)==i),find(I_top2(index_TC,2)==i)); 
     %TC around cluster i
     for j = 1:No_cluster
            p_TC(i,j) = length(find(C(index_TC(a))==j))/length(index_TC);
     end
 end
b = bar(p_TC,'stacked');
zzz = get(gca,'colororder');

mycolor = zeros(13,3);
mycolor(1:7,:) = zzz;
mycolor(8,:) = [0 0 0];
mycolor(9,:) = [0 0 0.803922];
mycolor(10,:) = [1 0 1];
mycolor(11,:) = [0.5 1 0];
mycolor(12,:) = [0, 0.75, 0.75];
mycolor(13,:) = [0, 0.5, 0];


for i =1:No_cluster
    b(i).FaceColor = mycolor(i,:);
    b(i).FaceAlpha = 0.8;
    b(i).EdgeAlpha = 0;
end
%yticks([0 0.1 0.2 0.3 0.4])
box off

xlabel('Clusters')
ylabel('TC%')

end

