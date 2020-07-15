function [ ] = state_state_TC(H,No_cluster,index_TC)
[I_top2,H_t1,H_t2] = toptwo(H);
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


number_transiton_p_new = zeros(No_cluster*(No_cluster-1)/2,No_cluster);
number_transiton_p= number_transiton/size(H,1);
c_order=zeros(No_cluster*(No_cluster-1)/2,2);
%get all possible state-state transition
a=1;
    for ii = 1:No_cluster
            
        for jj = ii+1:No_cluster
            c_order(a,1) = ii;
            c_order(a,2) = jj;
            a=a+1;
        end
    end


for i = 1:No_cluster*(No_cluster-1)/2
        number_transiton_p_new(i,c_order(i,1)) = number_transiton_p(c_order(i,1),c_order(i,2));
        number_transiton_p_new(i,c_order(i,2)) = number_transiton_p(c_order(i,2),c_order(i,1));
end
figure;
b = bar(number_transiton_p_new,'stacked');
box off
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
box off
end

