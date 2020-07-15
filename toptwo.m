function [I_top2,H_t1,H_t2] = toptwo(H)
%get the index of top two of each row
H_normailed =  H./sum(H,2);
I_top2 = zeros(size(H,1),2);
[H_t1,I_top2(:,1)] = max(H_normailed,[],2);
HH = H_normailed; HH(bsxfun(@eq, H, max(H, [], 2))) = -Inf;
[H_t2,I_top2(:,2)] = max(HH,[],2); 
for i = 1:size(H,1)
    if max(H_normailed(i,:))==1
        I_top2(i,:)=0;
    end
end
end

