function [result] = run_QuanTC(prodata,M,No_cluster,TC_cut)
result = {};
[H, C] = snmf(M,No_cluster);

% cpi value
cpi = zeros(size(prodata,1),1);
H_normailed =  H./sum(H,2);
for i = 1:length(cpi)
    for j = 1:No_cluster
        if H_normailed(i,j)>0
            cpi(i) = cpi(i) - H_normailed(i,j)*log(H_normailed(i,j));
        end
    end
    cpi(i) = cpi(i)/log(No_cluster);
end

% cluster-cluster
cluster_relation = H_normailed'*H_normailed;
cr = cluster_relation;

if size(H,2)>3
    cr(cr==min(min(cr(cr~=0))))=0;
    cr(cr==min(min(cr)))=0;%ignore smallest value
end


%plot potential transiton trajectory
G = digraph(cr,'OmitSelfLoops');
h = plot(G);
h.LineStyle = 'none';
h.NodeLabel =[];
h.MarkerSize=0.01;
hold on


%get the locaion of cluster centers & cells
cluster_location = zeros(No_cluster,2);
for i = 1:No_cluster
    cluster_location(i,1) = h.XData(i);
    cluster_location(i,2) = h.YData(i);
end
cell_location = zeros(size(prodata,1),2);
for i = 1:size(prodata,1)
    cell_location(i,1) = H_normailed(i,:)*cluster_location(:,1);
    cell_location(i,2) = H_normailed(i,:)*cluster_location(:,2);
end

close all
%get the weighted location of cells
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'MaxIterations',100,'OptimalityTolerance',1e-4,'Display','off');
[~,r] = pca(prodata);
ff = makefcn_weight(cluster_location,H_normailed,0.01/size(H_normailed,1),pdist(r)./max(pdist(r)));
x_weight = fminunc(ff,cell_location,options); %locations of cells



%threshold of CPI to find TC
C_TC = C;
C_TC(find(cpi>TC_cut)) = No_cluster + 1;












result{1} = H;
result{2} = C;
result{3} = cpi;
result{4} = cluster_location;
result{5} = cell_location;
result{6} = x_weight;
result{7} = C_TC;


















function fcn = makefcn_weight(a,p,lam,r)

fcn = @visu2_weight;
function [f,g] = visu2_weight(x)
% Calculate objective f
f = 0;
for jj = 1:size(p,2)
    f = f + sum(p(:,jj).*diag((x-ones(size(x,1),1).*a(jj,:))*(x-ones(size(x,1),1).*a(jj,:))'));
end
f = f - lam*sum(r*pdist(x,'squaredeuclidean')');

if nargout > 1 % gradient required
    g = zeros(size(x,1),2);

for jj = 1:size(p,2)
    g = g + 2*p(:,jj).*(x-ones(size(x,1),1).*a(jj,:));
end
g = g - 2*lam*(sum(squareform(r),2).*x-squareform(r)*x);

end
end
end

end

