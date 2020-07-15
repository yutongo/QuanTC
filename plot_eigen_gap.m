function [eigenvalues] = plot_eigen_gap(W)
%input: similarity matrix W

WB = W;
n = size(W,1);
D = diag(WB*ones(n,1));
Prw = eye(size(W)) - D^(-1/2)*WB*D^(-1/2);
if n>=1000
    No_eigs = 100;
    all_eigs = real(eigs(Prw,No_eigs,'sm'));
else
    all_eigs = real(eig(Prw));
end
    
zz = sort(abs(real(all_eigs)));        
No_cluster1 = length(find(zz<=0.01));
NC1 = No_cluster1;



if length(zz)>=2
    gap = zz(2:end) - zz(1:end-1);
    [~,No_cluster1] = max(gap);
end
tol = 0.01;
No_cluster2 = length(find(zz<=tol));
No_cluster = No_cluster1;
display('Number of cluster based on zero eigenvalues & Largest gap ');
display([No_cluster2 No_cluster]);

eigenvalues = zz;

scatter(1:min([30 size(eigenvalues,1)]),eigenvalues(1:min([30 size(eigenvalues,1)])),20,'filled');
box on;
set(gca,'LineWidth',1.5);
xlabel('i');
ylabel('Eigenvalue of graph Laplacian \lambda_i');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);
%print(['Results\' folder '\EigenGap'],'-dpdf','-r300'); 
end

