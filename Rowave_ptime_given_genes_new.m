function [fitdata_ptime,new_order_gg] = Rowave_ptime_given_genes_new(data,allgenes,rw_gene,ptime,gene_cluster,folder)
% Rolling wave plot for given genes
% addpath('Data/TF');
% rw_gene = importdata('NB134.txt');
% folder: save figs;
% folder1: save .csv

[~,~,gene_idxv] = intersect(rw_gene,allgenes,'stable');

datav = data(gene_idxv,:);
g_name = allgenes(gene_idxv);



%% Fit polynomia, then order genes, plot along pseudotime
% Order genes according to their expression peaks

Current_genes1 = g_name;
data_ptime1 = datav(:,ptime);

fitdata = zeros(size(data_ptime1));

for i = 1:length(Current_genes1)
    x = (1:length(ptime))./length(ptime);
    p = polyfit(x,data_ptime1(i,:),3); %%%%%%%%%%%%%%%%
    fitdata(i,:) = polyvalue(x,p);
end


% gene ordered w.r.t peak
% [~,pk_idx1] = max(fitdata');
% [~,g_order1] = sort(pk_idx1);
 fitdata_ptime = fitdata;%(g_order1,:);
Glabs_ptime = g_name;%(g_order1); 


% RW_Genes = Glabs_ptime;
% T = table(RW_Genes);
% writetable(T,['Results/ds5scaled/Reduced10k/RWplot/RW_Gene_list_' num2str(topn) '.csv'],'Delimiter',',');  



%% visualization of fitted data along ptime
idata1 = fitdata_ptime;
kk = 2;
center = mean(idata1,kk);
scale = std(idata1, 0,kk);

tscale = scale;
%=Check for zeros and set them to 1 so not to scale them.
scale(tscale == 0) = 1;
%== Center and scale the data
idata1 = bsxfun(@minus, idata1, center);
sdata1 = bsxfun(@rdivide, idata1, scale);

% figure;
% colormap jet;
% 
% clims = [-3 3];
% imagesc(sdata1,clims);
% % title('RowlingWave Along Pseudotime');
% 
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% if length(g_name) <= 200
%     yticks(1:length(g_name));
%     % yticklabels(g_name);
%     yticklabels(Glabs_ptime);
% end
% set(gca,'fontsize',6)
% cb = colorbar;
% ax = gca;
% axpos = ax.Position;
% cpos = cb.Position;
% cpos(3) = 0.5*cpos(3);
% cb.Position = cpos;
% ax.Position = axpos;

%print([ folder '\RW_' rw_gene{1}],'-dpdf','-r300','-fillpage'); %'-dpdf',


%% plot gene orders based on the clustering of genes
Kms_idx = gene_cluster;
% new_order_g = gene_cluster; %[];
% 
% for i = 1:length(ave_corder)
%     new_order_g = [new_order_g, find(Kms_idx==ave_corder(i))];
% end

figure;
clims = [-3 3];

colormap jet;



% gene ordered w.r.t peak in each cluster
new_order_gg = [];
ave_corder = unique(Kms_idx,'stable');
for i = 1:length(ave_corder)
    a =find(Kms_idx==ave_corder(i));
    [~,pk_idx1] = max(fitdata(a,:)');
    [~,g_order1] = sort(pk_idx1);
    new_order_gg = [new_order_gg,a(g_order1)];
end
fitdata_ptime = fitdata(new_order_gg,:);




idata1 = fitdata_ptime;
kk = 2;
center = mean(idata1,kk);
scale = std(idata1, 0,kk);

tscale = scale;
%=Check for zeros and set them to 1 so not to scale them.
scale(tscale == 0) = 1;
%== Center and scale the data
idata1 = bsxfun(@minus, idata1, center);
sdata1 = bsxfun(@rdivide, idata1, scale);


imagesc(sdata1,clims);
% title('RowlingWave Along Pseudotime');

set(gca,'xtick',[]);
set(gca,'ytick',[]);
if length(g_name) <= 200
    yticks(1:length(g_name));
    % yticklabels(g_name);
    yticklabels(g_name(new_order_gg));
end
set(gca,'fontsize',6)


cb = colorbar;
% ax = gca;
% axpos = ax.Position;
% cpos = cb.Position;
% cpos(3) = 0.5*cpos(3);
% cb.Position = cpos;
% ax.Position = axpos;

%print([ folder '\RW_ordered' rw_gene{1}],'-dpdf','-r600','-fillpage'); %'-depsc -dpdf'



end

function polyy = polyvalue(x,p)
polyy = zeros(size(x));
for ii = 1:length(x)
    polyy(ii) = p(1).*x(ii)^3 + p(2).*x(ii)^2 + p(3).*x(ii) + p(4);
end
end



