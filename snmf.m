function [H, C] = snmf(M,No_cluster)
addpath('NMFLibrary');
addpath('NNDSVD');
addpath('symnmf2');
[wh, ~] = symm_anls(M, No_cluster,[]);
H = wh.W;
[~,C] = max(H,[],2);
end

