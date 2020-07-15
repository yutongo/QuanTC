function [hbest,wbest,normbest] = nmf_knownH(a,k,h_known,lam,h0,w0)
% a is cell*gene, h should be kept as row-normalized


tolx = 1e-14;%statget(options,'TolX',defaultopt,'fast');
tolfun = 1e-3;%statget(options,'TolFun',defaultopt,'fast');
maxiter = 5000; %statget(options,'MaxIter',defaultopt,'fast');

% Special case, if K is full rank we know the answer
[n,m] = size(a);
if isempty(w0) && isempty(h0)
    if k==m
        w0 = a;
        h0 = eye(k);
    elseif k==n
        w0 = eye(k);
        h0 = a;
    end
end








sqrteps = sqrt(eps);
h=h0;w=w0;
for j=1:maxiter
        h00 = h;
        h = h00.*((a*w'+lam/2*h_known)./(h00*w*w'));
        % row normalize h each time
        h = h./sum(h,2);
        w = w.*((h00'*a)./(h00'*h00*w));
        
        
    % Get norm of difference and max change in factors
    %d = a - w*h;
    dnorm = (norm(a-h*w,'fro'))^2-lam*trace(h'*h_known);
    %sqrt(sum(sum(d.^2))/nm);
    %dnorm_first = (norm(a-h*w,'fro'))^2;
    dw = max(max(abs(w-w0) / (sqrteps+max(max(abs(w0))))));
    dh = max(max(abs(h-h0) / (sqrteps+max(max(abs(h0))))));
    delta = max(dw,dh);
    
    % Check for convergence
    if j>1
        if delta <= tolx
            break;
        elseif dnorm <= tolfun*max(1,dnorm0)
            dnorm = dnorm0;
            break;
        elseif j==maxiter
            break;
        elseif abs(dnorm-dnorm0) <=1e-8
            dnorm = dnorm0;
            break;
        elseif isnan(dnorm)
            dnorm = dnorm0;
            break;
        end
    end

    % Remember previous iteration results
    dnorm0 = dnorm;
    
end

hbest=h;
wbest=w;
normbest=dnorm;

end 
