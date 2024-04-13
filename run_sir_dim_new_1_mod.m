% Run dimensional SIR spreading model
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function [tt, nnS, nnI] = run_sir_dim_new_1_mod(N, k, D, alp, lam, gam, numS_0, numI_0, tt)

nnS = [];
nS = numS_0;
nnS = [nnS, nS];

nnI = [];
nI = numI_0;
nnI = [nnI, nI];

nR = N - (nS + nI);

aa = NaN;
rr = 0;

for t = tt(2:end)
    
    r = ((nI + nR - 1)/(alp/D) + (1/2)^D)^(1/D) - (1/2);%From central node to nodes on perimeter  
    if (r < 0)
        r = 0;%In case of numerical problems (r should always be nonnegative)
    end
    rr = [rr, r];
    
    if (r == 0)%Avoid division by 0
        a = 1;
    else
        a = min(1, nI/(alp*r^(D - 1)));%Estimate of the fraction of the boundary which has infected nodes
    end
    aa = [aa, a];
    
    b = 1 + (k - 1)*r^(D - 1)/(r^(D - 1) + (r + 1)^(D - 1) + (r + 2)^(D - 1));
    
    nu = lam*(1 - (1 - a)^b)*alp*(r + 1)^(D - 1);%Number of new infections
    
    nSNew = nS - nu;
    nINew = nI + nu - gam*nI;
    nRNew = nR + gam*nI;
    nS = nSNew;
    nnS = [nnS, nS];
    
    nI = nINew;
    nnI = [nnI, nI];
    
    nR = nRNew;
    
end
end