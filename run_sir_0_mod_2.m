% Run ground truth discrete-state SIR simulations
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2023
%
function [tt, nnS, nnI, nnSWithINeighbour, stateList] = run_sir_0_mod_2(A, lam, gam, indNonS_0, numI_0, numR_0)
N = size(A, 1);
states = zeros(N, 1);
states(indNonS_0(1:numI_0)) = 1;
states(indNonS_0(numI_0 + (1:numR_0))) = 2;

tt = [];
t = 0;
tt = [tt, t];

nnS = [];
nS = sum(states == 0);
nnS = [nnS, nS];

nnI = [];
nI = sum(states == 1);
nnI = [nnI, nI];

nR = N - (nS + nI);

nnSWithINeighbour = [];

stateList = NaN(N, 0);
stateList = [stateList, states];

G = graph(A); connComp = conncomp(G);
reachableComp = connComp(indNonS_0(1:numI_0));
NReachable = sum(ismember(connComp, reachableComp));

while (nI > 0) && (nR < NReachable)
    if (gam == 0)
        if (lam == 0) || (nI == NReachable)
            break;
        end
    end
    t = t + 1;
    oldStates = states;
    newStates = oldStates;
    %0 - susceptible, 1 - infected, 2 - recovered
    ILogicals = (oldStates == 1);
    hasINeighbourLogical = any(A(:, ILogicals), 2);
    
    SNodes = find(oldStates == 0); numS = numel(SNodes);
    SWithINeighbourLogical = hasINeighbourLogical(SNodes);
    SThreshProbs = lam.*SWithINeighbourLogical;
    SUnifRands = rand(numS, 1);
    newStates(SNodes(SUnifRands < SThreshProbs)) = 1;%S -> I
    
    INodes = find(ILogicals);
    IToRIndices = INodes(rand(numel(INodes), 1) < gam);
    newStates(IToRIndices) = 2;%I -> R
    states = newStates;
    
    tt = [tt, t];
    
    nS = sum(states == 0);
    nnS = [nnS, nS];
    
    nI = sum(states == 1);
    nnI = [nnI, nI];
    
    nR = N - (nS + nI);
    
    stateList = [stateList, states];
    
    nSWithINeighbour = sum(SWithINeighbourLogical);
    nnSWithINeighbour = [nnSWithINeighbour, nSWithINeighbour];
end
nnSWithINeighbour = [nnSWithINeighbour, 0];
end