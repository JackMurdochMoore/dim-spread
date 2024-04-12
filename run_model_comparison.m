% 
% On empirical or synthetic network, compare:
%     Monte Carlo  
%     Dimensional spreading model
%     Homogeneous mean field
%     Heterogeneous mean field
%     PDMC
% 
% Moore et al. (2024), "Network spreading from network dimension"
% 
% Save results in .mat file
% 
% Inputs:
%     networkFlag - Specifier of network
%     lam -         Rate of infection
%     gam -         Rate of recovery
%     nTrials -     Number of Monte Carlo trials, each with an
%                   independently randomly chosen initially infected node
%
% Example usage:
% networkFlag = 0; lam = 0.1; gam = 0.05; nTrials = 100; run_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = -3; lam = 0.1; gam = 0.05; nTrials = 100; run_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = -4; lam = 0.1; gam = 0.05; nTrials = 100; run_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = 15; lam = 0.1; gam = 0.05; nTrials = 100; run_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = 23; lam = 0.1; gam = 0.05; nTrials = 100; run_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = 46; lam = 0.1; gam = 0.05; nTrials = 100; run_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = 56; lam = 0.1; gam = 0.05; nTrials = 100; run_model_comparison(networkFlag, lam, gam, nTrials); 
% 
% lam = 0.1; gam = 0.05; nTrials = 100; for networkFlag = [15, 23, 46, -4, -3, 0]; run_model_comparison(networkFlag, lam, gam, nTrials); close all; end
% lam = 0.05; gam = 0.2; nTrials = 100; for networkFlag = [15, 23, 46, -4, -3, 0]; run_model_comparison(networkFlag, lam, gam, nTrials); close all; end
%
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function run_model_comparison(networkFlag, lam, gam, nTrials)

saveResultsFolder = 'results-time-series-new';

N0 = 10000;
D = 2; p = 0; k = 2*D;
r = Inf;

% networkFlag = -4; r = Inf;% Inclusivity
% networkFlag = -3;% Small world, Manhattan, LCC
% networkFlag = 0;% Small world, Manhattan
% networkFlag = 1;% TV show% N = 3,892, sMax = 4, diam. = 20
% networkFlag = 2;% Power grid% N = 4,941, sMax = 12, diam. = 46
% networkFlag = 3;% Politician% N = 5,908, sMax = 4, diam. = 14
% networkFlag = 12;% Computer science PhD % N = 1,025, sMax = 10, diam. = 28
% networkFlag = 13;% Erdos collaboration network % N = 4,991
% networkFlag = 15; % arXiv General Relativity % N = 4,158
% networkFlag = 16; % School friendship network # 27 % N = 1,152
% networkFlag = 17; % School friendship network # 67 % N = 439
% networkFlag = 18; % Mouse visual cortex % N = 193
% networkFlag = 19; % Yeast protein interactions % N = 1,458, sMax = 5, diam. = 19
% networkFlag = 20; % WormNet DM-LC% N = 483, sMax = 5, diam. = 18
% networkFlag = 21; % WormNet DM-HT% N = 2,831, sMax = 7, diam. = 19
% networkFlag = 22; % Human disease % N = 516, sMax = 4, diam. = 15
% networkFlag = 23; % WormNet CE-LC% N = 993, sMax = 7, diam. = 22
% networkFlag = 24; % WormNet CE-HT% N = 2,194, sMax = 7, diam. = 20
% networkFlag = 25; % Links% N = 4253
% networkFlag = 26; % Roads, Minnesota% N = 2640
% networkFlag = 27; % Roads, Europe% N = 1039
% networkFlag = 28; % Circuit I% N = 122
% networkFlag = 29; % Circuit II% N = 252
% networkFlag = 30; % Circuit III% N = 512
% networkFlag = 41; % Escorts %N = 15810, sMax = 4, diam. = 17
% networkFlag = 46; % Face-to-face, Firth-2020-Fri (day leading to largest LCC)% N = 308, sMax = 4, diam. = 12
% networkFlag = 56; % Soaplands% N = 62,917, sMax = 6, diam. = 33

DLims = [1, Inf]; iiDimEstMeth = 1;
[~, ~, ~, codeCell, descriptionCell, ~, ~, ~, ~] = est_corr_dim_new_1(1, 2, DLims);
dimEstMethCode = codeCell{iiDimEstMeth};

ticInit = tic;

fontSize = 22; lineWidth = 3; markerSize = 8;
legendBoxAlpha = 0.25;

rngSeed = 0;
rng(rngSeed);

sigma = 0; kappa = 1; omega = 0; rewireFlag = 1; lowerAndUpperQuantile = [-eps, 1 + eps];

saveCell = {...
    'rngSeed',...
    'saveCell', 'nameStr', 'totalTime',...
    'sigma', 'kappa', 'omega', 'rewireFlag', 'lowerAndUpperQuantile',...
    'N0', 'D', 'k', 'p', 'lam', 'gam', ...
    'A', ...
    'DCorr', 'alp', 'sMax', 'ss1', 'ss', 'nn1', 'nn', 'N', 'N1',...
    'errSIRDimMean', 'errSIRMFMean', 'errSIRPDMCMean', 'errSIRHetMFMean',...
    'R0DEst', 'R0MF', 'R0PDMC', 'R0Het',...
    'tt', 'tt1',...
    'nnnI', 'nnnR', 'nnnS',...
    'nnIMean', 'nnI1ModelMF', 'nnI1ModelMFHet', 'nnI1ModelPDMC', 'nnI1ModelDim',...
    'nnSMean', 'nnS1ModelMF', 'nnS1ModelPDMC', 'nnS1ModelMFHet', 'nnS1ModelDim',...
    'nnRMean', 'nnR1ModelMF', 'nnR1ModelPDMC', 'nnR1ModelMFHet', 'nnR1ModelDim',...
    };

f = figure; colOrder = get(gca, 'ColorOrder'); close(f);
trueCol = colOrder(1, :);%Blue
standardCol = colOrder(4, :);%Purple
standard2Col = colOrder(2, :);%Red
% standard2Col = 0.5*[1, 1, 1];%Grey
standard3Col = colOrder(5, :);%Green
proposedCol = colOrder(3, :);%Yellow

if ismember(networkFlag, [0, -3])
    %Determine which side length L will provide L^D closest to N:
    L = N0^(1/D);
    LL = [floor(L), ceil(L)];
    NN = LL.^D;
    [~, indMin] = min(abs(NN - N0));
    N = NN(indMin);
else
    N = N0;
end

if (networkFlag == 0)
    A = small_world_manhattan(N, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile);
    nameStrNetwork = 'Small world';
    nameStr = ['pred-SIR-new-1_small-world-manhattan', '_N-', num2str(N0), '_D-', num2str(D), '_k-', num2str(k), '_p-', num2str(p), '_lambda-', num2str(lam), '_gamma-', num2str(gam), '_dim-lims-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', dimEstMethCode, '_', num2str(nTrials), 'trials'];
elseif (networkFlag == -3)
    A = small_world_manhattan_lcc(N, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile);
    nameStrNetwork = 'Small world';
    nameStr = ['pred-SIR-new-1_small-world-manhattan-lcc', '_N-', num2str(N0), '_D-', num2str(D), '_k-', num2str(k), '_p-', num2str(p), '_lambda-', num2str(lam), '_gamma-', num2str(gam), '_dim-lims-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', dimEstMethCode, '_', num2str(nTrials), 'trials'];
elseif (networkFlag == -4)
    N = N0;
    m = k/2; m0 = 2*m + 1;
    if (r < Inf)
        A = inclusivity(N, m, m0, r);
    else
        A = BA_mod_2(N, m, m0);
    end
    nameStrNetwork = 'Inclusivity';
    nameStr = ['pred-SIR-new-1_inclusivity', '_N-', num2str(N0), '_k-', num2str(k), '_r-', num2str(r), '_lambda-', num2str(lam), '_gamma-', num2str(gam), '_dim-lims-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', dimEstMethCode, '_', num2str(nTrials), 'trials'];
elseif (networkFlag > 0)
    [A, nameStrNetwork] = load_network(networkFlag);
    nameStr = ['pred-SIR-new-1_', nameStrNetwork, '_lambda-', num2str(lam), '_gamma-', num2str(gam), '_dim-lims-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', dimEstMethCode, '_', num2str(nTrials), 'trials'];
end

N = size(A, 1);
k = full(sum(A(:))/N);

G = graph(A);

[ss, nn] = count_distances(A);
[DVec, sMaxVec, ~, ~, ~, ~, ~, bVec, bMat] = est_corr_dim_new_1(ss, nn, DLims);
DCorr = DVec(iiDimEstMeth); sMax = sMaxVec(iiDimEstMeth); b = bVec(iiDimEstMeth); alp = b*sum(nn(1:sMax))/N;
ss1 = ss(ss <= (sMax - 1));
nn1 = nn(ss <= (sMax - 1));

N1 = 1 + sum(nn1)/N; if isempty(N1); N1 = 1; end

[compID, compSize] = conncomp(G);

NCList = NaN(nTrials, 1);

nnnI = NaN(nTrials, 1); nnnS = NaN(nTrials, 1); maxLength = 1;
nnnSWithINeighbour = NaN(nTrials, 1);

indNonS_0Cell = cell(1, nTrials);

numI_0 = 1;%Number initially infected
numR_0 = 0;%Number initially recovered
numS_0 = N - (numI_0 + numR_0);
numNonS_0 = numI_0 + numR_0;%Number initially infected or recovered
for iiTrial = 1:nTrials
    indNonS_0 = datasample(1:N, numNonS_0, 'Replace', false);%Initially infected or recovered
    indNonS_0Cell{iiTrial} = indNonS_0;
    
    NC = sum(compSize(unique(compID(indNonS_0(1:numI_0)))));%Total size of the components which can be infected
    NCList(iiTrial) = NC;
    
    [tt, nnS, nnI, nnSWithINeighbour, stateList] = run_sir_0_mod_2(A, lam, gam, indNonS_0, numI_0, numR_0);%Run ground truth SIR simulation
    
    nnR = N - nnS - nnI;
    
    thisLength = numel(nnS);
    if (thisLength > maxLength)
        diffLength = thisLength - maxLength;
        nnnS = [nnnS, repmat(nnnS(:, end), [1, diffLength])];
        nnnS(iiTrial, :) = nnS;
        nnnI = [nnnI, repmat(nnnI(:, end), [1, diffLength])];
        nnnI(iiTrial, :) = nnI;
        nnnSWithINeighbour = [nnnSWithINeighbour, repmat(nnnSWithINeighbour(:, end), [1, diffLength])];
        nnnSWithINeighbour(iiTrial, :) = nnSWithINeighbour;
        maxLength = thisLength;
    else
        diffLength = maxLength - thisLength;
        nnS = [nnS, ones(1, diffLength)*nnS(end)];
        nnnS(iiTrial, :) = nnS;
        nnI = [nnI, ones(1, diffLength)*nnI(end)];
        nnnI(iiTrial, :) = nnI;
    end
    
end

nnnR = N - nnnS - nnnI;

nnSMean = mean(nnnS, 1); nnIMean = mean(nnnI, 1);
nnRMean = N - nnSMean - nnIMean;

tt = 0:(numel(nnSMean) - 1);

% Dimensional spreading model:  
[~, nnS1ModelDim, nnI1ModelDim] = run_sir_dim_new_1_mod(N, k, DCorr, alp, lam, gam, numS_0, numI_0, tt);
nnR1ModelDim = N - (nnS1ModelDim + nnI1ModelDim);

% Mean field:
[~, nnS1ModelMF, nnI1ModelMF] = run_sir_hom_mean_field(A, lam, gam, numS_0, numI_0, tt);
nnR1ModelMF = N - (nnS1ModelMF + nnI1ModelMF);

% Heterogeneous mean field:
[~, nnS1ModelMFHet, nnI1ModelMFHet] = run_sir_het_mean_field(A, lam, gam, numS_0, numI_0, tt);
nnR1ModelMFHet = N - (nnS1ModelMFHet + nnI1ModelMFHet);

% PDMC:
[~, nnS1ModelPDMC, nnI1ModelPDMC] = run_sir_pdmc_1(A, lam, gam, numS_0, numI_0, tt);
nnR1ModelPDMC = N - (nnS1ModelPDMC + nnI1ModelPDMC);

ii_tF = find(nnI1ModelDim + nnR1ModelDim <= N1, 1, 'last');
region1 = 1:ii_tF;
tt1 = tt(region1); nnS1 = nnSMean(region1); nnI1 = nnIMean(region1); nnR1 = nnRMean(region1);

nnI1ModelDim = nnI1ModelDim(region1); nnI1ModelMF = nnI1ModelMF(region1); nnI1ModelPDMC = nnI1ModelPDMC(region1); nnI1ModelMFHet = nnI1ModelMFHet(region1);
nnS1ModelDim = nnS1ModelDim(region1); nnS1ModelMF = nnS1ModelMF(region1); nnS1ModelPDMC = nnS1ModelPDMC(region1); nnS1ModelMFHet = nnS1ModelMFHet(region1);
nnR1ModelDim = nnR1ModelDim(region1); nnR1ModelMF = nnR1ModelMF(region1); nnR1ModelPDMC = nnR1ModelPDMC(region1); nnR1ModelMFHet = nnR1ModelMFHet(region1);

errSIRDimMean = sqrt(mean((nnS1ModelDim - nnS1).^2 + (nnI1ModelDim - nnI1).^2 + (nnR1ModelDim - nnR1).^2));
errSIRMFMean = sqrt(mean((nnS1ModelMF - nnS1).^2 + (nnI1ModelMF - nnI1).^2 + (nnR1ModelMF - nnR1).^2));
errSIRPDMCMean = sqrt(mean((nnS1ModelPDMC - nnS1).^2 + (nnI1ModelPDMC - nnI1).^2 + (nnR1ModelPDMC - nnR1).^2));
errSIRHetMFMean = sqrt(mean((nnS1ModelMFHet - nnS1).^2 + (nnI1ModelMFHet - nnI1).^2 + (nnR1ModelMFHet - nnR1).^2));

legCell = {'True', 'Hom. MF', 'Het. MF', 'PDMC', 'Dim.'};

disp(['Mean RMS Error over t ∈ {0, 1, ..., ', num2str(max(tt1)), '} [Hom. MF, Het. MF, PDMC, Dim.]:']);
disp([errSIRMFMean, errSIRPDMCMean, errSIRHetMFMean, errSIRDimMean]);

%Plot number infected:

f = figure; hold on;
pl = plot(tt, nnIMean/N, 'o', tt1, nnI1ModelMF/N, '-.', tt1, nnI1ModelMFHet/N, ':', tt1, nnI1ModelPDMC/N, '--', tt1, nnI1ModelDim/N, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
pl(1).Color = trueCol; pl(1).MarkerFaceColor = trueCol;
pl(2).Color = standardCol; pl(2).MarkerFaceColor = standardCol;
pl(3).Color = standard2Col; pl(3).MarkerFaceColor = standard2Col;
pl(4).Color = standard3Col; pl(4).MarkerFaceColor = standard3Col;
pl(5).Color = proposedCol; pl(5).MarkerFaceColor = proposedCol;
box on; xlim([min(tt) - eps, min(max(tt), max(tt1)  + 1)]);
xlabel('Time, $t$', 'Interpreter', 'LaTeX', 'FontSize', fontSize); ylabel('Frac. infected, $I(t)/N$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);
yLim = ylim;
set(gca, 'YScale', 'Log');
ylim([nnIMean(2)/N, yLim(2)]);
leg = legend(legCell, 'Interpreter', 'LaTeX', 'FontSize', fontSize, 'Location', 'Best'); leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);

%Plot number susceptible:

f = figure; hold on;
pl = plot(tt, nnSMean/N, 'o', tt1, nnS1ModelMF/N, '-.', tt1, nnS1ModelPDMC/N, ':', tt1, nnS1ModelMFHet/N, '--', tt1, nnS1ModelDim/N, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
pl(1).Color = trueCol; pl(1).MarkerFaceColor = trueCol;
pl(2).Color = standardCol; pl(2).MarkerFaceColor = standardCol;
pl(3).Color = standard2Col; pl(3).MarkerFaceColor = standard2Col;
pl(4).Color = standard3Col; pl(4).MarkerFaceColor = standard3Col;
pl(5).Color = proposedCol; pl(5).MarkerFaceColor = proposedCol;
box on; xlim([min(tt) - eps, min(max(tt), max(tt1)  + 1)]);
xlabel('Time, $t$', 'Interpreter', 'LaTeX', 'FontSize', fontSize); ylabel('Frac. susceptible, $S(t)/N$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);
leg = legend(legCell, 'Interpreter', 'LaTeX', 'FontSize', fontSize, 'Location', 'Best'); leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);

%Plot number recovered:

f = figure; hold on;
pl = plot(tt, nnRMean/N, 'o', tt1, nnR1ModelMF/N, '-.', tt1, nnR1ModelPDMC/N, ':', tt1, nnR1ModelMFHet/N, '--', tt1, nnR1ModelDim/N, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
pl(1).Color = trueCol; pl(1).MarkerFaceColor = trueCol;
pl(2).Color = standardCol; pl(2).MarkerFaceColor = standardCol;
pl(3).Color = standard2Col; pl(3).MarkerFaceColor = standard2Col;
pl(4).Color = standard3Col; pl(4).MarkerFaceColor = standard3Col;
pl(5).Color = proposedCol; pl(5).MarkerFaceColor = proposedCol;
box on; xlim([min(tt) - eps, min(max(tt), max(tt1)  + 1)]);
xlabel('Time, $t$', 'Interpreter', 'LaTeX', 'FontSize', fontSize); ylabel('Frac. recovered, $R(t)/N$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);
yLim = ylim;
set(gca, 'YScale', 'Log');
ylim([nnRMean(2)/N, yLim(2)]);
leg = legend(legCell, 'Interpreter', 'LaTeX', 'FontSize', fontSize, 'Location', 'Best'); leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);

%Plot number affected:

f = figure; hold on;
pl = plot(tt, (nnIMean + nnRMean)/N, 'o', tt1, (nnI1ModelMF + nnR1ModelMF)/N, '-.', tt1, (nnI1ModelPDMC + nnR1ModelPDMC)/N, ':', tt1, (nnI1ModelMFHet + nnR1ModelMFHet)/N, '--', tt1, (nnI1ModelDim + nnR1ModelDim)/N, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
pl(1).Color = trueCol; pl(1).MarkerFaceColor = trueCol;
pl(2).Color = standardCol; pl(2).MarkerFaceColor = standardCol;
pl(3).Color = standard2Col; pl(3).MarkerFaceColor = standard2Col;
pl(4).Color = standard3Col; pl(4).MarkerFaceColor = standard3Col;
pl(5).Color = proposedCol; pl(5).MarkerFaceColor = proposedCol;
box on; xlim([min(tt) - eps, min(max(tt), max(tt1)  + 1)]);
xlabel('Time, $t$', 'Interpreter', 'LaTeX', 'FontSize', fontSize); ylabel('Frac. affected, $R(t)/N$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);
yLim = ylim;
set(gca, 'YScale', 'Log');
ylim([(nnIMean(1) + nnRMean(1))/N, yLim(2)]);
leg = legend(legCell, 'Interpreter', 'LaTeX', 'FontSize', fontSize, 'Location', 'Best'); leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);

R0DEst = calc_R0_D_new_1(k, DCorr, lam, gam);
R0MF = calc_R0_hom(A, lam, gam);
R0PDMC = calc_R0_PDMC(A, lam, gam);
R0Het = calc_R0_het(A, lam, gam);

disp('R0 [Hom. MF, Het. MF, PDMC, Dim.]:');
disp([R0MF, R0Het, R0PDMC, R0DEst]);

totalTime = toc(ticInit);
save([saveResultsFolder, '\', nameStr, '.mat'], saveCell{:});
disp(['That took ', num2str(totalTime), ' s.']);
end

% Run dimensional SIR spreading model
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2023
%
function [tt, nnS, nnI] = run_sir_dim_new_1_mod(N, k, D, alp, lam, gam, numS_0, numI_0, tt)

% tt = [];
% t = 0;
% tt = [tt, t];

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
    
    nu = lam*(1 - (1 - a)^b)*alp*(r + 1)^(D - 1);
    
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

% Run homogeneous mean field SIR spreading model
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function [tt, nnS, nnI] = run_sir_hom_mean_field(A, lam, gam, numS_0, numI_0, tt)

N = size(A, 1);

% tt = [];
% t = 0;
% tt = [tt, t];

nnS = [];
nS = numS_0;
nnS = [nnS, nS];

nnI = [];
nI = numI_0;
nnI = [nnI, nI];

nR = 0;

kList = sum(A, 1);
mean_k = mean(kList);

for t = tt(2:end)
    
    nSNew = nS - lam*nS*(1 - (1 - nI/N)^mean_k);
    nINew = nI + lam*nS*(1 - (1 - nI/N)^mean_k) - gam*nI;
    nRNew = nR + gam*nI;
    
    nS = nSNew;
    nnS = [nnS, nS];
    
    nI = nINew;
    nnI = [nnI, nI];
    
    nR = nRNew;
    
end
end

% Run heterogeneous mean field SIR spreading model
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2023
%
function [tt, nnS, nnI] = run_sir_het_mean_field(A, lam, gam, numS_0, numI_0, tt)

N = size(A, 1);

kk = sum(A, 1);
kMean = mean(kk);
kMin = min(kk);
kMax = max(kk);
kList = kMin:kMax;
kCountList = histcounts(kk, (kMin - 0.5):(kMax + 0.5));
kList = kList(kCountList > 0);
kCountList = kCountList(kCountList > 0);
kNum = numel(kCountList);

nnSList = NaN(kNum, 0);
nSList = (numS_0/N)*kCountList';
nnSList = [nnSList, nSList];

nnIList = NaN(kNum, 0);
nIList = (numI_0/N)*kCountList';
nnIList = [nnIList, nIList];

nRList = zeros(kNum, 1);

for t = tt(2:end)
    
    pIList = nIList./(kCountList');
    Th = (kList.*kCountList)*pIList/(N*kMean);%Probability that a randomly chosen link points to an infected node
    
    nSListNew = nSList - lam*(nSList.*(1 - (1 - Th).^(kList')));
    nIListNew = nIList + lam*(nSList.*(1 - (1 - Th).^(kList'))) - gam*nIList;
    nRListNew = nRList + gam*nIList;
    
    nSList = nSListNew;
    nnSList = [nnSList, nSList];
    
    nIList = nIListNew;
    nnIList = [nnIList, nIList];
    
    nRList = nRListNew;
    
end

nnS = sum(nnSList, 1);
nnI = sum(nnIList, 1);

end

% Run PDMC SIR spreading model
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function [tt, nnS, nnI] = run_sir_pdmc_1(A, lam, gam, nS0, nI0, tt)

N = size(A, 1);

SS0 = (nS0/N)*ones(N, 1);
II0 = (nI0/N)*ones(N, 1);
RR0 = 1 - (SS0 + II0);

SS = SS0;
nnS = [];
nS = sum(SS);
nnS = [nnS, nS];

II = II0;
nnI = [];
nI = sum(II);
nnI = [nnI, nI];

RR = RR0;
nnR = [];
nR = sum(RR);
nnR = [nnR, nR];

for t = tt(2:end)
    qq1 = exp(A*log(1 - II));%Fast calculation when no element of II is equal to 1
    
    SSNew = SS - lam*SS.*(1 - qq1);
    IINew = II + lam*SS.*(1 - qq1) - gam*II;
    RRNew = RR + gam*II;
    
    SS = SSNew;
    II = IINew;
    RR = RRNew;
    
    nS = sum(SS);
    nnS = [nnS, nS];
    
    nI = sum(II);
    nnI = [nnI, nI];
    
    nR = sum(RR);
    nnR = [nnR, nR];
    
end
end

% Calculate R_0 for SIR on network A using PDMC. 
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function [R0, r] = calc_R0_PDMC(A, lam, gam)

A = 0.5*(A + A');%Make sure the adjacency matrix is symmetric

degs = unique(full(sum(A)));
if (numel(degs) == 1)
    r = degs;
else
    try
        %s = eig(A); s = sort(s, 'descend', 'ComparisonMethod', 'real'); r = s(1);
        s = svd(full(A)); r = s(1);
    catch
        %[~, r, ~] = eigs(A, 1, 'largestabs'); r = abs(r);
        [~, r, ~] = svds(A, 1);
    end
end

if isnan(r)
    '';
end

R0 = lam*r/gam;

end

% 
% Calculate R_0 for SIR on network A using mean field for homogeneous networks. 
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function R0 = calc_R0_hom(A, lam, gam)

N = size(A, 1);

kList = sum(A, 1);
max_k = max(kList);
kk = 0:max_k;
pp_k = histcounts(kList, -0.5:(max_k + 0.5))/N;
kMean = kk*(pp_k');

R0 = lam/gam*kMean;

end

%
% Calculate R_0 for SIR on network A using mean field for heterogeneous networks. 
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
% 
function [R0, kMean, k2Mean] = calc_R0_het(A, lam, gam)

N = size(A, 1);

kList = sum(A, 1);
max_k = max(kList);
kk = 0:max_k;
pp_k = histcounts(kList, -0.5:(max_k + 0.5))/N;
kMean = kk*(pp_k');
k2Mean = kk.^2*(pp_k');

R0 = lam/gam*k2Mean/kMean;

end

%
% Calculate R_0 for SIR on network with correlation dimension D and mean
% degree k. Start from r = 1, I = 0. 
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function R0 = calc_R0_D_new_1(k, D, lam, gam)
r = 1;
b = 1 + (k - 1)/(1 + (1 + r^-1)^(D - 1) + (1 + 2*r^-1)^(D - 1));
R0 = lam./gam.*b.*(1 + r^-1)^(D - 1);
end