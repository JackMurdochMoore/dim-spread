% On empirical or synthetic network, compare:
%     Monte Carlo  
%     Dimensional spreading model
%     Reduced effective degree model (slight modification of model outlined
%       in Sec. 2.2.3 of "Epidemic Spread in Networks: Existing Methods and
%       Current Challenges" by  Miller and Kiss (2014))
%     Homogeneous pair approximation ("The effects of local spatial structure on
%       epidemiological invasions" by Keeling (2011)) 
%     Pair-based model (“Deterministic epidemic models on contact networks:
%       correlations and unbiological terms” by Sharkey (2010)) 
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
% networkFlag = 0; lam = 0.1; gam = 0.05; nTrials = 100; run_alt_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = -3; lam = 0.1; gam = 0.05; nTrials = 100; run_alt_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = -4; lam = 0.1; gam = 0.05; nTrials = 100; run_alt_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = 2; lam = 0.1; gam = 0.05; nTrials = 100; run_alt_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = 5; lam = 0.1; gam = 0.05; nTrials = 100; run_alt_model_comparison(networkFlag, lam, gam, nTrials); 
% networkFlag = 7; lam = 0.1; gam = 0.05; nTrials = 100; run_alt_model_comparison(networkFlag, lam, gam, nTrials); 
% 
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function run_alt_model_comparison(networkFlag, lam, gam, nTrials)

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
    'errSIRDimMean', 'errSIRPAMean', 'errSIRPairBasedMean', 'errSIRREDMMean',...
    'tt', 'tt1',...
    'nnnI', 'nnnR', 'nnnS',...
    'nnIMean', 'nnI1ModelPA', 'nnI1ModelREDM', 'nnI1ModelPairBased', 'nnI1ModelDim',...
    'nnSMean', 'nnS1ModelPA', 'nnS1ModelPairBased', 'nnS1ModelREDM', 'nnS1ModelDim',...
    'nnRMean', 'nnR1ModelPA', 'nnR1ModelPairBased', 'nnR1ModelREDM', 'nnR1ModelDim',...
    };

f = figure; colOrder = get(gca, 'ColorOrder');
trueCol = colOrder(1, :);%Blue
standardCol = colOrder(7, :);%Burgundy/Maroon/Dark red
standard2Col = 0.5*[1, 1, 1];%Grey
standard3Col = colOrder(6, :);%Light blue
proposedCol = colOrder(3, :);%Yellow
close(f);

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
    nameStr = ['pred-SIR-new-1-alt_small-world-manhattan', '_N-', num2str(N0), '_D-', num2str(D), '_k-', num2str(k), '_p-', num2str(p), '_lambda-', num2str(lam), '_gamma-', num2str(gam), '_dim-lims-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', dimEstMethCode, '_', num2str(nTrials), 'trials'];
elseif (networkFlag == -3)
    A = small_world_manhattan_lcc(N, k, D, p, sigma, kappa, omega, rewireFlag, lowerAndUpperQuantile);
    nameStrNetwork = 'Small world';
    nameStr = ['pred-SIR-new-1-alt_small-world-manhattan-lcc', '_N-', num2str(N0), '_D-', num2str(D), '_k-', num2str(k), '_p-', num2str(p), '_lambda-', num2str(lam), '_gamma-', num2str(gam), '_dim-lims-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', dimEstMethCode, '_', num2str(nTrials), 'trials'];
elseif (networkFlag == -4)
    N = N0;
    m = k/2; m0 = 2*m + 1;
    if (r < Inf)
        A = inclusivity(N, m, m0, r);
    else
        A = BA_mod_2(N, m, m0);
    end
    nameStrNetwork = 'Inclusivity';
    nameStr = ['pred-SIR-new-1-alt_inclusivity', '_N-', num2str(N0), '_k-', num2str(k), '_r-', num2str(r), '_lambda-', num2str(lam), '_gamma-', num2str(gam), '_dim-lims-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', dimEstMethCode, '_', num2str(nTrials), 'trials'];
elseif (networkFlag > 0)
    [A, nameStrNetwork] = load_network(networkFlag);
    nameStr = ['pred-SIR-new-1-alt_', nameStrNetwork, '_lambda-', num2str(lam), '_gamma-', num2str(gam), '_dim-lims-', num2str(DLims(1)), '-', num2str(DLims(2)), '_', dimEstMethCode, '_', num2str(nTrials), 'trials'];
end

N = size(A, 1);
k = full(sum(A(:))/N);

G = graph(A);
[ss, nn] = count_distances(A);
[DVec, sMaxVec, ~, ~, ~, ~, ~, bVec, bMat] = est_corr_dim_new_1(ss, nn, DLims);
DCorr = DVec(iiDimEstMeth); sMax = sMaxVec(iiDimEstMeth); b = bVec(iiDimEstMeth); alp = b*sum(nn(1:sMax))/N;
D = DCorr;
ss1 = ss(ss <= (sMax - 1));
nn1 = nn(ss <= (sMax - 1));

N1 = 1 + sum(nn1)/N; if isempty(N1); N1 = 1; end

[compID, compSize] = conncomp(G);

NCList = NaN(nTrials, 1);

nnnI = NaN(nTrials, 1); nnnS = NaN(nTrials, 1); maxLength = 1;

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
    
    [~, nnS, nnI, ~, ~] = run_sir_0_mod_2(A, lam, gam, indNonS_0, numI_0, numR_0);%Run ground truth SIR simulation
    
    nnR = N - nnS - nnI;
    
    thisLength = numel(nnS);
    if (thisLength > maxLength)
        diffLength = thisLength - maxLength;
        nnnS = [nnnS, repmat(nnnS(:, end), [1, diffLength])];
        nnnS(iiTrial, :) = nnS;
        nnnI = [nnnI, repmat(nnnI(:, end), [1, diffLength])];
        nnnI(iiTrial, :) = nnI;
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

% Pair approximation:  
[~, nnS1ModelPA, nnI1ModelPA] = run_sir_pa(A, lam, gam, numS_0, numI_0, tt);
nnR1ModelPA = N - (nnS1ModelPA + nnI1ModelPA);

% Reduced effective degree model:  
[~, nnS1ModelREDM, nnI1ModelREDM] = run_sir_redm(A, lam, gam, tt);
nnR1ModelREDM = N - (nnS1ModelREDM + nnI1ModelREDM);

% Pair-based:  
[~, nnS1ModelPairBased, nnI1ModelPairBased] = run_sir_pair_based_model(A, lam, gam, tt);
nnR1ModelPairBased = N - (nnS1ModelPairBased + nnI1ModelPairBased);

ii_tF = find(nnI1ModelDim + nnR1ModelDim <= N1, 1, 'last');
region1 = 1:ii_tF;
tt1 = tt(region1); nnS1 = nnSMean(region1); nnI1 = nnIMean(region1); nnR1 = nnRMean(region1);

nnI1ModelDim = nnI1ModelDim(region1); nnI1ModelPA = nnI1ModelPA(region1); nnI1ModelPairBased = nnI1ModelPairBased(region1); nnI1ModelREDM = nnI1ModelREDM(region1);
nnS1ModelDim = nnS1ModelDim(region1); nnS1ModelPA = nnS1ModelPA(region1); nnS1ModelPairBased = nnS1ModelPairBased(region1); nnS1ModelREDM = nnS1ModelREDM(region1);
nnR1ModelDim = nnR1ModelDim(region1); nnR1ModelPA = nnR1ModelPA(region1); nnR1ModelPairBased = nnR1ModelPairBased(region1); nnR1ModelREDM = nnR1ModelREDM(region1);

errSIRDimMean = mean(sqrt((nnS1ModelDim - nnS1).^2 + (nnI1ModelDim - nnI1).^2 + (nnR1ModelDim - nnR1).^2));
errSIRPAMean = mean(sqrt((nnS1ModelPA - nnS1).^2 + (nnI1ModelPA - nnI1).^2 + (nnR1ModelPA - nnR1).^2));
errSIRPairBasedMean = mean(sqrt((nnS1ModelPairBased' - nnS1).^2 + (nnI1ModelPairBased' - nnI1).^2 + (nnR1ModelPairBased' - nnR1).^2));
errSIRREDMMean = mean(sqrt((nnS1ModelREDM' - nnS1).^2 + (nnI1ModelREDM' - nnI1).^2 + (nnR1ModelREDM' - nnR1).^2));

legCell = {'True', 'Hom. pair', 'Red. eff. deg.', 'Pair-based', 'Dim.'};

disp(['Mean error in (S(t)/N, I(t)/N, R(t)/N) over t ∈ {0, 1, ..., ', num2str(max(tt1)), '} [Hom. pair, Red. eff. deg., Pair-based, Dim.]:']);
disp([errSIRPAMean/N, errSIRREDMMean/N, errSIRPairBasedMean/N, errSIRDimMean/N]);

%Plot number infected:

f = figure; hold on;
pl = plot(tt, nnIMean/N, 'o', tt1, nnI1ModelPA/N, '-.', tt1, nnI1ModelREDM/N, ':', tt1, nnI1ModelPairBased/N, '--', tt1, nnI1ModelDim/N, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
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
pl = plot(tt, nnSMean/N, 'o', tt1, nnS1ModelPA/N, '-.', tt1, nnS1ModelPairBased/N, ':', tt1, nnS1ModelREDM/N, '--', tt1, nnS1ModelDim/N, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
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
pl = plot(tt, nnRMean/N, 'o', tt1, nnR1ModelPA/N, '-.', tt1, nnR1ModelPairBased/N, ':', tt1, nnR1ModelREDM/N, '--', tt1, nnR1ModelDim/N, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
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
pl = plot(tt, (nnIMean + nnRMean)/N, 'o', tt1, (nnI1ModelPA + nnR1ModelPA)/N, '-.', tt1, (nnI1ModelPairBased + nnR1ModelPairBased)/N, ':', tt1, (nnI1ModelREDM + nnR1ModelREDM)/N, '--', tt1, (nnI1ModelDim + nnR1ModelDim)/N, '-', 'LineWidth', lineWidth, 'MarkerSize', markerSize);
pl(1).Color = trueCol; pl(1).MarkerFaceColor = trueCol;
pl(2).Color = standardCol; pl(2).MarkerFaceColor = standardCol;
pl(3).Color = standard2Col; pl(3).MarkerFaceColor = standard2Col;
pl(4).Color = standard3Col; pl(4).MarkerFaceColor = standard3Col;
pl(5).Color = proposedCol; pl(5).MarkerFaceColor = proposedCol;
box on; xlim([min(tt) - eps, min(max(tt), max(tt1)  + 1)]);
xlabel('Time, $t$', 'Interpreter', 'LaTeX', 'FontSize', fontSize); ylabel('Frac. affected, $(I + R)/N$', 'Interpreter', 'LaTeX', 'FontSize', fontSize);
set(gca, 'TickLabelInterpreter', 'LaTeX', 'FontSize', fontSize);
yLim = ylim;
set(gca, 'YScale', 'Log');
ylim([(nnIMean(1) + nnRMean(1))/N, yLim(2)]);
leg = legend(legCell, 'Interpreter', 'LaTeX', 'FontSize', fontSize, 'Location', 'Best'); leg.BoxFace.ColorType='truecoloralpha'; leg.BoxFace.ColorData=uint8(255*[1; 1; 1; legendBoxAlpha]);

totalTime = toc(ticInit);
save([saveResultsFolder, '\', nameStr, '.mat'], saveCell{:});
disp(['That took ', num2str(totalTime), ' s.']);

end

% Implementation of SIR model in "The effects of local spatial structure on
% epidemiological invasions" by Keeling (2011)
% 
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function [tt, nnS, nnI] = run_sir_pa(A, lam, gam, numS_0, numI_0, tt)

N = size(A, 1); 
kMean = full(sum(A(:)))/N;%Mean degree
numTriangles = trace(A^3)/6;
B = A^2;
numConnectedTriples = (sum(B(:)) - trace(B))/2;
phi = numTriangles/numConnectedTriples;%Clustering coefficient

minInitVal = 10^-12;

%Choose initial conditions:
%First order (make sure all values are positive):
SIR0 = [numS_0, numI_0, N - (numS_0 + numI_0)];
initZeroInds = (SIR0 <= 0);
numInitZeroInds = sum(initZeroInds); numInitNonZeroInds = 3 - numInitZeroInds;
SIR0(initZeroInds) = minInitVal; SIR0(~initZeroInds) = SIR0(~initZeroInds) + minInitVal*numInitZeroInds/numInitNonZeroInds;
SIR0 = N*SIR0/sum(SIR0);%Normalise
S0 = SIR0(1);
I0 = SIR0(2);
R0 = SIR0(3);
%Second order:
SS0 = N*kMean*(S0/N)*(S0/N);
SI0 = N*kMean*(S0/N)*(I0/N);
SR0 = N*kMean*(S0/N)*(R0/N);
II0 = N*kMean*(I0/N)*(I0/N);
IR0 = N*kMean*(I0/N)*(R0/N);
RR0 = N*kMean*(R0/N)*(R0/N);
%Vector of initial conditions:
y0 = ...
    [S0; I0; R0;
    SS0; SI0; SR0;
    II0; IR0;
    RR0];

% %Check: should have the following (unless the tolerance tol is too small):
% tol = 10^-10;
% abs(sum(y0(1:3)) - N) < tol;
% abs(sum(y0(4:9)) - n*N/2) < tol;
% 
% %Check: should have the following (unless the tolerance tol is too small):
% tol = 10^-10;
% assert(abs((I0 + R0 + S0) - N) < tol);
% assert(abs((SS0 + 2*SI0 + 2*SR0 + II0 + 2*IR0 + RR0) - n*N) < tol);

oDEFun = @(t, y) deriv_SIR_moment_closure(y, N, kMean, phi, lam, gam);
[tt, y] = ode45(oDEFun, tt, y0);

%Extracting results:
%First order:
S = y(:, 1); nnS = S';
I = y(:, 2); nnI = I';
R = y(:, 3);
%Second order:
SS = y(:, 4);
SI = y(:, 5);
SR = y(:, 6);
II = y(:, 7);
IR = y(:, 8);
RR = y(:, 9);

end

function [dydt] = deriv_SIR_moment_closure(y, N, n, phi, lambda, gamma)

S = y(1);
I = y(2);
R = y(3);%Note to self: sum(y(1:3)) == N;

first_order = [S; I; R];

SS = y(4);
SI = y(5);
SR = y(6);
II = y(7);
IR = y(8);
RR = y(9);%Note to self: sum(y(4:9)) == N*n;

second_order = ...
    [SS, SI, SR;
    SI, II, IR;
    SR, IR, RR];

iiS = 1; iiI = 2; iiR = 3;

SSI = calc_third_order(iiS, iiS, iiI, first_order, second_order, N, n, phi);
ISI = calc_third_order(iiI, iiS, iiI, first_order, second_order, N, n, phi);
ISR = calc_third_order(iiI, iiS, iiR, first_order, second_order, N, n, phi);

dS = -lambda*SI;
dI = lambda*SI - gamma*I;
dR = gamma*I;

dSS = -2*lambda*SSI;
dSI = lambda*(SSI - ISI - SI) - gamma*SI;
dSR = -lambda*ISR + gamma*SI;
dII = 2*lambda*(ISI + SI) - 2*gamma*II;
dIR = lambda*ISR + gamma*(II - IR);
dRR = 2*gamma*IR;

dydt = ...
    [dS; dI; dR;
    dSS; dSI; dSR;
    dII; dIR;
    dRR];%Collect derivatives into a single column vector.
end

function XYZ = calc_third_order(iiX, iiY, iiZ, first_order, second_order, N, n, phi)
X = first_order(iiX);
Y = first_order(iiY);
Z = first_order(iiZ);
XY = second_order(iiX, iiY);
YZ = second_order(iiY, iiZ);
XZ = second_order(iiX, iiZ);
if (X < 10^(-6)) || (Y < 10^(-6)) || (Z < 10^(-6))%Avoid numerical issues
    XYZ = 0;
else
    XYZ = ((n - 1)/n)*(XY*YZ/Y)*((1 - phi) + (phi*N/n)*(XZ/(X*Z)));
end
end

% Reduced effective degree model
% Miller and Kiss (2014) Epidemic Spread in Networks: Existing Methods and Current Challenges
% Slightly modified for our case for which 1 or more infected neighbours
% provides the same chance of infection 
% 
%
% Moore et al. (2024), "Network spreading from network dimension"
%
% Jack Moore, 2024
%
function [tt, nnS, nnI] = run_sir_redm(A, lam, gam, tt)

N = size(A, 1);

kk = sum(A, 1);
kMean = mean(kk);
kEffMin = 0;
kEffMax = max(kk);
kCountList = histcounts(kk, (kEffMin - 0.5):(kEffMax + 0.5));
kNum = numel(kCountList);

n0 = kMean; %Number of effective partnerships between susceptible and infected individuals
xx0 = (N - 1)/N*kCountList';%Entry j + 1 is number of susceptible nodes with j infected partners
R0 = 0;

y0 = [n0; xx0; R0];

oDEFun = @(t, y) deriv_SIR_redm(y, N, lam, gam);
[tt, y] = ode45(oDEFun, tt, y0);

nnR = y(:, end);
xxList = y(:, 2:(end - 1));
nnS = sum(xxList, 2);
nnI = N - nnS - nnR;

end

function [dydt] = deriv_SIR_redm(y, N, lam, gam)

n = y(1);
xx = y(2:(end - 1));
R = y(end);
S = sum(xx);
I = N - S - R;
kEffList = 0:(numel(xx) - 1);%Effective degrees
IMean = n/(kEffList*xx);
dn = lam*(1 - 2*IMean)*((1 - (1 - IMean).^kEffList).*(kEffList - 1))*xx -(lam + gam)*n;
dxx = gam*IMean*((kEffList' + 1).*[xx(2:end); 0] - kEffList'.*xx) - lam*xx.*(1 - (1 - IMean).^kEffList');
dR = gam*I;

dydt = [dn; dxx; dR];%Collect derivatives into a single column vector.
end