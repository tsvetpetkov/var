%%%%% Restart program

close all;
clear;
clc;


%%%%% Load time-series data

% specify file name
Data.FileName = "sample_us_macro_data";
% import .csv file stored in the data folder and trim its labels
Data.Set = csvread(('../data/' + Data.FileName + '.csv'), 1, 0);

% specify the label of each data col
Data.Labels = {'Y', 'M', 'IP', 'UR', 'CRB', 'CPI', 'EBP', 'GS2'};
% specify the value order of each data col
Data.Values = [1, 2, 3, 4, 5, 6, 7, 8];

% map labels and values
Data.Map = containers.Map(Data.Labels, Data.Values);


%%%%% Specify VAR

VAR.tsFreq = 12;

VAR.monpolVar = {'GS1'};
VAR.monpolVarLbl = {'1-Year Rate (% pts)'};

VAR.macroVars = {'IP','CPI','CRB','EBP'};
VAR.macroVarsLbls = {'Industrial Production (%)','Consumer Prices (%)','Commodity Prices (%)','Excess Bond Premium (% pts)'};

VAR.laborVars = {'UR'};
VAR.laborVarsLbls = {'Unemployment Rate (% pts)'};

VAR.cholOrder = [1, 2, 3, 4, 5, 6];

VAR.startDate = [1992 1];
VAR.endDate = [2020 2];

VAR.lagLength = 12;
VAR.shockSize = 0.25;
VAR.irfHorizon = 61;


%%%%% Set up VAR

% retrieve sample start time
VAR.start = VAR.tsFreq*(VAR.startDate(1) - Data.Set(1, 1)) + VAR.startDate(2);
% retrieve sample end time
VAR.end = VAR.tsFreq*(VAR.endDate(1) - Data.Set(1, 1)) + VAR.endDate(2);

% construct VAR constant term
VAR.consTerm = ones((VAR.end - VAR.start - VAR.lagLength + 1), 1);

% construct set of variables
VAR.vars = [VAR.monpolVar, VAR.macroVars, VAR.laborVars];
% retrieve specified sample dataset
VAR.data = Data.Set(VAR.start:VAR.end, cell2mat(values(Data.Map, VAR.vars)));
% retrieve number of specified vars
VAR.numVars = numel(VAR.vars);


%%%%% Set up VAR outcome and explanatory variables

% construct Y values based on specified # of lags and Chol order
CholVAR.Y = VAR.data((VAR.lagLength + 1):end, VAR.cholOrder);

% construct X values based on specified # of vars and # of lags in Chol order
for ck = 1:VAR.lagLength
    % in blocks of (# of vars), add the k-th lag of Y
    CholVAR.lagData(:,((ck - 1) * VAR.numVars + 1):(ck * VAR.numVars)) = VAR.data((VAR.lagLength - ck + 1):(end - ck), VAR.cholOrder);
end
CholVAR.X = [CholVAR.lagData VAR.consTerm];


%%%%% Estimate Cholesky VAR

% retrieve reduced-form VAR beta coefficients
CholVAR.betas = CholVAR.X \ CholVAR.Y;
% retrieve reduced-form VAR residuals
CholVAR.residuals = CholVAR.Y - CholVAR.X * CholVAR.betas;

% retrieve VCV matrix of reduced-form residuals
CholVAR.residualsVCV = CholVAR.residuals' * CholVAR.residuals;

% retrieve Chol lower-triangular matrix form of the VCV matrix
CholVAR.VCVlower = chol(CholVAR.residualsVCV, 'lower');

% retrieve contemporaneous impact values
CholVAR.shockVals = CholVAR.VCVlower(:, find(VAR.cholOrder == 1));


%%%%% Create Cholesky IRFs

CholVAR.irfBuilder(VAR.lagLength + 1, :) = CholVAR.shockVals * (VAR.shockSize / CholVAR.shockVals(find(VAR.cholOrder == 1)));
for ch = 2:VAR.irfHorizon
    CholVAR.irfHelper = (CholVAR.irfBuilder((VAR.lagLength + ch - 1):-1:ch, :))';
    CholVAR.irfBuilder(VAR.lagLength + ch, :) = CholVAR.irfHelper(:)' * CholVAR.betas(1:(VAR.lagLength * VAR.numVars), :);
end
CholVAR.IRFs = CholVAR.irfBuilder((VAR.lagLength + 1):end, :);

%%%%% Plot IRFs

Fig.cols = 2;
Fig.numFigs = VAR.numVars;

Fig.cholIRFs = CholVAR.IRFs;
Fig.titles = VAR.varsLbls(VAR.cholOrder);

tiledlayout(ceil(Fig.numFigs / Fig.cols), Fig.cols, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for f = 1:Fig.numFigs
    nexttile
    set(gca, 'FontSize', 8.5)

    xticks(0:VAR.tsFreq:(VAR.irfHorizon - 1));
    grid minor

    yline(0, 'linewidth', 0.5, 'Color', "#A2142F", 'LineStyle', '-.');

    hold on
        plot(linspace(0, VAR.irfHorizon - 1, VAR.irfHorizon), Fig.cholIRFs(:, f), 'linewidth', 1, 'Color', "#4DBEEE", 'LineStyle', '-');
    hold off

    title(Fig.titles(f), 'FontSize', 11.5);
end
