%% Data Analysis Code for Torres et al. 2020
% "Long-term storage and age-biased export of fluvial organic carbon:
%  field evidence from West Iceland"

% This code will produce a figure similar to Figure 12 in the manuscript,
% but for only one randomly chosen value of tCut between 50 and 650 years

% THIS CODE USES FILES FROM:
% Lougheed, Bryan, and Stephen Obrochta. 
% "MatCal: Open source Bayesian 14 C age calibration in MatLab." 
% Journal of Open Research Software 4.1 (2016).
%
% Specifically, I've taken the intCal13data file and included it in this 
% repository to be used alongside this file

% This code also requires the temperlevyRAND.m function, which was based on
% a similar script written by Dr. Vamsi Ganti

% on a 2013 3.1 GHz Intel Core i7 (16 GB 1600 MHz DDR3), 
% this code runs in about 3 seconds

%% Define Model Parameters
% transit time distribution parameters
gamT = 1.2; % Pareto distribution parameter gamma
aT = 0.8; % Pareto distribution parameter alpha
bT = 120; % Pareto distribution parameter beta
p = 1E6; % number of random draws from age distribution
nt = 1; % max number of transport events 
tCut = random('uniform',50,650); %bend growth timescale (years) 
threshAge = 1.3E4; %Maximum allowable age
E_v = 2E-4; % aggradation rate, channel depths per year

%% Generate waiting time distributions
%dimensionless waiting distribution for exported sediments (ndWTDexprt)
longndWTDexprt = temperlevyRAND([p,nt], gamT, aT, bT); %tempered pareto variables
shortWTDexprt = exprnd(4,[p,nt]); % exponential variables (for t < tmin)
shortWTDexprt = shortWTDexprt(shortWTDexprt<min(longndWTDexprt)); %remove exp variables > tmin
ndWTDexprt = [longndWTDexprt;shortWTDexprt]; %combine pareto and exponential
WTDexprt = ndWTDexprt.*tCut; %dimensionalize

MWT = mean(WTDexprt); %mean waiting time
[Fwtd,WT] = ecdf(WTDexprt,'function','survivor'); %temperLevy S(x) generation
pWTstor = Fwtd./MWT;%PDF of waitings times in storage
[~,unind] = unique(WT); %unique values (for fit)
pWTstorfit = fit(WT(unind),pWTstor(unind),'linearinterp'); %linear interpolation
CDFwtd = integrate(pWTstorfit,linspace(0,max(WT),1E6),0); %integrate PDF for CDF
%Calculate inverse of CDF
invCDFwtd = fit(CDFwtd',linspace(0,max(WT),1E6)','linearinterp'); 
WTDstor = feval(invCDFwtd,rand(p,1)); %random values using inverse transform sampling

clearvars -except p tCut threshAge E_v WTDexprt WTDstor

%% Account for aggradation
WTDexprtAG = WTDexprt; %copy exported waiting times
% for those less than aggradation timescale, calculate age as time since
% initial deposition and a linear sedimentation rate
WTDexprtAG(WTDexprtAG<(1/E_v)) = WTDexprtAG(WTDexprtAG<(1/E_v)) + (-0.5*E_v.*WTDexprtAG(WTDexprtAG<(1/E_v)).^2);
% for those greater than aggradation timescale, mean age is fixed
WTDexprtAG(WTDexprtAG>=(1/E_v)) = 0.5*(1/E_v);
% Residence Time Distribution
WTDstorAG = WTDstor; %copy stored waiting times
% for those less than aggradation timescale, calculate age as time since
% initial deposition and a linear sedimentation rate
WTDstorAG(WTDstorAG<(1/E_v)) = WTDstorAG(WTDstorAG<(1/E_v)) + (-0.5*E_v.*WTDstorAG(WTDstorAG<(1/E_v)).^2);
% for those greater than aggradation timescale, mean age is fixed
WTDstorAG(WTDstorAG>=(1/E_v)) = 0.5*(1/E_v);

%% Calculate radiocarbon
% Convert age into Fm taking into account atmospheric variations
% Bomb spike dealt with separately
load('intCal13data.mat') % Load MatCal data 
%fit curve of time vs. atmospheric 14C (as Fm)
fmCalcurve = fit(hicurvecal',hicurvef14','linearinterp');
% NOTE: Matcal data in years BP (i.e., before 1950).
% Field data measured in 2017 

% Apply Calibration
youngerBPindWTDexportAG = find(WTDexprtAG<=67); %find ages younger than bomb spike
olderBPindWTDexportAG = find(WTDexprtAG>67); %find ages older than bomb spike
fmWTDexportAG = zeros(length(WTDexprtAG),1); %pre-allocation
%evaluate Fm for pre bomb spike
fmWTDexportAG(olderBPindWTDexportAG) = feval(fmCalcurve,WTDexprtAG(olderBPindWTDexportAG)-67);
fmWTDexportAG(youngerBPindWTDexportAG) = 1.5; %set to same high Fm for post bomb spike

youngerBPindWTDstorAG = find(WTDstorAG<=67);%find ages younger than bomb spike
olderBPindWTDstorAG = find(WTDstorAG>67);%find ages older than bomb spike
fmWTDstorAG = zeros(length(WTDstorAG),1);%pre-allocation
%evaluate Fm for pre bomb spike
fmWTDstorAG(olderBPindWTDstorAG) = feval(fmCalcurve,WTDstorAG(olderBPindWTDstorAG)-67);
fmWTDstorAG(youngerBPindWTDstorAG) = 1.5;%set to same high Fm for post bomb spike

% Save mean values
bulkStoredC = mean(fmWTDstorAG);
bulkExportedC = mean(fmWTDexportAG);

%% Apply threshold and re-calculate radiocarbon
%copy exported waiting times, but remove values greater than threshold
WTDexprtAGt = WTDexprt(WTDexprt<threshAge); 
WTDexprtAGt(WTDexprtAGt<(1/E_v)) = WTDexprtAGt(WTDexprtAGt<(1/E_v)) + (-0.5*E_v.*WTDexprtAGt(WTDexprtAGt<(1/E_v)).^2);
WTDexprtAGt(WTDexprtAGt>=(1/E_v)) = 0.5*(1/E_v);
WTDstorAGt = WTDstor(WTDstor<threshAge); %copy subset of waiting times
WTDstorAGt(WTDstorAGt<(1/E_v)) = WTDstorAGt(WTDstorAGt<(1/E_v)) + (-0.5*E_v.*WTDstorAGt(WTDstorAGt<(1/E_v)).^2);
WTDstorAGt(WTDstorAGt>=(1/E_v)) = 0.5*(1/E_v);

% Apply Calibration
youngerBPindWTDexportAGt = find(WTDexprtAGt<=67); %find ages younger than bomb spike
olderBPindWTDexportAGt = find(WTDexprtAGt>67); %find ages older than bomb spike
fmWTDexportAGt = zeros(length(WTDexprtAGt),1); %pre-allocation
fmWTDexportAGt(olderBPindWTDexportAGt) = feval(fmCalcurve,WTDexprtAGt(olderBPindWTDexportAGt)-67);
fmWTDexportAGt(youngerBPindWTDexportAGt) = 1.5; %set to same high Fm for post bomb spike

youngerBPindWTDstorAGt = find(WTDstorAGt<=67);%find ages younger than bomb spike
olderBPindWTDstorAGt = find(WTDstorAGt>67);%find ages older than bomb spike
fmWTDstorAGt = zeros(length(WTDstorAGt),1);%pre-allocation
fmWTDstorAGt(olderBPindWTDstorAGt) = feval(fmCalcurve,WTDstorAGt(olderBPindWTDstorAGt)-67);
fmWTDstorAGt(youngerBPindWTDstorAGt) = 1.5;%set to same high Fm for post bomb spike

% Save mean values
bulkStoredCt = mean(fmWTDstorAGt);
bulkExportedCt = mean(fmWTDexportAGt);

%% Moraine OC mixing
fgla = 0.25; %Fraction Moraine Carbon
fflu = 1-fgla; %Fraction Fluvial Carbon
bulkStoredCm = (fflu.*bulkStoredC)+(fgla.*0.27);
bulkStoredCtm = (fflu.*bulkStoredCt)+(fgla.*0.27);
bulkExportedCtm = (fflu.*bulkExportedCt)+(fgla.*0.27);
bulkExportedCm = (fflu.*bulkExportedC)+(fgla.*0.27);

%% Load field data and perform bootstrap re-sampling
data = readtable('NOSAMSallResults.csv'); %Load Data
depInd = find(strcmp(data.Type,'D')==1); % Find Deposit Samples
filtInd = find(strcmp(data.Type,'F')==1); %Find Suspended Sediment Samples
%Select RadioCarbon Data (as Fm)
dataS = data.Fm(depInd); 
dataF = data.Fm(filtInd);
dataSF = data([depInd;filtInd],:);
%Pick the same number of samples for filters and deposits
nSample = 8; %number of samples in each bootstrap
nBoot = 10000; % number of bootstrap calculations
F_ind = randi(length(dataF),[nSample,nBoot]); %random sampling (w/ replacement)
S_ind = randi(length(dataS),[nSample,nBoot]);
bootMeanF = mean(dataF(F_ind)); %mean of each resample
bootMeanS = mean(dataS(S_ind));

clearvars -except tCut WTDexprt WTDstor bootMeanF bootMeanS bulkStoredC ...
    bulkExportedC bulkStoredCm bulkExportedCm bulkStoredCtm ...
    bulkExportedCtm bulkStoredCt bulkExportedCt

%% Plot Results
figure(1)
clf
hold on
plot(bulkStoredC,bulkExportedC,'bo',...
    'markerfacecolor','w','linewidth',2,'markersize',14,'displayname','full')
plot(bulkStoredCt,bulkExportedCt,'ko',...
    'markerfacecolor','w','linewidth',2,'markersize',14,'displayname','truncated')
plot(bulkStoredCm,bulkExportedCm,'bd',...
    'markerfacecolor','w','linewidth',2,'markersize',14,'displayname','full-moraine')
plot(bulkStoredCtm,bulkExportedCtm,'kd',...
    'markerfacecolor','w','linewidth',2,'markersize',14,'displayname','truncated-moraine')
plot(linspace(0.6,1.2),linspace(0.6,1.2),'k-','displayname','1:1')
plot(mean(bootMeanS),mean(bootMeanF),'ro',...
    'markerfacecolor','r','linewidth',2,'markersize',16,'displayname','Field Data')
plot(prctile(bootMeanS,[5,95,95,5]),prctile(bootMeanF,[5,95,5,95]),'r+',...
    'markerfacecolor','r','linewidth',2,'markersize',16,'displayname','95% CI')
text(0.85,0.7,{'T_{cut}=',num2str(tCut)},'fontsize',16)
xlabel('Fluvial Deposit Fm')
ylabel('Suspended Sediment Fm')
set(gca,'fontsize',16)
box on
legend('location','northwest')
axis([0.6 1 0.6 1.2])
