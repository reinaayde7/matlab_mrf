clear all, close all, clc
addpath('./rr_dictionary/DictSimulation_NoSP');
addpath('./dependence');
addpath('./Colormap');
addpath('./Functions_Jesus');
addpath('./OpenSiemensRawData/');

% Path to IRT toolbox
irtPath = 'E:\POSTDOC_UoM\05_MATLAB\fessler-d2310\irt\';
addpath(irtPath); run setup;

% utilities + imagine toolbox
addpath('E:\POSTDOC_UoM\05_MATLAB\rr_utilities\');
codePath = 'E:\POSTDOC_UoM\05_MATLAB\Imagine_old\';
addpath(genpath(codePath))
%% load data
disp('open processed data')
RawData.folder = 'E:\scanner_data\twix_data\241011\';
RawData.name = 'SingleVialAnalysis_trueFISPvFISP_withB0map_21-Oct-2024.mat';
load([RawData.folder RawData.name]);

%% load nist prop
NISTdir = ['E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\_NIST_log\'];
NISTname = ['NIST_log01-Oct-2024.mat'];
load([NISTdir NISTname])
NISTplot(:,1) = NIST.T2array.T1.freemax';
NISTplot(:,2) = NIST.T2array.T2.freemax';

%% note
% up to here, NIST labels and data are in agreement with T1, T2 and b0
% prop.

%%
fprintf('Step 0: Setup Acquisition Parameters\n');

folderTXT = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp';
trtxtfile = [folderTXT '\FISP_TR.txt'];
tetxtfile = [folderTXT '\FISP_TE.txt'];
% fatxtfile = [folderTXT '\FISP_FA_orig.txt'];
fatxtfile = [folderTXT '\FISP_FA_body.txt'];
phtxtfile = [folderTXT '\FISP_PH.txt'];

%% dictionary generator 
isTRUEFISP = 1
log.delay = 3000;

vialN = 7;

log.t1series = NISTplot(vialN,1); % 0:50:500;% %[10:10:1000 1020:20:2000]; %[10:10:2000];
log.t2series = NISTplot(vialN,2); % 0:100:1000;%[2:2:100 105:5:250 250:50:1100]; %[2:2:100 105:5:300 310:10:500 520:20:800 850:50:1100];
% log.t1series = [10:10:1000 1020:20:2000 2050:50:3000]; %[10:10:2000];
% log.t2series = [2:2:100 105:5:300 310:10:500 500:50:1500 1500:100:2000];

if isTRUEFISP
    log.offres = 0; %[-40:1:40]; %[-30:1:-15, -14.5:0.5:-10.5, -10:0.1:10, 10.5:0.5:14.5, 15:1:30]; %[-20:0.1:20]; %[-15:0.5:15];% [-200:10:50 -45:5:45 50:10:200]; %[-200:10:200];
    % log.offres = b0map.avg_b0(vialN);
else
    log.offres = 0;
end

% EstMemSizeDict = size(log.t1series,2)*size(log.t2series,2)*size(log.offres,2)*log.rawinfo.Nex   *8/1024/1024/1024; %in GB
% fprintf('Estimated dictionary size: %.2f GB \n', EstMemSizeDict)
    

if isTRUEFISP
    tic
    [dict,r] = DEBUG_Calculate_MRF_TRUEFISP_DictwithDelays(dTRUEFISP.log.rawinfo,log.t1series,log.t2series,log.offres,...%[0.6:0.1:1.4] for B1
    tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,2);
    log.timeGen = toc;
else
    tic
    [dict,r] = DEBUG_Calculate_MRF_FISP_DictwithDelays(dFISP.log.rawinfo,log.t1series,log.t2series,0,...%[0.6:0.1:1.4] for B1
    tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,2);
    log.timeGen = toc;
end

%%
% dictO=dict;
% dict = j*imag(dict);
% dictnorm = sqrt(sum(dict.*conj(dict)));
%%
dictnorm = sqrt(sum(dict.*conj(dict)));
% vialN = 7;
% 
if isTRUEFISP
    fp = dTRUEFISP.avg_fp(vialN,:);
    % fp = fp*exp(-j*pi/2);
else
    fp = dFISP.avg_fp(vialN,:);
end
[t1,t2,b0] = patternmatchDEBUG(fp,1,r,0,dict,dictnorm,1,1);

%%
dictnorm = sqrt(sum(dict.*conj(dict)));
for vN = 1:14
    if isTRUEFISP
        fp = dTRUEFISP.avg_fp(vN,:);
        fp = fp*exp(-j*pi/2);
    else
        fp = dFISP.avg_fp(vN,:);
    end  
    [t1e(vN),t2e(vN),b0e(vN)] = patternmatchDEBUG(fp,1,r,0,dict,dictnorm,1,1);
end

%%
figure; 
subplot(321)
plot(flipud(NISTplot(:, 1)), '.', MarkerSize=10), hold on
plot(fliplr(t1e), '.', MarkerSize=10);
legend('NIST', 'MRF', 'Location', 'northwest')
title('T1 [ms]'), grid on
subplot(322)
plot(flipud(NISTplot(:, 2)), '.', MarkerSize=10), hold on
plot(fliplr(t2e), '.', MarkerSize=10);
legend('NIST', 'MRF', 'Location', 'northwest')
title('T2 [ms]'), grid on

subplot(323)
plot(NISTplot(:, 1), NISTplot(:, 1)), hold on
plot(NISTplot(:, 1), t1e, '.', MarkerSize=10)
xlabel('NIST');
ylabel('MRF');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('log(T1)'), grid on
subplot(324)
plot(NISTplot(:, 2), NISTplot(:, 2)), hold on
plot(NISTplot(:, 2), t2e, '.', MarkerSize=10);
xlabel('NIST');
ylabel('MRF');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('log(T2)'), grid on

subplot(325)
plot(flipud((t1e'-NISTplot(:, 1))./NISTplot(:, 1)), '.', MarkerSize=10);
ylim([-1;1])
title('error T1 [%]'), grid on
subplot(326)
plot(flipud((t2e'-NISTplot(:, 2))./NISTplot(:,2)), '.', MarkerSize=10);
ylim([-1;1])
title('error T2 [%]'), grid on

%%
saveFolder = 'E:\POSTDOC_UoM\08_Project_MRF\';
save([saveFolder '_rr_simulated_dictionaries\' 'rr_dict_NIST_trueFISP_t1t2B0v101_fov300.mat'],'dict','r','log','-v7.3');
% save([saveFolder '_rr_simulated_dictionaries\' 'rr_dict_NIST_FISP_t1t2v100.mat'],'dict','r','log','-v7.3');