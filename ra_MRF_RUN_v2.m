%
% MRF recon for 0.55T machine 
%
% Rudy Rizzo
% Based on previous scripts by Jesus Fajardo (jesuserf@umich.edu)
% Based on previous scripts by Yun Jiang (yunjiang@med.umich.edu)
clear;clc;close all;
tic
%% setup reconstruction options
% Rudy: we need SVD on otherwise dictionary matching phase is too big
lowrank_SVD             = 1; % 0 :-> 'original' MRF reconstruction
% 1 :-> compress dictionary using SVD and perform the
%       reconstuction in low-rank coefficient images.
isTRUEFISP              = 0;
    doB0regularization  = 0;

doCoilCompression       = 0; 
doPrewhitenData         = 0; %ra: it was set to 1 by Rudy

% --------------------------------------------------------------------
% 250224: PF works if run in SVD-compressed domain on cartesian gridded
% data. Considering the undersampling of MRF, the V^H decompression
% introduces non negligible errors for which PF won't be applicable for DIP
% recon. NB: PF algorithm cannot be applied to SOS traditional MRF in the
% time-domain dimenstion (ie without SVD compression) since there is no
% calibration area for 3D symmetry given projection undersampling of MRF

PartialFourier          = 0;
    mockPF              = 0;
    pfmock              = 6/8;
    TotalSlices         = 48;           % Slices without partial fourier
    pfmock_partitions = ceil(pfmock*TotalSlices);
% --------------------------------------------------------------------

GRAPPAz                 = 0;
    mockGRAPPAz         = 0;
    calibrationLines    = 25;
    Rz                  = 2;

CalculateDict           = 1;    % if 0 => loads a pregen dictionary
doFit                   = 1;    % do not run pm fit
        singleSliceFit  = 1;
        sliceToFit      = 35;

doSave                  = 1;
    tempfilename = ['recon_' date '_FISP'];
    saveFolder = fullfile(pwd, 'results', tempfilename); %'E:\POSTDOC_UoM\08_Project_MRF\';
    % ra:
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder)
    end
doRRdebug               = 0;

doNEXtruncation         = 0;
    truncF              = 0.6;
    dict_truncated_is_saved = 1;

nproj = 48;

%% setup path and irt
%activefilename = matlab.desktop.editor.getActiveFilename;
%[activepath,name,ext] = fileparts(activefilename);
%cd(activepath);
%cd('../MRF_recon_examples/');
addpath('./dictionary/DictSimulation_NoSP/');
% addpath('./dictionary/DictSimulation_SP');
% addpath('./data');
addpath('./dependence');
addpath('./Colormap');
addpath('./Functions_Jesus');
addpath('./OpenSiemensRawData/');

%addpath('E:/POSTDOC_UoM/05_MATLAB/rr_twix_recon') % ra: I commented momentarily because I don't have it now
%addpath('/data/Matlab_Github/MRF_spiral_3D');
%load('./Colormap/T1cm.mat');
%load('./Colormap/T2cm.mat');
% Path to IRT toolbox
irtPath = '/home/ayde/University of Michigan Dropbox/Reina Ayde/MatlabProjects/irt/irt/';
addpath(irtPath); run setup;


% addpath('E:\POSTDOC_UoM\05_MATLAB\rr_utilities\'); % ra: I commented momentarily
codePath = '/home/ayde/University of Michigan Dropbox/Reina Ayde/MatlabProjects/imagine/'; % ra: I commented momentarily
addpath(genpath(codePath)) % ra: I commented momentarily

NoiseDataFileName  = []; % measured noise data

%% Open Siemens 3D Raw data
disp('open Siemens Raw data')
% RawData.folder = 'E:\scanner_data\twix_data\250131\'; %NIST
% RawData.folder = 'E:\scanner_data\twix_data\250203\'; %HV1
% RawData.folder = 'E:\scanner_data\twix_data\250221\'; %HV2 + HV3
% RawData.folder = 'E:\scanner_data\twix_data\250225\'; %NIST
% RawData.folder = 'E:\scanner_data\twix_data\250228\'; %HV4 + NIST
% RawData.folder = 'E:\scanner_data\twix_data\250305\'; %HV5
% RawData.folder = 'E:\scanner_data\twix_data\250307\'; %P1
% RawData.folder = 'E:\scanner_data\twix_data\250328\'; %NIST+optFA
 %RawData.folder = 'C:\Users\ayde\University of Michigan Dropbox\Reina Ayde\brain_1mmISO';
%RawData.folder = 'C:\Users\ayde\University of Michigan Dropbox\Reina Ayde\HV3_forReina'; %P2
RawData.folder = '/home/ayde/University of Michigan Dropbox/Reina Ayde/data_for_Reina_trajectorydesign+brain1x1x3/';


if isTRUEFISP
    % RawData.name = 'meas_MID01026_FID112461_rr_MRF_v3_trueFISP_FAbody.dat'; % NIST
    
    %250221
    % RawData.name = 'meas_MID00564_FID01516_MRFbrain_trueFISP_FA10_np40o20.dat'; % HV3, FA10, Nex600, TRUEFISP
    % RawData.name = 'meas_MID00562_FID01514_MRFbrain_trueFISP_FA15_np40o20.dat'; % HV3, FA15, Nex600, TRUEFISP

    %250225
    % RawData.name = 'meas_MID00083_FID02393_NIST_p40_o20_te18_tr13_nex600_wait3000_FA10_bSSFP.dat'; % NIST, FA10, Nex600, wait3000
    % RawData.name = 'meas_MID00082_FID02392_NIST_p40_o20_te18_tr13_nex600_wait3000_FA15_bSSFP.dat'; % NIST, FA15, Nex600, wait3000

    %250228
    % RawData.name = 'meas_MID00040_FID03061_rr_MRF_NIST_bSSFP_FA15.dat'; % NIST, FA15, Nex600, wait3000
    % RawData.name = 'meas_MID00042_FID03063_rr_MRF_NIST_bSSFP_FA10.dat'; % NIST, FA10, Nex600, wait3000
    % RawData.name = 'meas_MID00026_FID03047_rr_MRF_brain_bSSFP_FA15.dat'; % brain, HV4, FA15, Nex600, wait3000

     % 250228
    % RawData.name = 'meas_MID00017_FID04539_rr_MRF_brain_bSSFP_FA15_600.dat'; % NIST, FA15, Nex600, wait3000

    % 250305
    % RawData.name = 'meas_MID00017_FID04539_rr_MRF_brain_bSSFP_FA15_600.dat'; % HV5, FA15, Nex600, wait3000

    % 250328
    % RawData.name ='meas_MID00201_FID10647_rr_MRF_brain_pSSFP_FA15_600_np36_o11'; %NIST, FA15, Nex600, wait1500
    % RawData.name ='meas_MID00203_FID10649_rr_MRF_brain_pSSFP_FAopt_600_np36_o11'; %NIST, FAopt, Nex600, wait1500
else

    %250131
    % RawData.name = 'meas_MID00026_FID12228_rrMRF_300_PF8_np60_FA15.dat'; % NIST, FA15
    % RawData.name = 'meas_MID00029_FID12231_rrMRF_300_PF8_np60_FA15_tilt.dat'; % NIST, FA15, tilted
    % RawData.name = 'meas_MID00030_FID12232_rrMRF_300_PF8_np60_FA20.dat'; % NIST, FA20
    % RawData.name = 'meas_MID00204_FID12687_rrMRF_300_PF8_np52_FA15_NIST.dat'; % NIST, FA15
    % RawData.name = 'meas_MID00206_FID12689_rrMRF_300_PF8_np52_FA20_NIST.dat'; % NIST, FA20
    
    %250203
    % RawData.name = 'meas_MID00189_FID12672_rrMRF_300_PF8_np52_FA20.dat'; % HV, FA20
    % RawData.name = 'meas_MID00185_FID12668_rrMRF_300_PF8_np52_FA15.dat'; % HV, FA15
    % RawData.name = 'meas_MID00188_FID12671_rrMRF_300_PF8_np48_tilted_FA15.dat'; % HV, FA15, titled
    
    %250221
    RawData.name = 'meas_MID00547_FID01499_MRFbrain_FISP_FA15_np40o20_nex600_FA15.dat'; % HV2, FA15, Nex600
    % RawData.name = 'meas_MID00548_FID01500_MRFbrain_FISP_FA15_np40o20_nex1000_FA15.dat'; % HV2, FA15, Nex1000
    % RawData.name = 'meas_MID00549_FID01501_MRFbrain_FISP_FA10_np40o20_nex600_FA10.dat'; % HV2, FA10, Nex600
    % RawData.name = 'meas_MID00550_FID01502_MRFbrain_FISP_FA10_np40o20_nex1000_FA10.dat'; % HV2, FA10, Nex1000
    % RawData.name = 'meas_MID00561_FID01513_MRFbrain_FISP_FA15_np40o20.dat'; % HV3, FA15, Nex600
      
    %250225
    % RawData.name = 'meas_MID00077_FID02387_NIST_p40_o20_te18_tr14_nex600_wait3000_FA15.dat'; % NIST, FA15, Nex600, wait3000
    % RawData.name = 'meas_MID00078_FID02388_NIST_p40_o20_te18_tr14_nex600_wait2000_FA15.dat'; % NIST, FA15, Nex600, wait2000
    % RawData.name = 'meas_MID00080_FID02390_NIST_p40_o20_te18_tr14_nex600_wait1500_FA15.dat'; % NIST, FA15, Nex600, wait1500
    % RawData.name = 'meas_MID00081_FID02391_NIST_p40_o20_te18_tr14_nex600_wait1000_FA15.dat'; % NIST, FA15, Nex600, wait1000
    % RawData.name = 'meas_MID00084_FID02394_NIST_p40_o20_te18_tr14_nex600_wait3000_FA10.dat'; % NIST, FA10, Nex600, wait3000

    %250228
    % RawData.name = 'meas_MID00024_FID03045_rr_MRF_brain_FISP_FA15.dat'; % HV4, FA15, Nex600, wait3000
    % RawData.name = 'meas_MID00025_FID03046_rr_MRF_brain_FISP_FA15_wait1500.dat'; % HV4, FA15, Nex600, wait1500
    % RawData.name = 'meas_MID00041_FID03062_rr_MRF_NIST_FISP_FA10.dat'; % NIST, FA10, Nex600, wait3000
    % RawData.name = 'meas_MID00038_FID03059_rr_MRF_NIST_FISP_FA15.dat'; % NIST, FA15, Nex600, wait3000

    %250305
    % RawData.name = 'meas_MID00014_FID04536_rr_MRF_brain_FISP_FA15_600.dat'; % HV5, FA15, Nex600, wait3000
    
    %250307
    % RawData.name = 'meas_MID00671_FID05198_rr_MRF_brain_FISP_FA15_600_W1500.dat'; % P1, FA15, Nex600, wait1500

    %250328
    % RawData.name = 'meas_MID00200_FID10646_rr_MRF_brain_FISP_FA15_600_np36_o11.dat'; % NIST, FA15, Nex600, wait1500, w1500
    % RawData.name = 'meas_MID00202_FID10648_rr_MRF_brain_FISP_FAopt_600_np36_o11.dat'; % NIST, FAopt, Nex600, wait1500, w1500
    %RawData.name = 'meas_MID00026_FID03728_rrMRF_f300_1mmISO_FISP_FA15.dat'; % BRAIN iso1mm 0.55T, FAopt, Nex600, wait3000, w1500

    %rudy's data
    %RawData.name = 'meas_MID00561_FID01513_MRFbrain_FISP_FA15_np40o20.dat'; % P2, FA15, Nex600, wait1500

end
RawDataFileName = fullfile(RawData.folder,RawData.name);

[path,name,ext] = fileparts(RawDataFileName);
[raw,noise,log.rawinfo,~] = loadSiemensRawData(RawDataFileName);
disp(['raw data size: ' num2str(size(raw))])

%% retrospective postprocessing
if doNEXtruncation %Image timepoints truncation
    t = size(raw,4)*truncF;
    raw = raw(:,:,:,1:t);
    log.rawinfo.Nex = t;
    log.rawinfo.nframes = t;
end

if PartialFourier %Mocking of Partial Fourier acquisition
    if mockPF
        raw = raw(:,:,1:pfmock_partitions,:);
        disp(['mocked PF raw data size: ' num2str(size(raw))])
    end
end
if ~mockPF
    TotalSlices = size(raw,3);
end

% if GRAPPAz
%     if mockGRAPPAz
%         if Rz == 2;
%             raw_noG = raw;
%             smp_mask = [1 0 1 0]';
%             mask = repmat(smp_mask,[ceil(60/4),1]);
%             mask = mask(1:TotalSlices);           
%             mask(TotalSlices/2+1 - (calibrationLines-1)/2:TotalSlices/2+1 + (calibrationLines-1)/2) = 1;
%             raw(:,:,not(mask),:) = 0;
%         end
%     end
%     maskZ = zeros(size(raw,3), size(raw,4));
%     maskZ(abs(raw(1,1,:,:))>eps) = 1;
%     figure,imagesc(maskZ); axis square;
%     xlim([0 20]);drawnow; xlabel('Time points'), ylabel('Kz')
% end

%% Prewhiten data
if doPrewhitenData

    disp('Load raw noise')
    if ~isempty(NoiseDataFileName)
        [noise_raw] = loadSiemensRawData(NoiseDataFileName);
        noise = noise_raw;
    end
    
    savedir = RawData.folder;

    disp('Prewhiten data')

    disp('calc noise decorrelation ')
    if ~isempty(noise)
        [dmtx] = calculate_noise_decorrelation_mtx(permute(noise,[1,3,4,2]));
    
    else
        dmtx = eye(log.rawinfo.ncoils);
    end

    if ndims(raw) == 4,
        raw = reshape(raw,[size(raw,1),size(raw,2),size(raw,3),1,size(raw,4)]);
    end

    raw = permute(raw,[1,3,4,5,2]);
    raw = apply_noise_decorrelation_mtx(raw,dmtx);
end

%% Coil Compression
% using SCC coil compression
if doCoilCompression
    fprintf('Coil Compression\n')
    raw = pcaCoilCompress(raw,0.90);
end
% ra: I think the following should be inside the if
% raw = permute(raw,[1,5,2,3,4]);
% raw = reshape(raw,[size(raw,1),size(raw,2),size(raw,4),size(raw,3),size(raw,5)]);
% disp(['raw data size (coil compressed): ' num2str(size(raw))])

%% Setup acquisition parameter
fprintf('Step 0: Setup Acquisition Parameters\n');

folderTXT = '/home/ayde/University of Michigan Dropbox/Reina Ayde/HV3_forReina/rudy_tom_MRF_bssfp';
trtxtfile = [folderTXT '/FISP_TR.txt'];
tetxtfile = [folderTXT '/FISP_TE.txt'];
%fatxtfile = [folderTXT '\FISP_FA_Body.txt'];
fatxtfile = [folderTXT '/FISP_FA_body_x1_5.txt'];
%fatxtfile = [folderTXT '\FISP_FA_body_x2.txt'];
%fatxtfile = [folderTXT '\FISP_FA_Body_shifted.txt'];
% fatxtfile = [folderTXT '\FISP_FAopt_nex600_TR14_TE1_8_SC_EPG_GMWMort.txt'];
% fatxtfile = [folderTXT '\pSSFP_FAopt_nex600_TR13_TE1_8_SC_EPG_GMWMort.txt'];
phtxtfile = [folderTXT '/FISP_PH.txt'];


log.FA = importdata(fatxtfile);
log.tr0 = log.rawinfo.TR(1)/1000/1000;
log.TR = importdata(trtxtfile)*1e-6 + log.tr0;
log.TE = importdata(tetxtfile);
log.PH = importdata(phtxtfile);

% plot the parameters
figure('name','Acquisition Parameters');
subplot(211);plot(log.FA(1:log.rawinfo.Nex,1),'LineWidth',2);
title('Flip Angles'); xlabel('Time Points');ylabel('Flip Angles (degree)');
subplot(212);plot(log.TR(1:log.rawinfo.Nex,1)*1000,'LineWidth',2);
title('Repetition Time'); xlabel('Time Points');ylabel('TR (ms)');

%% Dictionary handler: 
% 1. Calculate or import
% 2. SVD compressed or NOT
% dictionary_name = 'dict_trueFISP_TE1_8_TR13_FA15_Nex600_b0_0_wait1500';
%dictionary_name = 'dict_FISP_TE2_8_TR17_FAbody10_Nex1000_delay0_SVD';
% dictionary_name = 'dict_FISP_TE2_09_TR11_FAbody15_Nex600_wait1500_SVD';
%brain 1 x1 x3
dictionary_name = 'dict_FISP_TE1_8_TR14_FAbody15_Nex600'
if CalculateDict
    fprintf('Step 1: Calculate Dictionary\n');
    log.delay = 3000; % what is log.delay??
    %log.t1series = [10:10:1000 1020:20:2000 2050:50:3000]; %[10:10:2000];
    %log.t2series = [2:2:100 105:5:300 310:10:500 520:20:800 850:50:1500 1600:100:2000];
    %new fine dict for NIst
    %log.t1series = [10:10:1000 1000:2:1500 1520:20:2000 2050:50:3000]; %10:10:1000 1020:20:2000 
    %log.t2series = [10:10:1000 1002:2:1500 1520:20:2000 2050:50:3000];
    log.t1series = [10:2:1000 1020:20:2000 2050:50:3000];
    log.t2series = [2:2:100 105:5:300 310:10:500 520:20:800 850:50:1500 1600:100:2000];
    if isTRUEFISP
        % B0 brain = [-100:5:-45 -40:1:40 45:5:100]
        % B0 nist = [-40:1:40]
        log.offres = [0]; %[-100:5:-45 -40:1:40 45:5:100]; 
    else
        log.offres = [0];
    end

    EstMemSizeDict = size(log.t1series,2)*size(log.t2series,2)*size(log.offres,2)*log.rawinfo.Nex   *8/1024/1024/1024; %in GB
    fprintf('Estimated dictionary size: %.2f GB \n', EstMemSizeDict)
    

    if isTRUEFISP
        tic
        [dict,r] = Calculate_MRF_TRUEFISP_DictwithDelays(RawDataFileName,log.rawinfo,log.t1series,log.t2series,log.offres,...%[0.6:0.1:1.4] for B1
        tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,2);
        log.timeGen = toc;
        fprintf('Saving TRUEFISP dictionary \n')
        savedictFolder = fullfile(saveFolder, '_rr_simulated_dictionaries\');
        if ~exist(savedictFolder, 'dir')
            mkdir(savedictFolder)
        end
        save([savedictFolder dictionary_name '.mat'],'dict','r','log','-v7.3');
    else
        tic
        [dict,r] = Calculate_MRF_FISP_DictwithDelays(RawDataFileName,log.rawinfo,log.t1series,log.t2series,0,...%[0.6:0.1:1.4] for B1
        tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,2);
        log.timeGen = toc;
        fprintf('Saving FISP dictionary \n')
        savedictFolder = fullfile(saveFolder, '_rr_simulated_dictionaries\');
        if ~exist(savedictFolder, 'dir')
            mkdir(savedictFolder)
        end
        save([savedictFolder dictionary_name '.mat'],'dict','r','log','-v7.3');
    end

    fprintf('Compressing dictionary \n')
    [dictSVD,Vc,S] = svd_compress_dictionary(squeeze(dict),1e-2); %Rudy: 1e-2 originally, with offres == 0
    fprintf('Saving SVD compressed dictionary \n')
    savedictFolder = fullfile(saveFolder, '_rr_simulated_dictionaries\');
    if ~exist(savedictFolder, 'dir')
        mkdir(savedictFolder)
    end
    save([savedictFolder dictionary_name '.mat'],'dict','r','log','-v7.3');

    if lowrank_SVD
        dict = dictSVD;
    end

else
    fprintf('Step 1: Load Dictionary\n');   
    if lowrank_SVD           
            if ~dict_truncated_is_saved && doNEXtruncation 
                fprintf('>> non-compressed version  -> to RE-compress...\n');
                % we need to re-compress it! import non-compressed dictionary!
                 %dfile = load(['C:\Users\ayde\University of Michigan Dropbox\Reina Ayde\HV3_forReina\' dictionary_name '.mat'],'dict','r','log');
                 dfile = load(['C:\Users\ayde\University of Michigan Dropbox\Reina Ayde\data_for_Reina\dict_FISP_TE1_8_TR14_FAbody15_Nex600_SVD.mat'])
                 dict = dfile.dict; r = dfile.r; dict_log = dfile.log; 
                 dict = dict(1:log.rawinfo.nframes,:); %truncation
                 fprintf('RE-Compressing dictionary... \n')
                 [dictSVD,Vc,S] = svd_compress_dictionary(squeeze(dict),1e-2); %Rudy: 1e-2 originally, with offres == 0
                 dict=dictSVD;
                 fprintf('Saving SVD compressed dictionary \n')
                 save([saveFolder '_rr_simulated_dictionaries\' dictionary_name '_t' num2str(t) '_compressed.mat'],'dictSVD','r','log','Vc','-v7.3');
            else
                fprintf('>> compressed version... \n');
                %dfile = load(['E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\' dictionary_name '.mat'],'dictSVD','r','log','Vc');
                %dfile = load(['C:\Users\ayde\University of Michigan Dropbox\Reina Ayde\HV3_forReina\' dictionary_name '.mat'],'dictSVD','r','log','Vc');
                dfile = load(['C:\Users\ayde\University of Michigan Dropbox\Reina Ayde\brain_1mmISO\Brain_HV_251217_1mmISO_s5_f300_Nex1000_FAbody15_TE2_8_TR17\dictCompressed.mat'])
                dict = dfile.dictCompressed; r = dfile.r;  Vc = dfile.Phi; %dict_log = dfile.log;
                clearvars dfile
            end
    else    
            dfile = load(['C:\Users\ayde\University of Michigan Dropbox\Reina Ayde\data_for_Reina\' dictionary_name '.mat'],'dict','r','log');
            dict = dfile.dict; r = dfile.r; dict_log = dfile.log; 
            clearvars dfile        
    end
end

if lowrank_SVD
    %dip setup
    dip.dictCompressed = dict;
    dip.r = r;
    dip.Phi = Vc;
end
%% load spiral trajectory
fprintf('Step 2: Prepare NUFFT \n');

%load('./data/Spiral_Traj.mat','kxall','kyall');
SpiralMeasFileName.folder = '/home/ayde/University of Michigan Dropbox/Reina Ayde/MatlabProjects/rudy_matlab_script/rudy_tom_MRF_bssfp/';
% SpiralMeasFileName.name = 'spiral_055T_MRF_FOV400mtx400_jesus_241118.mat'; %original from Jesus but acquired by Rudy
% SpiralMeasFileName.name = 'SpiralHT.mat';
% SpiralMeasFileName.name = 'spiral_055T_MRF_FOV400mtx400_rudy_241118.mat'; %optimized 400/400 setting from Rudy
SpiralMeasFileName.name = 'spiral_055T_MRF_FOV300mtx300_rudy_241118.mat'; %optimized 300/300 setting from Rudy
% SpiralMeasFileName.name = 'spiral_055T_MRF_FOV250mtx250_np96_M0_rudy_tra_250407'; %250/250 with np96
SpMeasFileName = fullfile(SpiralMeasFileName.folder,SpiralMeasFileName.name);
% [kxall,kyall] = GradTrajMeas(SpMeasFileName);

load(SpMeasFileName)

% Prepare NUFFT

% --------------------------------------------------------------------
% JESUS
% --------------------------------------------------------------------
% Calcuate kmax and the corresponding point
[kmax_measured,uplimit]=max(sqrt(kxall(:,1).^2+kyall(:,1).^2));
% remove additional points
adcpad=10; % the ADC acquires 20 points before and after spiral: NOTE that the way I (rudy) measure spirals, has a pad of 10 points only! this affects the density compensation wrongly, otherwise!
% account for the adcpad
kx = kxall(1+adcpad:end-adcpad,:);
ky = kyall(1+adcpad:end-adcpad,:);
kx = kx(1:uplimit,:)/10; %1/cm->1/mm
ky = ky(1:uplimit,:)/10; %1/cm->1/mm

%plot
figure;plot(kxall(1:adcpad,1),kyall(1:adcpad,1),'r',...
    kxall(adcpad+1:uplimit,1),kyall(adcpad+1:uplimit,1),'b',...
    kxall(uplimit+1:end,1),kyall(uplimit+1:end,1),'r','LineWidth',2);
axis([-5,5,-5,5]);axis square;
title('1 spiral arm'); xlabel('kx (1/cm)');ylabel('ky (1/cm)');



% Rudy: check for DIP
% the DIP has a 4 depth Unet that squeeze by factor 2 the input size. The
% input size is twice the original size, therefore here we only look for a
% factor 2 three times, 2^3 = 8

% for example: FOV:300 -> Unet layers: 600->300->150->75->36.5(error)
% we need to round it up to 304!
fprintf(['Check FOV for DIP - original: ' num2str(log.rawinfo.matrix) ' mm \n'])
if (log.rawinfo.matrix/8) ~= 0
    f = ceil(log.rawinfo.matrix/8);
    log.rawinfo.matrix = f*8;
end
fprintf(['Check FOV for DIP - corrected: ' num2str(log.rawinfo.matrix) ' mm \n'])

% Normalize the kspace
resolution = log.rawinfo.fieldofview(1,1)/log.rawinfo.matrix; % spatial resolution
kmax = 1/(2*resolution);
resolution_measured = 1/(2*(kmax_measured/10));
% normalized spiral trajector to [-0.5 0.5]

kx = kx./(2*kmax_measured/10);
ky = ky./(2*kmax_measured/10);
k = [kx(:) ky(:)];

mask=true([log.rawinfo.matrix log.rawinfo.matrix]);
N=size(mask);
nufft_args = {N, [5 5], 2*N, N/2, 'table', 2^12, 'minmax:kb'};

% Calculate DCF
G_fid_dcf = Gmri(k, mask, 'fov', N(1), 'basis', {'dirac'}, 'nufft', nufft_args);
dcf_fid = reshape(abs(mri_density_comp_v2(k, 'pipe', 'G', G_fid_dcf.arg.Gnufft)),uplimit,size(kx,2));
%dcf_fid = reshape(abs( mri_density_comp(k, 'pipe', 'G', G_fid_dcf.arg.Gnufft)),uplimit,size(kx,2));

adcpad = 20; %MRF data are padded by 20 points

% dip setup
j=1;
for i=1:log.rawinfo.nframes   
    if j>nproj
        j=1;
    end
        dip.idproj(i) = j;
        dip.k(:,i) = complex(kx(:,j),ky(:,j));  
        dip.w(:,i) = dcf_fid(:,j);
        j = j+1;
end
    
% Calculate nufft
if lowrank_SVD,
    k= double([kx(:) ky(:)]);
    G_gridding = Gmri(k, mask, 'fov', log.rawinfo.fieldofview(1), 'basis', {'dirac'}, 'nufft', nufft_args);
else
    G_gridding = cell(size(kx,2),1);

    for s =1:size(kx,2), %Rudy: build up a G_gridding function for each spiral arm (= n. 48)
        k = double([kx(:,s) ky(:,s)]);
        G_gridding{s} = Gmri(k, mask, 'fov', log.rawinfo.fieldofview(1), 'basis', {'dirac'}, 'nufft', nufft_args);
        %nufft_st{ii,1} = nufft_init(2*pi*tempk,N,J,K,N/2,'minmax:kb');% Fessler's toolbox require the k in radian --YunJiang-01.27.16
    end
end

size(raw)
raw = raw(:,:,:,1:log.rawinfo.Nex); % for saving time, we will only process first 1000 images
% Remove adc padding in the raw data
raw = raw(1+adcpad:end-adcpad,:,:,:); % for 1proj
raw = raw(1:uplimit,:,:,:);
raw_orig = raw;

%% phase rolling to center FOV for inplane shifts (estimated at the scanner in mm)
delta_x = 0;%-104; %read [mm] - vertical offset with images convention in matlab (- moves up, + moves down)
delta_y = 0;%+56; %phase [mm] - horizonatl offset with images convention in matlab (- moves left, + moves right)
% phase_corr = exp(1i*(kx.*delta_x + ky.*delta_y));
phase_corr = exp(-1i*(kx.*2*pi*delta_x + ky.*2*pi*delta_y));
kk=1;
for i =1:size(raw_orig,4)
    s(:,:,:,i) = raw_orig(:,:,:,i).*phase_corr(:,kk);
    kk=kk+1;
    if kk>nproj;
        kk=1;
    end
end
raw = s;
%%
if GRAPPAz
    if mockGRAPPAz
        if Rz == 2
            raw_noG = raw;
            smp_mask = [1 0 1 0]';
            mask_G = repmat(smp_mask,[ceil(60/4),1]);
            mask_G = mask_G(1:TotalSlices);       
            mask_CAL = TotalSlices/2+1 - (calibrationLines-1)/2:TotalSlices/2+1 + (calibrationLines-1)/2;
            % mask_G(mask_CAL) = 1;
            raw(:,:,not(mask_G),:) = 0;
        end
    end
    maskZ = zeros(size(raw,3), size(raw,4));
    maskZ(abs(raw(1,1,:,:))>eps) = 1;
    figure,imagesc(maskZ); axis square;
    xlim([0 20]);drawnow; xlabel('Time points'), ylabel('Kz')
end

%%
% %load extra calib
% CalData.name ='meas_MID00020_FID04542_rr_MRF_brain_FISP_calFA50_wait100.dat';
% CalDataFileName = fullfile(RawData.folder,CalData.name);
% 
% [path,name,ext] = fileparts(CalDataFileName);
% [cal,noise,log.calinfo,~] = loadSiemensRawData(CalDataFileName);
% 
% cal = cal(1+adcpad:end-adcpad,:,:,:); 
% cal = cal(1:uplimit,:,:,:);
%%

if GRAPPAz
    codePath = 'E:\POSTDOC_UoM\05_MATLAB\parallel-master-d2404';
    addpath(genpath(codePath))
    raw_G = zeros(size(raw));
    % 
    % padS = 3;
    % raw_pad = zeros([size(raw,1)+padS*2,size(raw,2), size(raw,3), size(raw,4)]);
    % raw_pad(padS+1:end-padS,:,:,:) = raw;
    % 
    if mockGRAPPAz
        disp('mocked GRAPPA recon...')
        f = waitbar(0, 'mocked GRAPPA recon...'); tic;
        % a=1;
        for tp = 1:1%size(raw,4)
            waitbar(tp/size(raw,4), f, sprintf('Time Point: %d/%d', tp, size(raw,4)));
            
        %--------------------------------------------------------------------
            % Spiral through-time GRAPPA 
            % % % data = permute(raw_noG(:,:,:,tp),[1,3,2]);    
            % % % data_u = data;
            % % % data_u(:,2:2:end,:) = 0;
            acs_lines = 13:37;
            data_u = permute(raw(:,:,:,tp),[1,3,2]);
        
            pad=3;
            kx = 3;
            ky=2;
            
            cal_idx=1;
            if tp<nproj+1
                doCal=1;
                wsc=[];  
                v=cal_idx:nproj:600
                cal = permute(raw_noG(:,:,acs_lines,[cal_idx:nproj:end]),[1,3,2,4]);     
                [data_r, ws{cal_idx}] = RR_spiralGRAPPA_MRF_parundersmp(data_u,pad,kx,ky,doCal,cal,wsc);
                cal_idx=cal_idx+1;
            else %calibration already done!
                doCal=0;
                cal=[];
                [data_r, ~] = RR_spiralGRAPPA_MRF_parundersmp(data_u,pad,kx,ky,doCal,cal,ws{cal_idx});
                cal_idx=cal_idx+1;
            end
        
            if cal_idx>nproj
                cal_idx=1;
            end    
            
            raw_G(:,:,:,tp) = permute(data_r,[1,3,2]);
           
        end
        close(f); time = toc; GRAPPA_time = time2clock(time);
    end
    % raw_G = permute(raw_G,[1,3,2,4]);
    % %plotting for debug
    % idx = 1:50;
    idx = 50:150;
    % idx = 150:300;
    % idx = 500:1000;
    % idx = 1000:2000;
    % idx = 2000:3000;

    coil = 10;
    timep = 1;
    figure, 
    subplot(222),imagesc(angle(squeeze(raw(idx,coil,:,timep)))), title('raw undersmp'), colorbar, 
    subplot(221),imagesc(angle(squeeze(raw_G(idx,coil,:,timep)))), title('raw GRAPPA'), colorbar, ylabel('RO')
    subplot(223),imagesc(angle(squeeze(raw_noG(idx,coil,:,timep)))), title('raw fullysmp'), colorbar, ylabel('RO'),xlabel('Kz')
    subplot(224),imagesc( angle( squeeze(raw_noG(idx,coil,:,timep)) - squeeze(raw_G(idx,coil,:,timep)))), title('raw fullysmp - raw GRAPPA'), colorbar, xlabel('Kz')
    title('kspace KZ - ABS')

    figure, 
    subplot(222),imagesc(abs(squeeze(raw(idx,coil,:,timep)))), title('raw undersmp'), colorbar, 
    subplot(221),imagesc(abs(squeeze(raw_G(idx,coil,:,timep)))), title('raw GRAPPA'), colorbar, ylabel('RO')
    subplot(223),imagesc(abs(squeeze(raw_noG(idx,coil,:,timep)))), title('raw fullysmp'), colorbar, ylabel('RO'),xlabel('Kz')
    subplot(224),imagesc( abs( squeeze(raw_noG(idx,coil,:,timep)) - squeeze(raw_G(idx,coil,:,timep)))), title('raw fullysmp - raw GRAPPA'), colorbar, xlabel('Kz')
    title('kspace KZ - PHASE')

    figure, 
    subplot(222),imagesc(abs(fft(squeeze(raw(idx,coil,:,timep)),[],2))), title('raw undersmp'), colorbar,
    subplot(221),imagesc(abs(fft(squeeze(raw_G(idx,coil,:,timep)),[],2))), title('raw GRAPPA'), colorbar, ylabel('RO')
    subplot(223),imagesc(abs(fft(squeeze(raw_noG(idx,coil,:,timep)),[],2))), title('raw fullysmp'), colorbar, ylabel('RO'),xlabel('Z')
    subplot(224),imagesc( abs( fft(squeeze(raw_noG(idx,coil,:,timep)),[],2) - fft(squeeze(raw_G(idx,coil,:,timep)),[],2) )), title('raw fullysmp - raw GRAPPA'),  colorbar, xlabel('Z')
    title('image space Z - ABS')

    
    % raw=raw_G;
    % raw(:,:,acs_lines,:) = raw_noG(:,:,acs_lines,:);  %restore the acquired ACS lines where they belong
end
   
%%
if lowrank_SVD
    rawsvd = zeros(size(raw,1),nproj,size(raw,3),size(raw,2),size(dict,1)); %48 spiral rotations
    disp(size(raw))
    disp('SVD-compressing raw data...')
    f = waitbar(0, 'SVD-compressing raw data...'); tic;
    % idtxtfile = [folderTXT '\FISP_ID.txt'];
    % projID = textread(idtxtfile,'%d\b');
    % idproj = projID(1:size(raw_orig,4))+1;

    projID = repmat(1:nproj, 1, ceil(size(raw_orig,4)/nproj)); %Rudy 2504 - linear spiral sampling independent from txt file
    idproj = projID(1:size(raw_orig,4))';
    
    for nn = 1:size(raw,3)
        waitbar(nn/size(raw,3), f, sprintf('Partition: %d/%d', nn, size(raw,3)));       

        raw_svd_1s = lowrank_2Dksp( squeeze(raw(:,:,nn,:)) ,idproj, Vc, nproj);
        rawsvd(:,:,nn,:,:) = raw_svd_1s;

    end
    close(f); time = toc; % ra I commented the following cause I don't have this function time2clock %SVD_comp_time = time2clock(time);
    raw = rawsvd;
    [SpiralReadout,~,NumSlices,NumCoils,NumFrames] = size(raw);
    disp('Data SVD-compressed!')
    disp(['data size (SVD-compressed): ' num2str(size(raw))])
else
    [SpiralReadout,NumCoils,NumSlices,NumFrames] = size(raw);
end

% % % %% saving for GRIDDING + PF on full Nex to be run on cluster
% % % tempfilename = ['data_GRIDDING_PF_cluster_' date '_' RawData.name(1:13)];
% % % save([RawData.folder tempfilename], 'NumFrames', 'TotalSlices', 'NumCoils', 'log', 'G_gridding', 'kx', 'dcf_fid', 'raw', 'mask', '-v7.3')

%% NUFFT gridding
disp('K-space GRIDDING...')
% Perform NUFFT

tempindex =  zeros(1, NumFrames);
img_XYkz = squeeze(single(zeros(log.rawinfo.matrix,log.rawinfo.matrix,TotalSlices,NumCoils,NumFrames)));
tic
f = waitbar(0, 'K-space GRIDDING...'); tic;
for iFrame= 1:NumFrames %
    waitbar(iFrame/NumFrames, f, sprintf('Time-Point OR SVD-frame: %d/%d', iFrame, NumFrames));
    if lowrank_SVD %idproj becomes the SVD dimension (5) and you use only 1 gridding element instead of all 48
        G = G_gridding;
        dcf = dcf_fid;
        for inc=1:NumCoils
            for p=1:size(raw,3)
                k_temp = squeeze(raw(:,:,p,inc,iFrame));  % fullly sampled data
                image_temp = G'*(dcf(:).*k_temp(:));
                img_XYkz(:,:,p,inc,iFrame) = embed(image_temp,mask);
            end
        end
    else
        tempindex(:,iFrame) = mod(iFrame - 1,size(kx,2))+1; %variable that matches Nex to spiral arm
        G = G_gridding{tempindex(:,iFrame)};
        dcf = dcf_fid(:,tempindex(:,iFrame));
        for inc=1:NumCoils
            for p=1:size(raw,3)
                k_temp = squeeze(raw(:,inc,p,iFrame));  % single shot
                image_temp = G'*(dcf(:).*k_temp(:));
                img_XYkz(:,:,p,inc,iFrame) = embed(image_temp,mask);
            end
        end
    end
end
ksp_cartesian = transformImageToKspace(transformImageToKspace(img_XYkz,1),2);

close(f); time = toc; %ra I commented the follwing % gridding_time = time2clock(time);
disp('K-space GRIDDED!')
disp(['data size (gridded): ' num2str(size(ksp_cartesian))])
% 60partition, 600 frames (NO SVD) = 48min

%%
if PartialFourier
    disp('Partial Fourier reconstruction...')
    % ksp_cartesian_full = ksp_cartesian;
    f = waitbar(0, 'run POCS...'); tic;
    tic
    for svd_i = 1:size(raw,5)
        waitbar(svd_i/size(raw,5), f, sprintf('Time-Point OR SVD-frame: %d/%d', svd_i, size(raw,5)));
        ksp = permute(squeeze(ksp_cartesian(:,:,:,:,svd_i)), [4,1,2,3]);
        iter = 20;
        watchProgress = 0;
        verbose = 1;
        [~, kspFull] = pocs( ksp, iter, watchProgress, verbose );
        pippo(:,:,:,:,svd_i) = kspFull;
        % pippou(:,:,:,:,svd_i) = ksp;
        clearvars ksp kspFull
    end
    ksp_cartesian_full = permute(pippo,[2,3,4,1,5]);

    figure,
    subplot(221),  imagesc(abs(squeeze(ksp_cartesian(160,:,:,1,1)))), colorbar, title('Partial Fourier'), ylabel('Kx')
    subplot(222),  imagesc(abs(squeeze(ksp_cartesian_full(160,:,:,1,1)))), colorbar, title('Fully sampled')
    subplot(223),  imagesc(angle(squeeze(ksp_cartesian(160,:,:,1,1)))), colorbar, ylabel('Kx'), xlabel('Kz')
    subplot(224),  imagesc(angle(squeeze(ksp_cartesian_full(160,:,:,1,1)))), colorbar, xlabel('Kz')

    close(f); time = toc; PF_POCS_time = time2clock(time);
    disp('Partial Fourier reconstruction DONE!')

    ksp_cartesian_orig = ksp_cartesian;
    ksp_cartesian = ksp_cartesian_full;
end

%% CSM estimation

image_uncombined = transformKspaceToImage(transformKspaceToImage(transformKspaceToImage(ksp_cartesian,3),2),1);
image_uncombined = flip(image_uncombined, 3); % This is to match the B0 and B1 ordering - flip in Kz domain

fprintf('3D CSM estimation... \n');
%average through Nex (coils do not depend on the eigenvalue/ex)
image_proxy_coil = sum(image_uncombined,5);

tic
smoothing = 200; %200 by rudy
chunks = size(ksp_cartesian,3)/6;
csm = ismrm_estimate_csm_walsh_3D(single(image_proxy_coil), smoothing, chunks);
time = toc; %csm_estimate_time = time2clock(time);

% image_proxy_combined=squeeze(sum( conj(csm).*image_proxy_coil, 4 ));
for eig = 1:size(image_uncombined,5)
    image_combined_3D(:,:,:,eig) = squeeze(sum( conj(csm).*squeeze(image_uncombined(:,:,:,:,eig)), 4 ));
end

%imagine(image_combined_3D) %visual check
%imagine(csm)

proxy_img = sum(image_combined_3D,3);
dip.coilmap = csm;

%% DIP saving
destFolder = fullfile(saveFolder, 'DIP/');
Folder = 'FISP/';
what = 'Brain_1x1x3_T1_2ms1000';
%what = 'NIST_T2L_250328';
% what = 'Brain_P1_250307_tilt';

extraN = 'f300_Nex1000_FA15';
if ~PartialFourier
    tmp = flip(transformKspaceToImage(raw_orig,3),3);
    %tmp = transformKspaceToImage(raw_orig,3);
else
    tmp = flip(transformKspaceToImage(rawd_o,3));
end
for nn=1:TotalSlices

    subFolderName =[what '_s' num2str(nn) '_' extraN];
    mkdir(fullfile(destFolder,Folder,subFolderName))
    
    % coilmaps
    coilmap = squeeze(dip.coilmap(:,:,nn,:));
    save([destFolder Folder subFolderName '/coilmap.mat'], 'coilmap', '-v7.3'); %[Nx, Ny, coils]
    
    %data
    DATA = squeeze(tmp(:,:,nn,:)); % [4004 4 38 600]
    save([destFolder Folder subFolderName '/DATA.mat'],'DATA','-v7.3'); %[Nro, coils, Nex];
    
    dictCompressed = dict;
    Phi = Vc;
    save([destFolder Folder subFolderName '/dictCompressed.mat'], 'Phi', 'dictCompressed', 'r', '-v7.3'); %Phi: [Nex,EIG], dictCompressed: [EIG,DicSIZE], r: [DicSIZE, DicVAR]
    
    clearvars k w
    k(:,1,:) = dip.k;
    w(:,1,:) = dip.w;
    save([destFolder Folder subFolderName '/trajectory.mat'], 'k', 'w', '-v7.3'); %k: [Nro,1,Nex], w:[Nro,1,Nex]
       
end %slices
disp('DIP savings DONE!')

%% truncate image to optimize the fit to only area of interest

% phantom 250
% s_vert = 30:240;
% s_hor = 50:256;

% phantom 300
% s_vert = 50:250;
% s_hor = 80:270;

% brain 300
s_vert = 40:260;
s_hor = 60:300;

a = squeeze(image_combined_3D(s_vert,s_hor,:,1));

%imagine(a)
mask_fit = mask-1;
mask_fit(s_vert,s_hor) = 1;


%% LLR settings
addpath(genpath('/home/ayde/University of Michigan Dropbox/Reina Ayde/MatlabProjects/CardiacMRF_2D/Reconstruction Code/Gridding'));
addpath(genpath('/home/ayde/University of Michigan Dropbox/Reina Ayde/MatlabProjects/CardiacMRF_2D/Reconstruction Code/Dictionary/LowRankRecon'));
addpath(genpath('/home/ayde/University of Michigan Dropbox/Reina Ayde/MatlabProjects/CardiacMRF_2D/CardiacMRF_2D/Reconstruction Code/Miscellaneous'));
raw_orig_Z = flip(transformKspaceToImage(raw_orig,3),3);
readOSFactor=1;
numSpiralArms=nproj;
nr = size(raw_orig,1);
nex = size(raw_orig,4);
numCoils = size(raw_orig,2); 
% idproj = projID(1:nex)+1;
w = repmat(dcf_fid, [1,ceil(nex/numSpiralArms)]);
w = w(:,1:size(raw_orig,4));
wi = dcf_fid;
Phi = Vc;
use_gpu=0;

% I don't have the LowRankRecon
params = setupParameters_LowrankMRF2D();
params.block_dim = [6,6];           % locally low-rank patch size
params.block_step = 6;              % overlap between local patches (if equal to patch size above, then patches are non-overlapping)
params.lambdaLLR = 0.03;            % locally low-rank regularization
params.lambdaSpatialTV = 0;%0.003;     % spatial TV regularization
params.lambdaWav = 0.01;               % wavelet regularization
params.betaMethod = 'Dai-Yuan';
params.beta = 0.6;
params.alpha = 0.01;
params.numIter = 20;                    % max number of iterations
params.stopThresh = 1e-3;
params.updateFreq = 0;                  % how often to show intermediate results

%% pattern matching
if doFit
   if (isTRUEFISP && doB0regularization)  
        % load b0 map

        b0mat.name= 'meas_MID01030_FID112465_b0_mapping_rudy_OFFRES.mat'; %phantom
        % b0mat.name= 'meas_MID00343_FID113798_b0_mapping_rudy_200x200_st3_5_OFFRES.mat'; %vol1
        load([RawData.folder b0mat.name]);
        LPF = flip(LPF,3);
        % LPF_2D = LPF(:,:,11); %Rudy 241015 when LPF is only for a slice %vol
        LPF_2D = LPF*2;
        figure, imagesc(LPF_2D), colorbar
        figure, hist(LPF_2D(:),100)

        %pattern match with B0from b0map
        tic
        fprintf('Step 5: pattern match... \n');
    
        %1. b0map binning to dictionary descrete b0 value 
        
        b0values_in_dict = unique(r(:,3));
        % Function to find the closest value in roundVector
        for ix = 1:size(LPF_2D,1)
            for iy = 1:size(LPF_2D,2)
                [v,index] = min(abs(b0values_in_dict - LPF_2D(ix,iy)));
                roundedLPF_2D(ix,iy) = b0values_in_dict(index);
            end
        end
        % imagine(cat(3, LPF_2D, roundedLPF_2D))
    
        %2. b0map masks generation
        for ib0 = 1:length(b0values_in_dict)
            b0mask(:,:,ib0) = roundedLPF_2D == b0values_in_dict(ib0);
        end
        % imagine(b0mask)
        b0mask = b0mask.*repmat(mask_fit,[1,1,size(b0mask,3)]); %fitting region of interest only
    
        %3. iterative mask pm 
        t1map_sum = zeros(size(proxy_img)); t2map_sum = t1map_sum; m0map_sum = t1map_sum;
    
        for ib0 = 1:length(b0values_in_dict)
            fprintf(['b0mask #: ' num2str(ib0) '\n'])
            % select sub-portion of dictionary
            i_dict_b0mask = find(r(:,3)==b0values_in_dict(ib0));
            sub_dict = dict(1:NumFrames,i_dict_b0mask);
            sub_r = r(i_dict_b0mask,:);
    
            [t1map,t2map,~,m0map] = patternmatch(image_combined,b0mask(:,:,ib0),sub_r,0,sub_dict,16);
    
            t1map_sum = t1map_sum + t1map;
            t2map_sum = t2map_sum + t2map;
            m0map_sum = m0map_sum + m0map;
        end
        
        t1map3d(:,:,nn) = t1map_sum;
        t2map3d(:,:,nn) = t2map_sum;
        m0map3d(:,:,nn) = m0map_sum;
        b0map3d(:,:,nn) = roundedLPF_2D;

        % extra for plot and saving single slice
        t1map = t1map_sum;
        t2map = t2map_sum;
        m0map = m0map_sum;
        b0map = roundedLPF_2D.*mask_fit;

        time_pm = toc
   else
        %standard pattern match
        fprintf('Step 5: pattern match... \n');
        
        tic
        
        if singleSliceFit
            fprintf(['Single Slice Fit - Slice n:' num2str(sliceToFit) '\n']);
            % -------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf('LR recon\n');
            image_combined = squeeze(image_combined_3D(:,:,sliceToFit,:));
            [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors

            % -------------------------------------------------------------------------------------------------------------------------------------------------------
            fprintf('LLR recon \n');
            % See these papers for more details...
            % Gastao Cruz, MRM 2019. "Sparsity and locally low rank regularization for MR fingerprinting".
            % Jesse Hamilton, NMR Biomed 2019. "Simultaneous multislice cardiac magnetic resonance fingerprinting using low rank reconstruction".

            DATA = squeeze(raw_orig_Z(:,:,sliceToFit,:)); % [4004 4 38 600]
            DATA = DATA .* permute(repmat(sqrt(squeeze(w)),[1 1 numCoils]),[1 3 2]);
            coilmap = squeeze(dip.coilmap(:,:,sliceToFit,:));

            if use_gpu
                FT = gpuNUFFT([col(kx) col(ky)]',wi(:)/max(wi(:)),readOSFactor,3,8,[N(1)*readOSFactor N(1)*readOSFactor],coilmap);
            else
                FT = NUFFT(kx+1i*ky,wi/max(wi(:)),[0 0],[N(1)*readOSFactor N(1)*readOSFactor]);
            end
    
            % adjoint operator (spiral k-space to image domain)
            Et = @(x)lowrankMRF2D_adjoint(x,coilmap,idproj,Phi,FT,numSpiralArms);
    
            % forward operator (image domain to spiral k-space)
            E = @(x)lowrankMRF2D_forward(x,Phi,FT,idproj,coilmap,[nr numSpiralArms]);
            
            % --------------- sanity check
            % rudy = Et(DATA);
            % itest = squeeze(image_combined_3D(:,:,sliceToFit,:));
            % imagine(cat(4,itest,rudy))
            % --------------- sanity check

            y0 = E(Et(DATA));
            unitv = sum(abs(DATA(:)))/sum(abs(y0(:))); % initial step size
            clear y0
    
            fprintf('computing initial guess\n');
            tic; x0 = Et(DATA); timeAdjoint=toc;
            fprintf('%.2f seconds\n',timeAdjoint);
            scaling = max(abs(x0(:))); % normalize image, so we can use the same regularization parameters for different datasets
            x0 = x0/scaling;
            DATA = DATA/scaling;    
            params.t0 = unitv;                  % initial step size     

            [images_lowrank,t0,dx,obj,update] = nonlinearCGDescent(x0,[],E,Et,DATA,params);
              
            fprintf('matching to dictionary\n');
            [t1mapLLR,t2mapLLR,b0mapLLR,m0mapLLR] = patternmatch(images_lowrank,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors
            
        else %fit whole 3D volume
            
            tic
            for nn=1:size(image_combined_3D,3) %loop over slices
                fprintf(['Slice n:' num2str(nn) '\n']);
                % -------------------------------------------------------------------------------------------------------------------------------------------------------
                fprintf('LR recon \n');
                image_combined = squeeze(image_combined_3D(:,:,nn,:));
                [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors
            
                t1map3d(:,:,nn) = t1map;
                t2map3d(:,:,nn) = t2map;
                m0map3d(:,:,nn) = m0map;
                b0map3d(:,:,nn) = b0map;
                
                % -------------------------------------------------------------------------------------------------------------------------------------------------------
                fprintf('LLR recon \n');            
                
                % See these papers for more details...
                % Gastao Cruz, MRM 2019. "Sparsity and locally low rank regularization for MR fingerprinting".
                % Jesse Hamilton, NMR Biomed 2019. "Simultaneous multislice cardiac magnetic resonance fingerprinting using low rank reconstruction".
       
                DATA = squeeze(raw_orig_Z(:,:,nn,:)); % [4004 4 38 600]
                DATA = DATA .* permute(repmat(sqrt(squeeze(w)),[1 1 numCoils]),[1 3 2]);
                coilmap = squeeze(dip.coilmap(:,:,sliceToFit,:));
    
                if use_gpu
                    FT = gpuNUFFT([col(kx) col(ky)]',wi(:)/max(wi(:)),readOSFactor,3,8,[N(1)*readOSFactor N(1)*readOSFactor],coilmap);
                else
                    FT = NUFFT(kx+1i*ky,wi/max(wi(:)),[0 0],[N(1)*readOSFactor N(1)*readOSFactor]);
                end
        
                % adjoint operator (spiral k-space to image domain)
                Et = @(x)lowrankMRF2D_adjoint(x,coilmap,idproj,Phi,FT,numSpiralArms);
        
                % forward operator (image domain to spiral k-space)
                E = @(x)lowrankMRF2D_forward(x,Phi,FT,idproj,coilmap,[nr numSpiralArms]);
                
                % --------------- sanity check
                % rudy = Et(DATA);
                % itest = squeeze(image_combined_3D(:,:,sliceToFit,:));
                % imagine(cat(4,itest,rudy))
                % --------------- sanity check
    
                y0 = E(Et(DATA));
                unitv = sum(abs(DATA(:)))/sum(abs(y0(:))); % initial step size
                clear y0
        
                fprintf('computing initial guess\n');
                x0 = Et(DATA); timeAdjoint=toc;
                fprintf('%.2f seconds\n',timeAdjoint);
                scaling = max(abs(x0(:))); % normalize image, so we can use the same regularization parameters for different datasets
                x0 = x0/scaling;
                DATA = DATA/scaling;    
                params.t0 = unitv;                  % initial step size     
    
                [images_lowrank,t0,dx,obj,update] = nonlinearCGDescent(x0,[],E,Et,DATA,params);
                
                image_combined_3D_SLLR(:,:,nn,:) =  images_lowrank;
                fprintf('matching to dictionary\n');
                [t1mapLLR,t2mapLLR,b0mapLLR,m0mapLLR] = patternmatch(images_lowrank,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors

                t1map3dLLR(:,:,nn) = t1mapLLR;
                t2map3dLLR(:,:,nn) = t2mapLLR;
                m0map3dLLR(:,:,nn) = m0mapLLR;
                b0map3dLLR(:,:,nn) = b0mapLLR;
            end %for loop over slices
        end %single slice OR 3D volume fit
        %time2clock(toc);
   end
end %doFit

%% 3D recon visual
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 3000]; %1500/3000
climst2 = [0, 1500]; %250/1200
% idx1 = [100:200];
% idx2 = [100:200];
idx1 = s_vert;
idx2 = s_hor;


toSliceViewer = rot90(cat(1,t1map3d(idx1,idx2,:),t1map3dLLR(idx1,idx2,:),t1map3d(idx1,idx2,:)-t1map3dLLR(idx1,idx2,:)));
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

toSliceViewer = rot90(cat(1,t2map3d(idx1,idx2,:),t2map3dLLR(idx1,idx2,:),t2map3d(idx1,idx2,:)-t2map3dLLR(idx1,idx2,:)));
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar

m00 = abs(m0map3d)/max(abs(m0map3d(:)));
m00LR = abs(m0map3dLLR)/max(abs(m0map3dLLR(:)));

toSliceViewer = rot90(cat(1,m00(idx1,idx2,:),m00LR(idx1,idx2,:),m00(idx1,idx2,:)-m00LR(idx1,idx2,:)));
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 0.6], 'ScaleFactors', [1,1,1]), colorbar

%% 3D volume saving
%tempfilename = ['3D_MRF_fit_' date '_' RawData.name(1:13) '_brain_P2_FISP_FA15_nex600_w1500_3Dfit_SVD_LR_SLLR'];
save([RawData.folder], 'image_combined_3D', 't1map3d', 't2map3d', 'm0map3d', 'b0map3d', ...
    'image_combined_3D_SLLR','t1map3dLLR', 't2map3dLLR', 'm0map3dLLR', 'b0map3dLLR');

%% single slice plots
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500]; %3000/1500
climst2 = [0, 300];  %2000/250
% idx1 = [100:200];
% idx2 = [100:200];
idx1 = s_vert;
idx2 = s_hor;
figure, 
%ax1 = subplot(311); imagesc(rot90(cat(1,t1map(idx1,idx2),t1mapLLR(idx1,idx2),t1map(idx1,idx2)-t1mapLLR(idx1,idx2))), climst1), colormap(ax1, T1colormap), colorbar
%ax2 = subplot(312); imagesc(rot90(cat(1,t2map(idx1,idx2),t2mapLLR(idx1,idx2),t2map(idx1,idx2)-t2mapLLR(idx1,idx2))), climst2), colormap(ax2, T2colormap), colorbar

ax1 = subplot(211); imagesc(rot90(cat(1,t1map(idx1,idx2))), climst1), colormap(ax1, T1colormap), colorbar
ax2 = subplot(212); imagesc(rot90(cat(1,t2map(idx1,idx2))), climst2), colormap(ax2, T2colormap), colorbar

m00 = abs(m0map)/max(abs(m0map(:)));
m00LR = abs(m0mapLLR)/max(abs(m0mapLLR(:)));
ax3 = subplot(313); imagesc(rot90(cat(1,m00(idx1,idx2),m00LR(idx1,idx2),m00(idx1,idx2)-m00LR(idx1,idx2))),[0 0.6]), colormap(ax3, gray), colorbar

figure, imagesc(rot90(cat(1,b0map(idx1,idx2),b0mapLLR(idx1,idx2),b0map(idx1,idx2)-b0mapLLR(idx1,idx2))), [-100 100]), colormap(gray), colorbar

%% single slice saving
image_combined = single(image_combined);
tempfilename = ['MRF_fit_' date '_' RawData.name(1:13) '_NIST_s' num2str(sliceToFit) '_SVDdict_Nex600_FA15_FISP_w1500_LR_LLR'];
% tempfilename = ['MRF_fit_' date '_' RawData.name(1:13) '_CAROTID_s10_SVDdict_Nex1000'];
save([RawData.folder tempfilename], 'image_combined', 't1map', 't2map', 'm0map', 'b0map', ...
    'images_lowrank','t1mapLLR', 't2mapLLR', 'm0mapLLR', 'b0mapLLR');











%% denoising traditional

%noise calibration area: pick the indexes where pure noise is present
if singleSliceFit       
    figure, imagesc(abs(squeeze(image_combined_3D(:,:,sliceToFit,1))))
else
   %imagine(image_combined_3D)
end
%% noise ref
idNoiseref = 25:50;
if singleSliceFit     
    noise_area = squeeze(image_combined_3D(idNoiseref,idNoiseref,sliceToFit,:));
    for i = 1:size(noise_area,3) %number of SVD components
        tmp = noise_area(:,:,i);
        noise_ref(:,i) = tmp(:);
        clearvars temp
    end
    std_noise = std(noise_ref);
    avg_std_noise= mean(std_noise);
else
    for s = 1:size(image_combined_3D,3) %slices
        noise_area = squeeze(image_combined_3D(idNoiseref,idNoiseref,s,:));
        for i = 1:size(noise_area,3) %number of SVD components
            tmp = noise_area(:,:,i);
            noise_ref(:,i,s) = tmp(:);
            clearvars temp
        end
        std_noise(s,:) = std(squeeze(noise_ref(:,:,s)));
        avg_std_noise(s)= mean(squeeze(std_noise(s,:)));
    end
end

%% denoising (time consuming - especially 3D)
addpath('C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\matlab_scripts\3D_MRF_FISP_Prostate-main_d240531-main\PCAtensordenoiser');
tic
if singleSliceFit    
    coeff_images = squeeze(image_combined_3D(:,:,sliceToFit,:));

    coeff_images_bm4d = (BM4D(real(squeeze(coeff_images)), avg_std_noise) + 1i*BM4D(imag(squeeze(coeff_images)), avg_std_noise));
    coeff_images_bm3d = (BM3D(real(squeeze(coeff_images)), avg_std_noise) + 1i*BM3D(imag(squeeze(coeff_images)), avg_std_noise));    
    % coeff_images_MPPCA= denoise_recursive_tensor(coeff_images,[5 5], 'sigma', avg_std_noise);
    coeff_images_MPPCA= denoise_recursive_tensor(coeff_images,[5 5]); %Marchenko-Pastour modeling 
    % imagine(cat(4,coeff_images, coeff_images_MPPCA,  coeff_images-coeff_images_MPPCA))
    % imagine(cat(4,coeff_images, coeff_images_bm4d, coeff_images_bm3d, coeff_images-coeff_images_bm4d, coeff_images-coeff_images_bm3d))
else
     for s = 1:size(image_combined_3D,3) %slices
         fprintf(['Denoising slice ' num2str(s) ' - '])
         coeff_images = squeeze(image_combined_3D(:,:,s,:));
         coeff_images_bm4d(:,:,s,:) = (BM4D(real(squeeze(coeff_images)), avg_std_noise(s)) + 1i*BM4D(imag(squeeze(coeff_images)), avg_std_noise(s))); fprintf('BM4D done - ')
         % coeff_images_bm3d(:,:,s,:)  = (BM3D(real(squeeze(coeff_images)), avg_std_noise(s)) + 1i*BM3D(imag(squeeze(coeff_images)), avg_std_noise(s))); fprintf('BM3D done - ')    
         coeff_images_MPPCA(:,:,s,:) = denoise_recursive_tensor(coeff_images,[5 5]); %Marchenko-Pastour modeling 
         fprintf('MPPCA done - \n')    
     end
end
time2clock(toc);

%% fitting
if doFit   
    if singleSliceFit    
        %BM3D
        image_combined = coeff_images_bm3d;
        [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors

        image_combined = single(image_combined);
        tempfilename = ['MRF_fit_' date '_' RawData.name(1:13) '_Brain_HV3_s' num2str(sliceToFit) '_SVDdict_Nex600_FA15_FISP_bm3d'];
        save([RawData.folder tempfilename], 'image_combined', 't1map', 't2map', 'm0map', 'b0map', 'proxy_img')

        %BM4D
        image_combined = coeff_images_bm4d;
        [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors

        image_combined = single(image_combined);
        tempfilename = ['MRF_fit_' date '_' RawData.name(1:13) '_Brain_HV3_s' num2str(sliceToFit) '_SVDdict_Nex600_FA15_FISP_bm4d'];
        save([RawData.folder tempfilename], 'image_combined', 't1map', 't2map', 'm0map', 'b0map', 'proxy_img')

         %MPPCA
        image_combined = coeff_images_MPPCA;
        [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors

        image_combined = single(image_combined);
        tempfilename = ['MRF_fit_' date '_' RawData.name(1:13) '_Brain_HV3_s' num2str(sliceToFit) '_SVDdict_Nex600_FA15_FISP_MPPCA'];
        save([RawData.folder tempfilename], 'image_combined', 't1map', 't2map', 'm0map', 'b0map', 'proxy_img')
    else
        % clearvars t1map3d t2map3d m0map3d b0map3d
        % for s = 1:size(image_combined_3D,3) %slices
        %     %BM3D
        %     % clearvars t1map t2map b0map m0map image_combined
        %     image_combined = squeeze(coeff_images_bm3d(:,:,s,:));
        %     [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors
        %     t1map3d(:,:,s) = t1map;
        %     t2map3d(:,:,s) = t2map;
        %     m0map3d(:,:,s) = m0map;
        %     b0map3d(:,:,s) = b0map;
        % end
        % % % 3D volume saving
        % tempfilename = ['3D_MRF_fit_' date '_' RawData.name(1:13) '_brain_HV3_FISP_FA15_nex600_3Dfit_SVD_bm3d'];
        % save([RawData.folder tempfilename], 't1map3d', 't2map3d', 'm0map3d', 'b0map3d','image_combined_3D')
        % 
         for s = 1:size(image_combined_3D,3) %slices
            %BM4D
            image_combined = squeeze(coeff_images_bm4d(:,:,s,:));
            [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors
            t1map3d(:,:,s) = t1map;
            t2map3d(:,:,s) = t2map;
            m0map3d(:,:,s) = m0map;
            b0map3d(:,:,s) = b0map;
         end
        % 3D volume saving
        tempfilename = ['3D_MRF_fit_' date '_' RawData.name(1:13) '_brain_HV3_FISP_FA15_nex600_3Dfit_SVD_bm4d'];
        save([RawData.folder tempfilename], 't1map3d', 't2map3d', 'm0map3d', 'b0map3d','image_combined_3D')
        % 
         for s = 1:size(image_combined_3D,3) %slices
            %MPPCA
            image_combined = squeeze(coeff_images_MPPCA(:,:,s,:));
            [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors 
            t1map3d(:,:,s) = t1map;
            t2map3d(:,:,s) = t2map;
            m0map3d(:,:,s) = m0map;
            b0map3d(:,:,s) = b0map;
         end
         % 3D volume saving
        tempfilename = ['3D_MRF_fit_' date '_' RawData.name(1:13) '_brain_HV3_FISP_FA15_nex600_3Dfit_SVD_MPPCA'];
        save([RawData.folder tempfilename], 't1map3d', 't2map3d', 'm0map3d', 'b0map3d','image_combined_3D')       
    end
end %doFit


%%
% % % %% pattern matching recursive for multiple times
% % % nIter = 20;
% % % if doFit   
% % %     for iter = 3:nIter
% % %         %standard pattern match
% % %         fprintf(['Step 5: pattern match iter ' num2str(iter) '\n']);
% % % 
% % % 
% % %         if singleSliceFit
% % % 
% % %             image_combined = squeeze(image_combined_3D(:,:,23,:));
% % %             [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors
% % % 
% % %             image_combined = single(image_combined);
% % %             tempfilename = ['MRF_fit_' date '_' RawData.name(1:13) '_Brain_HV3_s23_SVDdict_Nex600_FA15_FISP_iter' num2str(iter)];
% % %             save([RawData.folder tempfilename], 'image_combined', 't1map', 't2map', 'm0map', 'b0map', 'proxy_img')
% % %         end
% % %    end
% % % end %doFit
