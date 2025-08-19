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
doCoilCompression       = 1;
doPrewhitenData         = 1;


TotalSlices = 60; %Slices without partial fourier

CalculateDict           = 0; %if 0 => loads a pregen dictionary
doFit                   =0; %do not run pm fit

doSave                  = 1;
    tempfilename = ['recon_' date '_trueFISP'];
    saveFolder = 'E:\POSTDOC_UoM\08_Project_MRF\';
doRRdebug               = 0;

doNEXtruncation         = 0;
    truncF              = 0.6;
    dict_truncated_is_saved = 1;


%% setup path and irt
%activefilename = matlab.desktop.editor.getActiveFilename;
%[activepath,name,ext] = fileparts(activefilename);
%cd(activepath);
%cd('../MRF_recon_examples/');
addpath('./rr_dictionary/DictSimulation_NoSP');
% addpath('./dictionary/DictSimulation_SP');
% addpath('./data');
addpath('./dependence');
addpath('./Colormap');
addpath('./Functions_Jesus');
addpath('./OpenSiemensRawData/');
addpath('E:/POSTDOC_UoM/05_MATLAB/rr_twix_recon')
%addpath('/data/Matlab_Github/MRF_spiral_3D');
%load('./Colormap/T1cm.mat');
%load('./Colormap/T2cm.mat');
% Path to IRT toolbox
irtPath = 'E:\POSTDOC_UoM\05_MATLAB\fessler-d2310\irt\';
addpath(irtPath); run setup;

addpath('E:\POSTDOC_UoM\05_MATLAB\rr_utilities\');
codePath = 'E:\POSTDOC_UoM\05_MATLAB\Imagine_old\';
addpath(genpath(codePath))

NoiseDataFileName  = []; % measured noise data

%% Open Siemens 3D Raw data
disp('open Siemens Raw data')
% RawData.folder = 'E:\scanner_data\twix_data\241011\'; %NIST
% RawData.folder = 'E:\scanner_data\twix_data\241018\'; %v01
% RawData.folder = 'E:\scanner_data\twix_data\241028\'; %v02
% RawData.folder = 'E:\scanner_data\twix_data\241031\'; %v03 and v04
% RawData.folder = 'E:\scanner_data\twix_data\241119\';
% RawData.folder = 'E:\scanner_data\twix_data\241126\'; %
% RawData.folder = 'E:\scanner_data\twix_data\241205\'; % NIST, VA60, 1mm^3 ISO
% RawData.folder = 'E:\scanner_data\twix_data\241217\'; % BRAIN, VA60, 1mm^3 ISO
% RawData.folder = 'E:\scanner_data\twix_data\250109\';
RawData.folder = 'E:\scanner_data\twix_data\250124\';

if isTRUEFISP
    % RawData.name = 'meas_MID01026_FID112461_rr_MRF_v3_trueFISP_FAbody.dat'; % NIST
    % RawData.name = 'meas_MID00342_FID113797_rr_MRF_v3_trueFISP_st3_5.dat'; %v01 
    % RawData.name = 'meas_MID00052_FID115801_rr_MRF_v3_trueFISP.dat'; %v02
    % RawData.name = 'meas_MID00275_FID116938_rr_MRF_v3_trueFISP.dat'; %v03
    % RawData.name = 'meas_MID00288_FID116951_rr_MRF_v3_trueFISP.dat'; %v04
    % RawData.name = 'meas_MID00051_FID122609_rr_MRF_f300_v2_bSSFP_FA1_5_TR13.dat'; % 
    % RawData.name = 'meas_MID00130_FID00376_rrMRF_f300_1mmISO_FA1_bSSFP.dat'; % NIST, VA60, 1mm^3 ISO, FA1.0 
    % RawData.name = 'meas_MID00027_FID03729_rrMRF_f300_1mmISO_trueFISP_FA10.dat'; % BRAIN, VA60, 1mm^3 ISO, FA1.0 
    % RawData.name = 'meas_MID01291_FID07636_rrMRF_f300_1mmISO_FA1_bSSFP.dat'; % 
    % RawData.name = 'meas_MID01294_FID07639_rrMRF_f300_1mmISO_FA15_bSSFP.dat'; % 
    % RawData.name = 'meas_MID01296_FID07641_rrMRF_f300_1mmISO_FA1_shifted_Nex400_bSSFP.dat'; % 
  
else
    % RawData.name = 'meas_MID01025_FID112460_rr_MRF_v3_FISP.dat'; %NIST
    % RawData.name = 'meas_MID00341_FID113796_rr_MRF_v3_FISP_st3_5.dat'; %v01
    % RawData.name = 'meas_MID00051_FID115800_rr_MRF_v3_FISP.dat'; %v02
    % RawData.name = 'meas_MID00274_FID116937_rr_MRF_v3_FISP.dat'; %v03
    % RawData.name = 'meas_MID00287_FID116950_rr_MRF_v3_FISP.dat'; %v04
    % RawData.name = 'meas_MID00290_FID116953_rr_MRF_v3_FISP_600nex.dat'; 
    % RawData.name = 'meas_MID00166_FID121155_rr_MRF_f300_v2_TR14_FISP_x20.dat'; % 
    % RawData.name = 'meas_MID00162_FID121151_rr_MRF_f300_v2_TR14_FISP.dat'; %
    % RawData.name = 'meas_MID00164_FID121153_rr_MRF_f300_v2_TR14_FISP_x15.dat'; %
    % RawData.name = 'meas_MID00129_FID00375_rrMRF_f300_1mmISO_FA1.dat'; % NIST, VA60, 1mm^3 ISO, FA1.0 
    % RawData.name = 'meas_MID00131_FID00377_rrMRF_f300_1mmISO_FA15.dat'; % NIST, VA60, 1mm^3 ISO, FA1.5
    % RawData.name = 'meas_MID00026_FID03728_rrMRF_f300_1mmISO_FISP_FA15.dat'; % BRAIN, VA60, 1mm^3 ISO, FA1.5
    % RawData.name = 'meas_MID00043_FID03745_rrMRF_f300_ST3_FISP_FA15.dat'; % carotid v1 volume 1, FA1.5
    % RawData.name = 'meas_MID00045_FID03747_rrMRF_f300_ST3_FISP_FA15.dat'; % carotid v1 volume 2, FA1.5
    % RawData.name = 'meas_MID01289_FID07634_rrMRF_f300_1mmISO_FA1.dat'; %
    % RawData.name = 'meas_MID01293_FID07638_rrMRF_f300_1mmISO_FA15.dat'; %
    % RawData.name = 'meas_MID01289_FID07634_rrMRF_f300_1mmISO_FA1.dat'; %
    RawData.name = 'meas_MID00238_FID10664_rr_3D_MRF_pf8_FA15.dat'; %
    
end
RawDataFileName = fullfile(RawData.folder,RawData.name);

[path,name,ext] = fileparts(RawDataFileName);
[raw,noise,log.rawinfo,~] = loadSiemensRawData(RawDataFileName);
size(raw)

%%
if doNEXtruncation
    t = size(raw,4)*truncF;
    raw = raw(:,:,:,1:t);
    log.rawinfo.Nex = t;
    log.rawinfo.nframes = t;
end

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
size(raw) % for nist: 5092 20 1 1000 6
%% Coil Compression
% using SCC coil compression
if doCoilCompression
    fprintf('Coil Compression\n')
    raw = pcaCoilCompress(raw,0.90);
end
raw = permute(raw,[1,5,2,3,4]);
raw = reshape(raw,[size(raw,1),size(raw,2),size(raw,3),size(raw,5)]);
size(raw)

%% Setup acquisition parameter
fprintf('Step 0: Setup Acquisition Parameters\n');

folderTXT = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp';
trtxtfile = [folderTXT '\FISP_TR.txt'];
tetxtfile = [folderTXT '\FISP_TE.txt'];
% fatxtfile = [folderTXT '\FISP_FA_Body.txt'];
fatxtfile = [folderTXT '\FISP_FA_body_x1_5.txt'];
% fatxtfile = [folderTXT '\FISP_FA_body_x2.txt'];
% fatxtfile = [folderTXT '\FISP_FA_Body_shifted.txt'];
phtxtfile = [folderTXT '\FISP_PH.txt'];


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
dictionary_name = 'dict_FISP_TE1_8_TR16_FAbody15_Nex600_SVD';
% dictionary_name = 'dict_FISP_TE1_9_TR14_FAbody_shifted_Nex400_SVD';
if CalculateDict
    fprintf('Step 1: Calculate Dictionary\n');
    log.delay = 3000;
    log.t1series = [10:10:1000 1020:20:2000 2050:50:3000]; %[10:10:2000];
    log.t2series = [2:2:100 105:5:300 310:10:500 520:20:800 850:50:1500 1600:100:2000];
    if isTRUEFISP
        % B0 brain = [-100:5:-45 -40:1:40 45:5:100]
        % B0 nist = [-40:1:40]
        log.offres = [-40:1:40]; 
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
        save([saveFolder '_rr_simulated_dictionaries\' dictionary_name '.mat'],'dict','r','log','-v7.3');
    else
        tic
        [dict,r] = Calculate_MRF_FISP_DictwithDelays(RawDataFileName,log.rawinfo,log.t1series,log.t2series,0,...%[0.6:0.1:1.4] for B1
        tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,2);
        log.timeGen = toc;
        fprintf('Saving FISP dictionary \n')
        save([saveFolder '_rr_simulated_dictionaries\' dictionary_name '.mat'],'dict','r','log','-v7.3');
    end

    fprintf('Compressing dictionary \n')
    [dictSVD,Vc,S] = svd_compress_dictionary(squeeze(dict),1e-2); %Rudy: 1e-2 originally, with offres == 0
    fprintf('Saving SVD compressed dictionary \n')
    save([saveFolder '_rr_simulated_dictionaries\' dictionary_name '_SVD.mat'],'dictSVD','r','log','Vc','-v7.3');

    if lowrank_SVD
        dict = dictSVD;
    end

else
    fprintf('Step 1: Load Dictionary\n');   
    if lowrank_SVD
                        
            if ~dict_truncated_is_saved && doNEXtruncation 
                fprintf('>> non-compressed version  -> to RE-compress...\n');
                % we need to re-compress it! import non-compressed dictionary!
                 dfile = load(['E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\' dictionary_name '.mat'],'dict','r','log');
                 dict = dfile.dict; r = dfile.r; dict_log = dfile.log; 
                 dict = dict(1:log.rawinfo.nframes,:); %truncation
                 fprintf('RE-Compressing dictionary... \n')
                 [dictSVD,Vc,S] = svd_compress_dictionary(squeeze(dict),1e-2); %Rudy: 1e-2 originally, with offres == 0
                 dict=dictSVD;
                 fprintf('Saving SVD compressed dictionary \n')
                 save([saveFolder '_rr_simulated_dictionaries\' dictionary_name '_t' num2str(t) '_compressed.mat'],'dictSVD','r','log','Vc','-v7.3');
            else
                fprintf('>> compressed version... \n');
                dfile = load(['E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\' dictionary_name '.mat'],'dictSVD','r','log','Vc');
                dict = dfile.dictSVD; r = dfile.r; dict_log = dfile.log; Vc = dfile.Vc;
                clearvars dfile
            end
    else    
            dfile = load(['E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\' dictionary_name '.mat'],'dict','r','log');
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
SpiralMeasFileName.folder = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp\';
% SpiralMeasFileName.name = 'spiral_055T_MRF_FOV400mtx400_jesus_241118.mat'; %original from Jesus but acquired by Rudy
% SpiralMeasFileName.name = 'SpiralHT.mat';
SpiralMeasFileName.name = 'spiral_055T_MRF_FOV400mtx400_rudy_241118.mat'; %optimized 400/400 setting from Rudy
% SpiralMeasFileName.name = 'spiral_055T_MRF_FOV300mtx300_rudy_241118.mat'; %optimized 300/300 setting from Rudy
SpMeasFileName = fullfile(SpiralMeasFileName.folder,SpiralMeasFileName.name);
% [kxall,kyall] = GradTrajMeas(SpMeasFileName);

load(SpMeasFileName)

%% Prepare NUFFT

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
dcf_fid = reshape(abs( mri_density_comp_v2(k, 'pipe', 'G', G_fid_dcf.arg.Gnufft)),uplimit,size(kx,2));
adcpad = 20; %MRF data are padded by 20 points

%%
% dip setup
j=1;
for i=1:log.rawinfo.nframes   
    if j>48
        j=1;
    end
        dip.idproj(i) = j;
        dip.k(:,i) = complex(kx(:,j),ky(:,j));  
        dip.w(:,i) = dcf_fid(:,j);
        j = j+1;
end
    
%%
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

%% load raw data

fprintf('Step 3: Gridding spiral data... \n');
size(raw)
%load('./data/MRF_rawdata.mat','raw');
raw = raw(:,:,:,1:log.rawinfo.Nex); % for saving time, we will only process first 1000 images
% Remove adc padding in the raw data
raw = raw(1+adcpad:end-adcpad,:,:,:); % for 1proj
raw = raw(1:uplimit,:,:,:);

disp(size(raw))
disp('Nonzero Slices:')
nsl = size(raw,3) %number of nonzero slices

%%
% Rudy: if lowrank_SVD is off, we can ignore this portion of code!
% Rudy: you cannot run SVD on one partition only here because you are
% considering the full kspace data!
raw_orig = raw;
if lowrank_SVD
    rawsvd = zeros(size(raw,1),48,size(raw,3),size(raw,2),size(dict,1)); %48 spiral rotations
    disp(size(raw))
    
    %rawbackup = raw;
    disp('SVD-compressing raw data...')
    idtxtfile = [folderTXT '\FISP_ID.txt'];
    for nn = 1:size(raw,3),
        disp('SVD-compressing slice number:')
        disp(nn)
        %raw = rawbackup(:,:,nn,:); %JESUS
    

        projID = textread(idtxtfile,'%d\b');


        %Rudy: compressing data following the compressing rule defined
        %for the dictionary (Vc)
        raw1s = compress_kspace_svd(permute(squeeze(raw(:,:,nn,:)),[1,3,2]),Vc,projID,(size(kx,2)));
        raw1s = permute(raw1s,[1,2,4,3]);
        [SpiralReadout,~,NumCoils,NumFrames] = size(raw1s);
        rawsvd(:,:,nn,:,:) = raw1s;
    end
    raw = rawsvd;
else
    [SpiralReadout,NumCoils,NumSlices,NumFrames] = size(raw);
end
disp(size(raw))


%%
raw = transformImageToKspace(raw,3); %Rudy: FFT on kz dimension (from kspace to image space)
size(raw)
raw = flip(raw, 3); % This is to match the B0 and B1 ordering - flip in Kz domain
raw_backup = raw;% JESUS


%%

for nn = 1:size(raw,3),
    % nn = 31;
    disp('Slice number:')
    nn
    
   
    % Perform NUFFT
    image_uncombined = squeeze(single(zeros(log.rawinfo.matrix,log.rawinfo.matrix,NumCoils,NumFrames)));
    tempindex =  zeros(1, NumFrames);

    % f = waitbar(0, 'Gridding');
    for iproj= 1:NumFrames %
        if lowrank_SVD %idproj becomes the SVD dimension (5) and you use only 1 gridding element instead of all 48
            G = G_gridding;
            dcf = dcf_fid;
            for inc=1:NumCoils

                k_temp = squeeze(raw_backup(:,:,nn,inc,iproj));  % fullly sampled data
                image_temp = G'*(dcf(:).*k_temp(:));
                image_uncombined(:,:,inc,iproj) = embed(image_temp,mask);
                % index_progress = (iproj-1)*NumCoils+inc;
                % waitbar(index_progress/(NumFrames*NumCoils), f,...
                %     sprintf('Progress: %d %%', floor(index_progress/(NumFrames*NumCoils)*100)));
            end
        else %idproj here is the num of spiral
            tempindex(:,iproj) = mod(iproj - 1,size(kx,2))+1; %variable that matches Nex to spiral arm
            G = G_gridding{tempindex(:,iproj)};
            dcf = dcf_fid(:,tempindex(:,iproj));
            for inc=1:NumCoils

                k_temp = squeeze(raw_backup(:,inc,nn,iproj));  % single shot
                image_temp = G'*(dcf(:).*k_temp(:));
                image_uncombined(:,:,inc,iproj) = embed(image_temp,mask);
                % index_progress = (iproj-1)*NumCoils+inc;
                %waitbar(index_progress/(NumFrames*NumCoils), f,...
                %JESUS commented
                    %sprintf('Progress: %d %%',
                    %floor(index_progress/(NumFrames*NumCoils)*100)));
                    %%JESUS commented
            end
        end
    end
    
    % close(f)
   
    %
    % % % % % % ksp= transformImageToKspace(transformImageToKspace(image_uncombined,1),2);
    % % % % % % dx = 83.2; %mm
    % % % % % % dy = 105.6; %mm
    % % % % % % [nx,ny,~,~] = size(ksp);
    % % % % % % kx = (-floor(nx/2):ceil(nx/2)-1) / nx; % Normalized kx (frequency domain) range
    % % % % % % ky = (-floor(ny/2):ceil(ny/2)-1) / ny; % Normalized ky (frequency domain) range
    % % % % % % [KX, KY] = meshgrid(ky, kx); % Meshgrid for kx and ky
    % % % % % % % Compute the phase ramp
    % % % % % % phaseRamp = -exp(1i * 2 * pi * (dx * KX + dy * KY));
    % % % % % % % Apply the phase ramp to k-space
    % % % % % % ksp_shifted = ksp .* phaseRamp;
    % % % % % % image_uncombined = transformImageToKspace(transformImageToKspace(ksp_shifted,1),2);
    
    %%
    % CSM estimation
    fprintf('Step 4: Estimate coil sensitivities and perform the coil combination... \n');

    if lowrank_SVD,
        %figure('name','1st coefficient image');
        %    mr_imshow(abs(squeeze(image_uncombined(:,:,:,1))),[],[4 4]);
        %    title('1st coefficient image')

        csm = estimate_csm_walsh(image_uncombined(:,:,:,1));
    else
        avg_image = sum(image_uncombined,4);
        %figure('name','Averaged Image along Time');
        % mr_imshow(abs(avg_image),[],[4 4]);
        csm = estimate_csm_walsh(avg_image);
    end
    %figure('name','Coil Sensitivity Maps');
    %mr_imshow(abs(csm),[],[4,4]);

    % Perform coil combination
    image_combined=single(zeros(log.rawinfo.matrix,log.rawinfo.matrix,NumFrames));

    for iproj = 1:NumFrames
        image_combined(:,:,iproj)=squeeze(sum(conj(csm).*squeeze(image_uncombined(:,:,:,iproj)),3));
    end
    
    % check on PROXY image
    % figure, sliceViewer(abs(image_combined)) %single value image
    proxy_img = sum(image_combined,3);
    dip.coilmap = csm;
    %%


    figure, 
    subplot(121), imagesc(abs(proxy_img)), colormap(gray), colorbar, title('abs'), axis square
    subplot(122), imagesc(angle(proxy_img)), colormap(gray), colorbar, title('phase'), axis square
    

    pixels = [175, 200;
              277, 208;
              204, 223];

    figure, 
    subplot(311), plot(squeeze(abs(image_combined(pixels(1,1), pixels(1,2),:)))),title('fingerprints (ABS)')
    subplot(312), plot(squeeze(abs(image_combined(pixels(2,1), pixels(2,2),:))))
    subplot(313), plot(squeeze(abs(image_combined(pixels(3,1), pixels(3,2),:))))
    xlabel('Nex')


    %% truncate image to optimize the fit to only area of interest
    % phantom 400
    % s_vert = 100:300;
    % s_hor = 130:310;

    % phantom 300
    s_vert = 50:250;
    s_hor = 80:270;

    % brain v01
    % s_vert = 120:300;
    % s_hor = 140:350;

    % brain v02
    % s_vert = 110:300;
    % s_hor = 110:350;

    %brain 300
    % s_vert = 60:240;
    % s_hor = 80:300;

    %carotid 300
    % s_vert = 30:270;
    % s_hor = 50:300;

    a = image_combined(s_vert,s_hor,:);
    imagine(a)
    mask_fit = mask-1;
    mask_fit(s_vert,s_hor) = 1;
    
    ic_out(:,:,:,nn) = a;
% % image_combined = real(image_combined) - j*imag(image_combined);
   %%
    % if isTRUEFISP
    %     image_combined = image_combined*exp(j*pi); %to compensate for error in simulation wehre dictionary is conj
    % end

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
           
            [t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask_fit,r,0,dict(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors
            
            % REAL pattern match -> excecute a complex match but using real signal modulation: it gives slightly different but overall consistent results to complex matching
            % dictnorm = sqrt(sum(dict.*conj(dict)));
            % [t1map,t2map,b0map,m0map] = patternmatch_real(image_combined,mask_fit,r,0,dict,16,dictnorm); %when offres are in blocks 16 helps to avoid RAM outofbound errors
            
            if nn == 1,
                t1map3d = zeros(size(t1map,1),size(t1map,1),size(raw_backup,3));    
                t2map3d = zeros(size(t1map,1),size(t1map,1),size(raw_backup,3));
                m0map3d = zeros(size(m0map,1),size(m0map,1),size(raw_backup,3));
                b0map3d = zeros(size(b0map,1),size(b0map,1),size(raw_backup,3));
            end
            
            t1map3d(:,:,nn) = t1map;
            t2map3d(:,:,nn) = t2map;
            m0map3d(:,:,nn) = m0map;
            b0map3d(:,:,nn) = b0map;
            time_pm = toc
       end
        
       raw = raw_backup;%JESUS
    end %doFit
end % for loop


%% 3D recon visual
t1map3d = rot90((reshape(t1map3d,[size(t1map3d,1) size(t1map3d,1) 1 size(t1map3d,3)])));
t2map3d = rot90((reshape(t2map3d,[size(t2map3d,1) size(t2map3d,1) 1 size(t2map3d,3)])));
m0map3d = rot90((abs(reshape(m0map3d,[size(m0map3d,1) size(m0map3d,1) 1 size(m0map3d,3)]))));

addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500];
climst2 = [0, 300]; %uplimit:1200 for NIST T2 at 0.55T 
sliceSTART = 1;
sliceEND = 20;
% 
% 
toSliceViewer = squeeze(t1map3d);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

toSliceViewer = squeeze(t2map3d);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar

toSliceViewer = squeeze(abs(m0map3d));
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 1], 'ScaleFactors', [1,1,1]), colorbar

toSliceViewer = squeeze(b0map3d);
figure,
sliceViewer(toSliceViewer, "DisplayRange", [-30, 30],'ScaleFactors', [1,1,1]), colorbar

%% 3D volume saving
tempfilename = ['3D_MRF_fit_' date '_' RawData.name(1:13) '_brain_FISP_FA15_nex1000'];
save([RawData.folder tempfilename], 't1map3d', 't2map3d', 'm0map3d', 'b0map3d')

%%
tempfilename = ['3D_SVDimage_' date '_' RawData.name(1:13) '_brain_FISP_FA15_nex600_simPF38_POCS_Rudy'];
save([RawData.folder tempfilename], 'ic_out')

%% single slice plots
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');
climst1 = [0, 3000]; %phantom
% climst1 = [0, 1500];
climst2 = [0, 1200];  %phantom
% climst2 = [0, 300]; 
figure,
ax1 = subplot(131); imagesc(t1map, climst1), colormap(ax1, T1colormap), colorbar
ax2 = subplot(132); imagesc(t2map, climst2), colormap(ax2,T2colormap), colorbar
ax3 = subplot(133); imagesc(abs(m0map)/max(abs(m0map(:)))), colormap(ax3,gray), colorbar

figure, imagesc(b0map), colormap(gray), colorbar
%% single slice saving
image_combined = single(image_combined);
tempfilename = ['MRF_fit_' date '_' RawData.name(1:13) '_NIST_T2L_s14_SVDdict_Nex1000'];
% tempfilename = ['MRF_fit_' date '_' RawData.name(1:13) '_CAROTID_s10_SVDdict_Nex1000'];
save([RawData.folder tempfilename], 'image_combined', 't1map', 't2map', 'm0map', 'b0map', 'proxy_img')



%% save for DIP
% destFolder = 'E:\POSTDOC_UoM\08_Project_MRF\DIP_data\';
% Folder = 'TRUEFISP\';
% subFolderName ='NIST_s11_nex300_FA1_241104';
% mkdir(fullfile(destFolder,Folder,subFolderName))
% 
% coilmap = dip.coilmap;
% save([destFolder Folder subFolderName '/coilmap.mat'], 'coilmap', '-v7.3'); %[Nx, Ny, coils]
% 
% DATA = dip.kspace_data;
% save([destFolder Folder subFolderName '/DATA.mat'],'DATA','-v7.3'); %[Nro, coils, Nex];
% 
% dictCompressed = dict;
% Phi = Vc;
% save([destFolder Folder subFolderName '/dictCompressed.mat'], 'Phi', 'dictCompressed', 'r', '-v7.3'); %Phi: [Nex,EIG], dictCompressed: [EIG,DicSIZE], r: [DicSIZE, DicVAR]
% 
% clearvars k w
% k(:,1,:) = dip.k;
% w(:,1,:) = dip.w;
% save([destFolder Folder subFolderName '/trajectory.mat'], 'k', 'w', '-v7.3'); %k: [Nro,1,Nex], w:[Nro,1,Nex]
% 
% disp('DIP savings DONE!')


%% save for SPIJN and PV-MRF
destFolder = 'E:\POSTDOC_UoM\08_Project_MRF\PV-MRF\';
FolderName ='brain_v01_241018_s11';
mkdir(fullfile(destFolder,FolderName))

xxNorm = sqrt(sum(image_combined.*conj(image_combined),3));
dictnorm = sqrt(sum(dict(1:NumFrames,:).*conj(dict(1:NumFrames,:))));


Ic = zeros(size(image_combined,1),size(image_combined,2),size(image_combined,3)*2);
Ic(:,:,1:size(image_combined,3)) = imag(image_combined./xxNorm);
Ic(:,:,size(image_combined,3)+1:end) = real(image_combined./xxNorm);
a = Ic(s_vert,s_hor,:);
Ic = a;
% T1_list = r(:,1);
% T2_list = r(:,2);


D_norm_factors = dictnorm;
Dc = zeros(size(dict,1)*2, size(dict,2));
Dc(1:size(dict,1),:) = real(dict./dictnorm); %compressed and normalised dictionary
Dc(size(dict,1)+1:end,:) = imag(dict./dictnorm);

% Compr = Vc;

save([destFolder FolderName '/spijn_data.mat'], 'Ic', 'r', 'D_norm_factors', 'Dc'); %k: [Nro,1,Nex], w:[Nro,1,Nex]

%% single slice debug saving
image_combined = single(image_combined);
tempfilename = ['CCdebug_fit_' date '_' RawData.name(1:13) '_NIST_T2L_SVDdict_Nex1000'];
save([RawData.folder tempfilename], 'image_combined', 'proxy_img')