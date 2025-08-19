%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MRF itSENSE recon with BM4D denoising in ADMM or separated (after
% itSENSE)
% v1.0 240807
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rudy Rizzo    (rurizzo@med.umich.edu)
% Jesus Fajardo (jesuserf@umich.edu)
% Yun Jiang     (yunjiang@med.umich.edu)

% Ideal BM4D (simultaneous denoising par-by-par and EIG-by-EIG) is desired.
% BM4D current implementation is a BM3D multichannel where denoising is
% iteratively run through volumetric 3D data [Nx, Ny, Nz, ch] iterating on
% channel ch. This implementation it more efficient than the relative
% BM4D_multichannel presented in the library.

clear;clc;close all;

%% setup reconstruction options
lowrank_SVD             = 1; % 0 :-> 'original' MRF reconstruction
doCoilCompression       = 1;
doPrewhitenData         = 1;
CalculateDict           = 0; %if 0 => loads a pregenerated dictionary
doKTundersmp            = 1;
doSave                  = 1;
    tempfilename = ['recon_' date '_trueFISP_3Tvol'];
    saveFolder = 'E:\POSTDOC_UoM\08_Project_MRF\Jesus_recons\';


%% setup path
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



%% Open Siemens 3D Raw data
fprintf('Step 0: Data load\n');

% phantom data
% RawData.folder = 'E:\POSTDOC_UoM\08_Project_MRF\FreeMax_MRF_Prostate\PhantomData\';
% RawData.name = 'meas_MID00666_FID53511_MRFFISP3DFlexThickness_Sinc.dat';

% prostate data - volunteer 
RawData.folder = 'E:\POSTDOC_UoM\08_Project_MRF\FreeMax_MRF_Prostate\VolunteerData\';
RawData.name = 'meas_MID01730_FID48462_MRFFISP3D_3mm_2pi_Jesus.dat';

% prostate data - patient
% RawData.folder = 'E:\POSTDOC_UoM\08_Project_MRF\FreeMax_MRF_Prostate\PatientData\';
% RawData.name = 'meas_MID00044_FID56331_MRFFISP3D_3mm_2pi_Jesus.dat';

% prostate data - volunteer 3T
% RawData.folder = 'E:\POSTDOC_UoM\08_Project_MRF\FreeMax_MRF_Prostate\InVivoData3T\';
% RawData.name = 'meas_MID02705_FID104846_MRFFISP3D_1mm_2pi_jesus.dat';

tic
RawDataFileName = fullfile(RawData.folder,RawData.name);
[path,name,ext] = fileparts(RawDataFileName);
[raw,noise,log.rawinfo,~] = loadSiemensRawData(RawDataFileName);
time2clock(toc);
size(raw)

%% k-t undersmp
if doKTundersmp
    rawp =permute(raw,[1,3,4,2]); % raw original dim: RO, coils, par, Nex
    [rawUS, smp_mask] = RetroUS(rawp,2,1);
    raw_fs = raw; %store orginal
    raw=permute(rawUS,[1,4,2,3]); %re-permute and keep it for further processing
    clearvars rawp
end

%% Prewhiten data
if doPrewhitenData  
    savedir = RawData.folder;
    fprintf('Step 0.1: Data whitening \n');
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
size(raw)
%% Coil Compression
% using SCC coil compression
if doCoilCompression
    fprintf('Step 0.2: Coil compression \n');
    raw = pcaCoilCompress(raw,0.90);
end
raw = permute(raw,[1,5,2,3,4]);
raw = reshape(raw,[size(raw,1),size(raw,2),size(raw,3),size(raw,5)]);
size(raw)

%% Setup acquisition parameter
fprintf('Step 0.3: Setup Acquisition Parameters\n');

folderTXT = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp';
trtxtfile = [folderTXT '\FISP_TR.txt'];
tetxtfile = [folderTXT '\FISP_TE.txt'];
fatxtfile = [folderTXT '\FISP_FA.txt'];
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

%% Calculate dictionary based on input text file

if CalculateDict
    fprintf('Step 1: Calculate Dictionary \n');   

    log.delay = 3000;
    log.t1series = [10:10:2000 2020:20:3000];% 3050:50:3500 4000:500:5000];
    log.t2series = [2:2:100 105:5:300 310:10:500 520:20:800 850:50:1500 1600:100:2000];

    log.offres = [0];

    EstMemSizeDict = size(log.t1series,2)*size(log.t2series,2)*size(log.offres,2)*log.rawinfo.Nex   *8/1024/1024/1024; %in GB
    fprintf('Estimated dictionary size: %.2f GB \n', EstMemSizeDict)

    tic
    [dict,r] = Calculate_MRF_FISP_DictwithDelays(RawDataFileName,log.rawinfo,log.t1series,log.t2series,0,...%[0.6:0.1:1.4] for B1
    tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,2);
    log.timeGen = toc;
    save([saveFolder '_rr_simulated_dictionaries\' 'rr_dict_NIST_FISP_noSPcorr_noB1.mat'],'dict','r','log','-v7.3');

else
    fprintf('Step 1: Load Dictionary \n');     

    % dfile = load('E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\rr_dict_NIST_FISP_t1t2.mat','dict','r', 'log');
    dfile = load('E:\POSTDOC_UoM\08_Project_MRF\DictionaryFreeMax_Jesus_Aug29_2024\Dict_FISP_noSPcorr.mat','dict','r'); %0.55T dictionary
    % dfile = load('E:\POSTDOC_UoM\08_Project_MRF\FreeMax_MRF_Prostate\InVivoData3T\Dict_FISP_noSPcorr.mat','dict','r'); %3T dictionary
    dict = dfile.dict; r = dfile.r; %dict_log = dfile.log; 
    clearvars dfile
end

%%
if lowrank_SVD
    fprintf('Step 1.1: SVD compressing dictionary\n')
    [dictSVD,Vc,S] = svd_compress_dictionary(squeeze(dict),1e-2); %Rudy: 1e-2 originally, with offres == 0
end

dict = dictSVD;
clearvars dictSVD

%% load spiral trajectory + NUFFT prep
fprintf('Step 2: Prepare NUFFT \n');

%spiral for 0.55T
SpiralMeasFileName.folder = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp\';
SpiralMeasFileName.name = 'SpiralHT.mat';

%spiral for 3T
% SpiralMeasFileName.folder = 'E:\POSTDOC_UoM\08_Project_MRF\FreeMax_MRF_Prostate\InVivoData3T\';
% SpiralMeasFileName.name = 'spiral_3T_MRF_YunJesus_240911.mat';

SpMeasFileName = fullfile(SpiralMeasFileName.folder,SpiralMeasFileName.name);
% [kxall,kyall] = GradTrajMeas(SpMeasFileName);

load(SpMeasFileName)

% Calcuate kmax and the corresponding point
[kmax_measured,uplimit]=max(sqrt(kxall(:,1).^2+kyall(:,1).^2));
% remove additional points
adcpad=20; % the ADC acquires 20 points before and after spiral
% account for the adcpad
kx = kxall(1+adcpad:end-adcpad,:);
ky = kyall(1+adcpad:end-adcpad,:);
kx = kx(1:uplimit,:)/10; %1/cm->1/mm
ky = ky(1:uplimit,:)/10; %1/cm->1/mm

%plot
figure(3);plot(kxall(1:adcpad,1),kyall(1:adcpad,1),'r',...
    kxall(adcpad+1:uplimit,1),kyall(adcpad+1:uplimit,1),'b',...
    kxall(uplimit+1:end,1),kyall(uplimit+1:end,1),'r','LineWidth',2);
axis([-5,5,-5,5]);axis square;
title('1 spiral arm'); xlabel('kx (1/cm)');ylabel('ky (1/cm)');

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
G_fid_dcf = Gmri(k, mask, 'fov', log.rawinfo.fieldofview(1), 'basis', {'dirac'}, 'nufft', nufft_args);
dcf_fid = reshape(abs( mri_density_comp_v2(k, 'pipe', 'G', G_fid_dcf.arg.Gnufft)),uplimit,size(kx,2));

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

%% Gridding

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
    % plot kspace from one coil
    %     figure('name','raw data');imagesc(log(abs(squeeze(raw(:,1,1,:)))));
    %     title('raw data from one coil');
    %     xlabel('Time Points'); ylabel('Data along Spiral');
end

disp(size(raw))

raw = transformImageToKspace(raw,3); %Rudy: FFT on kz dimension
size(raw)
raw = flip(raw, 3); % This is to match the B0 and B1 ordering - flip in Kz domain

%raw = circshift(raw,16,3);%This is just for some patients
raw_backup = raw;% JESUS
image_uncombined = squeeze(single(zeros(log.rawinfo.matrix,log.rawinfo.matrix,NumCoils,nsl,NumFrames)));
% b1index = 0;
for nn = 1:size(raw,3),
    % nn = 11;
    disp('Slice number:')
    nn
    
   
    % Perform NUFFT
    
    tempindex =  zeros(1, NumFrames);

    % f = waitbar(0, 'Gridding');
    for iEx= 1:NumFrames % across excitation OR eigenvalues
        if lowrank_SVD %idproj becomes the SVD dimension (5) and you use only 1 gridding element instead of all 48
            G = G_gridding;
            dcf = dcf_fid;
            for inc=1:NumCoils

                k_temp = squeeze(raw_backup(:,:,nn,inc,iEx));  % fullly sampled data
                image_temp = G'*(dcf(:).*k_temp(:));
                image_uncombined(:,:,inc,nn,iEx) = embed(image_temp,mask);
                index_progress = (iEx-1)*NumCoils+inc;
                % waitbar(index_progress/(NumFrames*NumCoils), f,...
                %     sprintf('Progress: %d %%', floor(index_progress/(NumFrames*NumCoils)*100)));
            end
        else %idproj here is the num of spiral
            tempindex(:,iEx) = mod(iEx - 1,size(kx,2))+1; %variable that matches Nex to spiral arm
            G = G_gridding{tempindex(:,iEx)};
            dcf = dcf_fid(:,tempindex(:,iEx));
            for inc=1:NumCoils

                k_temp = squeeze(raw_backup(:,inc,nn,iEx));  % single shot
                image_temp = G'*(dcf(:).*k_temp(:));
                image_uncombined(:,:,inc,nn,iEx) = embed(image_temp,mask);
                index_progress = (iEx-1)*NumCoils+inc;
                %waitbar(index_progress/(NumFrames*NumCoils), f,...
                %JESUS commented
                    %sprintf('Progress: %d %%',
                    %floor(index_progress/(NumFrames*NumCoils)*100)));
                    %%JESUS commented
            end
        end
    end
end
    % close(f)
%% 3D Coil Sensitivity Map (CSM) estimation
fprintf('Step 4: 3D CSM estimation... \n');
%average through Nex (coils do not depend on the eigenvalue/ex)
image_proxy_coil = permute(sum(image_uncombined,5),[1,2,4,3]);

tic
smoothing = 200;
chunks = size(raw_backup,3);
csm = ismrm_estimate_csm_walsh_3D(single(image_proxy_coil), smoothing, chunks);
toc

image_proxy_combined=squeeze(sum( conj(csm).*image_proxy_coil, 4 ));
% imagine(rot90(image_proxy_combined,1)) %proxy image coil-by-coil visualization

for eig = 1:size(image_uncombined,5)
    image_combined(:,:,:,eig) = squeeze(sum( conj(csm).*permute(squeeze(image_uncombined(:,:,:,:,eig)), [1,2,4,3]), 4 ));
end
% imagine(rot90(image_combined,1)) %coil combined eigenvalue images 





































%% iterative SENSE recon
fprintf('Step 5: itSENSE recon - param definition... \n');
addpath(genpath('E:\POSTDOC_UoM\02_Project_PE_lung_MRI_055T\Jesse-RealTimeCine_DeepImagePrior-main-d2310\Reconstruction_Code\3D\MoCom_recon\code'))

data = permute(raw_backup,[1,2,3,5,4]); % [readout, spirals, slices, EIG, coils]
data = fftshift(fft(ifftshift(data),[],3)); %move back to kspace after SVD compression
dataCS = data/max(abs(data(:)));
csmSENSE = permute(repmat(csm,[1 1 1 1 size(data,4)]), [1,2,3,5,4]);
siz = size(csmSENSE); siz = siz(1:4); %siz(4) = size(data,5);
Ksiz = size(data);

% fully sampled data
masks = ones([Ksiz(2) Ksiz(3)]);

%simulate undersampled Caipirihna-like R=2 
% for p=1:20
%     if mod(p,2)
%         masks(:,p) = repmat([1,0], [1,24]);
%     else
%         masks(:,p) = repmat([0,1], [1,24]);
%     end
% end

ksp_dcfXY = dcf_fid;
spirals = kx+1i*ky;
ksp_loc = spirals;

use_gpu=1;
SENSE_op = itSENSE_SoS_GPU_MRF(masks,csmSENSE,ksp_loc,ksp_dcfXY,siz,Ksiz,use_gpu);
% quick zero-fill test
sense_zf = SENSE_op'*(data.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4) Ksiz(5)])));
%% 2D recon of only 1 slice 

data = permute(raw_backup,[1,2,3,5,4]); % [readout, spirals, slices, EIG, coils]
% % data = fftshift(fft(ifftshift(data),[],3)); %move back to kspace after SVD compression
dataCS = data/max(abs(data(:)));
dataCS2D = squeeze(dataCS(:,:,16,:,:));
csmSENSE = permute(repmat(csm,[1 1 1 1 size(data,4)]), [1,2,3,5,4]);
csmSENSE2D = squeeze(csmSENSE(:,:,16,:,:));
siz = size(csmSENSE2D); siz = siz(1:3); %siz(4) = size(data,5);
Ksiz = size(dataCS2D);

% fully sampled data
masks = ones([Ksiz(2) Ksiz(3)]);

%simulate undersampled Caipirihna-like R=2 
% for p=1:20
%     if mod(p,2)
%         masks(:,p) = repmat([1,0], [1,24]);
%     else
%         masks(:,p) = repmat([0,1], [1,24]);
%     end
% end

ksp_dcfXY = dcf_fid;
spirals = kx+1i*ky;
ksp_loc = spirals;

use_gpu=1;
SENSE_op2D = itSENSE_SoS_GPU_MRF_2D(masks,csmSENSE2D,ksp_loc,ksp_dcfXY,siz,Ksiz,use_gpu);
% quick zero-fill test
sense_zf = SENSE_op2D'*(dataCS2D.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4)])));


%% (1) 3D iterative sense reconstruction 
fprintf('Step 5.1: 3D itSENSE recon ... \n');
% pretty fast: matrix size [400,400,32,5] and 4 sense iteration done in ~6min
tic
itSense.maxit = 10; 
itSense.limit = 3E-3;
[itSENSE_recon,residuals_sense,sense_its] = Conjugate_Gradient((dataCS.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4) Ksiz(5)]))),...
    SENSE_op,itSense.maxit,itSense.limit);
time2clock(toc);

%plot
% imagine(rot90(itSENSE_recon,1))
% imagine(rot90(image_combined,1))

%% (1) 2D iterative sense reconstruction 
fprintf('Step 5.1: 2D itSENSE recon ... \n');
% pretty fast: matrix size [400,400,32,5] and 4 sense iteration done in ~6min
tic
itSense.maxit = 10; 
itSense.limit = 3E-3;
[itSENSE_recon,residuals_sense,sense_its] = Conjugate_Gradient((dataCS2D.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4)]))),...
    SENSE_op2D,itSense.maxit,itSense.limit);
time2clock(toc);

%plot
imagine(rot90(cat(4,itSENSE_recon,squeeze(image_combined(:,:,16,:))),1))

%% (2) 3D itSENSE reconstruction + BM4D denoiser
fprintf('Step 5.2: itSENSE recon + BM4D denoising ... \n');
% it uses CS reconstruction as warm start
% for a dataset sized [400, 400, 32, 5] and 2:2 (ADMM:CG) iterations takes ~26min 
% it requires optimization by testing multiple maxit values (ADMM) and
% lambdas for both

itSenseBM4D.ADMM_maxit = 2; % number of iterations in outer loop of ADMM recon
itSenseBM4D.ADMM_CG_lambda = 0.05;
itSenseBM4D.ADMM_CG_maxit = 2; % number of CG iterations (don't need many since it has warm start)
itSenseBM4D.ADMM_CG_limit = 1E-8; % can also stop iterations by residual, but using a fixed number is easier
itSenseBM4D.bm4d_lambda = 0.01; %lambda denoiser regularization (the higher the smoother)

dataADMM = data/max(abs(data(:)));
itSENSE_recon = itSENSE_recon/ max(abs(itSENSE_recon(:)));

tic
[ADMM_recon,bm3d_its,res_it,ADMM_its] = ADMM_itSENSE_BM4D_MRF_solver((dataADMM.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4) Ksiz(5)]))),...
SENSE_op,itSenseBM4D.ADMM_maxit,itSenseBM4D.ADMM_CG_lambda,itSenseBM4D.ADMM_CG_maxit,itSenseBM4D.ADMM_CG_limit,itSENSE_recon,itSenseBM4D.bm4d_lambda);
time2clock(toc);

% plots
% imagine(cat(4, rot90(ADMM_recon(:,:,:,1),1),rot90(itSENSE_recon(:,:,:,1),1), rot90(image_combined(:,:,:,1),1))) %comparison

%% (2) 2D itSENSE reconstruction + BM3D denoiser
fprintf('Step 5.2: 2D itSENSE recon + BM3D denoising ... \n');
% it uses CS reconstruction as warm start
% for a dataset sized [400, 400, 32, 5] and 2:2 (ADMM:CG) iterations takes ~26min 
% it requires optimization by testing multiple maxit values (ADMM) and
% lambdas for both

itSenseBM3D.ADMM_maxit = 5; % number of iterations in outer loop of ADMM recon
itSenseBM3D.ADMM_CG_lambda = 0.05;
itSenseBM3D.ADMM_CG_maxit = 2; % number of CG iterations (don't need many since it has warm start)
itSenseBM3D.ADMM_CG_limit = 1E-8; % can also stop iterations by residual, but using a fixed number is easier
itSenseBM3D.bm4d_lambda = 0.01; %lambda denoiser regularization (the higher the smoother)

dataADMM = dataCS2D/max(abs(dataCS2D(:)));
itSENSE_recon = itSENSE_recon/ max(abs(itSENSE_recon(:)));

tic
[ADMM_recon,bm3d_its,res_it,ADMM_its] = ADMM_itSENSE_BM4D_MRF_solver_2D((dataADMM.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4)]))),...
SENSE_op2D,itSenseBM3D.ADMM_maxit,itSenseBM3D.ADMM_CG_lambda,itSenseBM3D.ADMM_CG_maxit,itSenseBM3D.ADMM_CG_limit,itSENSE_recon,itSenseBM3D.bm4d_lambda);
time2clock(toc);

% plots
imagine(cat(4, rot90(ADMM_recon,1),rot90(itSENSE_recon,1), rot90(squeeze(image_combined(:,:,16,:)),1))) %comparison

%% alternatively to ADMM denoising, one can run denoising AFTER CS recon - running only (1) without (2)
% % % % % % fprintf('Step 5.3: BM4D denoising after itSENSE ... \n');
% % % % % % % for a dataset sized [400, 400, 32, 5] takes ~19min
% % % % % % t = itSENSE_recon;
% % % % % % norm = max(abs(t(:)));
% % % % % % t = t/norm;
% % % % % % bm4d_lambda = 0.01; %denoising factor
% % % % % % 
% % % % % % tic
% % % % % % tden = zeros(size(t));
% % % % % % for e= 1:size(t,4)
% % % % % %     tden(:,:,:,e) = (BM4D(real(squeeze(t(:,:,:,e))),bm4d_lambda) + 1i*BM4D(imag(squeeze(t(:,:,:,e))),bm4d_lambda));
% % % % % % end
% % % % % % time2clock(toc);
% % % % % % 
% % % % % % %plot
% % % % % % % imagine(rot90(cat(4,t(:,:,16,:),tden(:,:,16,:)),1)) %same partition different EIG

%%
fprintf('Step 6: par-by-par pattern matching ... \n');
% pattern matching par-by-par
t1map3d = zeros(size(image_combined,1),size(image_combined,2),size(image_combined,3));t2map3d=t1map3d;m0map3d=t1map3d;
t1map3dCS =t1map3d;t2map3dCS=t1map3d;m0map3dCS=t1map3d;

% for nn = 6:30
nn=16;
    tic
    %---------------------------------------------------------------------
    % pm for standard reconstructed image
    im = squeeze(image_combined(:,:,nn,:));
    [t1map,t2map,~,m0map] = patternmatch(im,mask,r,0,dict(1:NumFrames,:),32); %when offres are in blocks 16 helps to avoid RAM outofbound errors    
    t1map3d(:,:,nn) = t1map;
    t2map3d(:,:,nn) = t2map;
    m0map3d(:,:,nn) = m0map;

    %---------------------------------------------------------------------
    % pm for compressed sensing reconstructed image
    imCS = squeeze(t(:,:,nn,:));
    % imCS=tden(:,:,nn,:);
    [t1map,t2map,~,m0map] = patternmatch(imCS,mask,r,0,dict(1:NumFrames,:),32); %when offres are in blocks 16 helps to avoid RAM outofbound errors    
    t1map3dCS(:,:,nn) = t1map;
    t2map3dCS(:,:,nn) = t2map;
    m0map3dCS(:,:,nn) = m0map;

    
    time2clock(toc);
% end


%% plotting of maps
fprintf('Step 7: plotting of maps ... \n');

addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');
%%

climst1 = [0, 2000];
climst2 = [0, 600]; %uplimit:1200 for NIST T2 at 0.55T 
sliceSTART = 14;
sliceEND = 14;

% Standard reconstruction PLOTS

toSliceViewer = flip(permute(t1map3d,[2,1,3]),1);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

toSliceViewer = flip(permute(t2map3d,[2,1,3]),1);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar

toSliceViewer = flip(permute(abs(m0map3d),[2,1,3]),1);
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 1], 'ScaleFactors', [1,1,1]), colorbar

% Compressed Sensing PLOTS

toSliceViewer = flip(permute(t1map3dCS,[2,1,3]),1);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

toSliceViewer = flip(permute(t2map3dCS,[2,1,3]),1);
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar

toSliceViewer = flip(permute(abs(m0map3dCS),[2,1,3]),1);
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 1], 'ScaleFactors', [1,1,1]), colorbar


%% saving if needed
fprintf('Step 8: savings ... \n');

if lowrank_SVD
    tempfilename = strcat(tempfilename,'_SVDcompressed');
end
% if SliceProfileCorrection
%     tempfilename = strcat(tempfilename,'_SPcorrected');
% end
if PCA_denoising
    tempfilename = strcat(tempfilename,'_PCAdenoising');
end

if B1correction
    tempfilename = strcat(tempfilename,'_B1Corr');
end
tempfilename = strcat(tempfilename,'_Nex',num2str(log.rawinfo.Nex));
tempfilename = strcat([saveFolder '_recon\' tempfilename 'slice_11_T2array_b0v3_SVD2e-2_26comp'],'.mat');
if doSave
    save(tempfilename,'t1map3d','t2map3d','m0map3d','b0map3d','-v7.3'); %JESUS modified end
    disp('Total .mat file saved');
end
