%
% MRF recon for 0.55T machine 
% single debug on PM for specific vials in NIST phantom
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
isTRUEFISP              = 1;
    doB0regularization  = 0;
doCoilCompression       = 1;
doPrewhitenData         = 1;
CalculateDict           = 1; %if 0 => loads a pregen dictionary
doSave                  = 0;
    tempfilename = ['recon_' date '_trueFISP'];
    saveFolder = 'E:\POSTDOC_UoM\08_Project_MRF\';
doRRdebug               = 0;


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
%addpath('/data/Matlab_Github/MRF_spiral_3D');
%load('./Colormap/T1cm.mat');
%load('./Colormap/T2cm.mat');
% Path to IRT toolbox
irtPath = 'E:\POSTDOC_UoM\05_MATLAB\fessler-d2310\irt\';
addpath(irtPath); run setup;

addpath('E:\POSTDOC_UoM\05_MATLAB\rr_utilities\');
codePath = 'E:\POSTDOC_UoM\05_MATLAB\Imagine_old\';
addpath(genpath(codePath))


%% Open Siemens 3D Raw data
disp('open Siemens Raw data')
RawData.folder = 'E:\scanner_data\twix_data\241011\';

if isTRUEFISP
    RawData.name = 'meas_MID01026_FID112461_rr_MRF_v3_trueFISP_FAbody.dat';
else
    RawData.name = 'meas_MID01025_FID112460_rr_MRF_v3_FISP.dat';
end
RawDataFileName = fullfile(RawData.folder,RawData.name);

[path,name,ext] = fileparts(RawDataFileName);
[raw,noise,log.rawinfo,~] = loadSiemensRawData(RawDataFileName);
size(raw)

%% Prewhiten data
if doPrewhitenData 
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
size(raw)
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
% fatxtfile = [folderTXT '\FISP_FA_orig.txt'];
fatxtfile = [folderTXT '\FISP_FA_body.txt'];
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

if CalculateDict,

    log.delay = 3000;
    log.t1series = [10:10:1000 1020:20:2000 2050:50:3000]; %[10:10:2000];
    log.t2series = [2:2:100 105:5:300 310:10:500 520:20:800 850:50:1500 1600:100:2000];
    if isTRUEFISP
        log.offres = [-40:1:40]; %[-30:1:-15, -14.5:0.5:-10.5, -10:0.1:10, 10.5:0.5:14.5, 15:1:30]; %[-20:0.1:20]; %[-15:0.5:15];% [-200:10:50 -45:5:45 50:10:200]; %[-200:10:200];
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
        save([saveFolder '_rr_simulated_dictionaries\' 'rr_dict_NIST_trueFISP_t1t2B0_NIST_v001.mat'],'dict','r','log','-v7.3');
    else
        tic
        [dict,r] = Calculate_MRF_FISP_DictwithDelays(RawDataFileName,log.rawinfo,log.t1series,log.t2series,0,...%[0.6:0.1:1.4] for B1
        tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,2);
        log.timeGen = toc;
        save([saveFolder '_rr_simulated_dictionaries\' 'rr_dict_NIST_FISP_noSPcorr_noB1.mat'],'dict','r','log','-v7.3');
    end

else

    fprintf('Step 1: Load Dictionary\n');

      
    if ~isTRUEFISP %FISP
        dfile = load('E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\rr_dict_NIST_FISP_t1t2.mat','dict','r', 'log');
        dict = dfile.dict; r = dfile.r; dict_log = dfile.log; 
        clearvars dfile
    else           %true FISP
        dfile = load('E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\rr_dict_NIST_trueFISP_t1t2B0v12.mat','dict','r','log');
        dict = dfile.dict; r = dfile.r; dict_log = dfile.log;
        clearvars dfile
    end


end

%%
if lowrank_SVD,
    fprintf('Step 1.1: Compressing dictionary\n')
    [dictSVD,Vc,S] = svd_compress_dictionary(squeeze(dict),1e-2); %Rudy: 1e-2 originally, with offres == 0
end


%%
dict = dictSVD;
clearvars dictSVD

%% load spiral trajectory
fprintf('Step 2: Prepare NUFFT \n');

%load('./data/Spiral_Traj.mat','kxall','kyall');
SpiralMeasFileName.folder = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\rudy_tom_MRF_bssfp\';
SpiralMeasFileName.name = 'SpiralHT.mat';
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
%raw = flip(raw,3);
% rawzeropad = zeros(size(raw,1),size(raw,2),nsl,size(raw,4));
% rawzeropad(:,:,1:size(raw,3),:) = raw;
% raw = rawzeropad;
% size(raw)
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
    % plot kspace from one coil
    %     figure('name','raw data');imagesc(log(abs(squeeze(raw(:,1,1,:)))));
    %     title('raw data from one coil');
    %     xlabel('Time Points'); ylabel('Data along Spiral');
end

disp(size(raw))
%%
raw = transformImageToKspace(raw,3); %Rudy: FFT on kz dimension
size(raw)


%%

% JESUS start
% for loop iterating over N slices JESUS
% tt_1 = zeros(size(raw,3),size(dict,1)); %each of these variables is to save a through-time measurement in certain pixel JESUS
% tt_2 = zeros(size(raw,3),size(dict,1));
% tt_3 = zeros(size(raw,3),size(dict,1));
% tt_4 = zeros(size(raw,3),size(dict,1));

raw = flip(raw, 3); % This is to match the B0 and B1 ordering - flip in Kz domain

%%
%raw = circshift(raw,16,3);%This is just for some patients
raw_backup = raw;% JESUS

nn = 11;
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
            index_progress = (iproj-1)*NumCoils+inc;
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
            index_progress = (iproj-1)*NumCoils+inc;
            %waitbar(index_progress/(NumFrames*NumCoils), f,...
            %JESUS commented
                %sprintf('Progress: %d %%',
                %floor(index_progress/(NumFrames*NumCoils)*100)));
                %%JESUS commented
        end
    end
end
% close(f)


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

%% save for comparison
image_combined = single(image_combined);
tempfilename = ['DataStored_img_rawdatainfo_' date '_' RawData.name(1:13)];
RawData.folder = 'E:\scanner_data\twix_data\241011\';
save([RawData.folder tempfilename], 'image_combined', 'proxy_img', 'log')

%% draw ROIs
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
figure;
imagesc(abs(proxy_img));
axis image;
rois = get_rois;


%% collect fingerprints and avg fingerprint in ROI

% short T2, short T1
[ix,iy] = find(rois{1}==1);
fp1 = reshape(image_combined(ix,iy,:),[size(ix,1)*size(iy,1),NumFrames]);
avg_fp1 = mean(fp1);
std_fp1 = std(abs(fp1));

% long T2, long T1
[ix,iy] = find(rois{2}==1);
fp2 = reshape(image_combined(ix,iy,:),[size(ix,1)*size(iy,1),NumFrames]);
avg_fp2 = mean(fp2);
std_fp2 = std(abs(fp2));


x = 1:NumFrames;
figure, 
subplot(211)
fill([x fliplr(x)], [abs(avg_fp1)+std_fp1 fliplr(abs(avg_fp1)-std_fp1)], [0.7 0.7 0.7], 'EdgeColor', 'none'), hold on
plot(abs(avg_fp1), 'k'), ylabel('T1:   233 , T2:   52')
subplot(212)
fill([x fliplr(x)], [abs(avg_fp2)+std_fp2 fliplr(abs(avg_fp2)-std_fp2)], [0.7 0.7 0.7], 'EdgeColor', 'none'), hold on
plot(abs(avg_fp2), 'k'), ylabel('T1:  1247 , T2:  417')
xlabel('Nex')


%%
% collect fingerprint dictionary close to measured signal
t2_targ =420;
t1_targ =1250;
[~,index] = min(abs((t1_targ+t2_targ) - (r(:,1)+r(:,2))));

[~,index] = min(abs(t2_targ - r(:,2)));
vt2 = r(index,2);
[~,index] = min(abs(t1_targ - r(:,1)));
vt1 = r(index,1);

row_indexes = find( (r(:,1) == vt1) & (r(:,2)== vt2));
v_row_indexes = r(row_indexes,:);

%%
dictnorm = sqrt(sum(dict(1:NumFrames,:).*conj(dict(1:NumFrames,:))));
%%
[t1,t2,b0] = patternmatchDEBUG(avg_fp1,1,r,0,dict,dictnorm,1);




%% multiple ROIS
codePath = 'E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\matlab_roi\';
addpath(genpath(codePath))
% 
figure;
imagesc(abs(proxy_img));
axis image;
rois = get_rois;
%%
avg_fp = zeros(size(rois,2),NumFrames);
for roiIDX = 1: size(rois,2)
    [ix,iy] = find(rois{roiIDX}==1);
    fp = reshape(image_combined(ix,iy,:),[size(ix,1)*size(iy,1),NumFrames]);
    avg_fp(roiIDX,:) = mean(fp);
end

%%
dictnorm = sqrt(sum(dict(1:NumFrames,:).*conj(dict(1:NumFrames,:))));
%%
% avg_fp = avg_fp';

proxyPattern = zeros(1,1,1000);
doPlot=1;

for roiIDX = 1:size(rois,2)
    signal = avg_fp(roiIDX,:);
    % for true fisp, wrong dictionary with a conj(). Either fit to
    % conj(dict) or add conj on signal. and B0 has opposite sign.
    signal = conj(signal);

    % signal = imag(signal) + j*real(signal);

    [t1e,t2e,b0e] = patternmatchDEBUG(signal,1,r,0,dict,dictnorm,1,doPlot);
    t1(roiIDX)=t1e; t2(roiIDX)=t2e; b0(roiIDX) = b0e; 
    % %%%% to use the patternmatch original function
    % proxyPattern(1,1,:) = avg_fp(roiIDX,:);
    % [t1s(roiIDX),t2s(roiIDX),b0s(roiIDX)] = patternmatch(proxyPattern,1,r,0,dict,1);
end

%% plots
NISTdir = ['E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\_NIST_log\'];
NISTname = ['NIST_log01-Oct-2024.mat'];
load([NISTdir NISTname])
NISTplot(:,1) = NIST.T2array.T1.freemax';
NISTplot(:,2) = NIST.T2array.T2.freemax';
%%

figure; 
subplot(321)
plot(flipud(NISTplot(:, 1)), '.', MarkerSize=10), hold on
plot(fliplr(t1), '.', MarkerSize=10);
legend('NIST', 'MRF', 'Location', 'northwest')
title('T1 [ms]'), grid on
subplot(322)
plot(flipud(NISTplot(:, 2)), '.', MarkerSize=10), hold on
plot(fliplr(t2), '.', MarkerSize=10);
legend('NIST', 'MRF', 'Location', 'northwest')
title('T2 [ms]'), grid on

subplot(323)
plot(NISTplot(:, 1), NISTplot(:, 1)), hold on
plot(NISTplot(:, 1), t1, '.', MarkerSize=10)
xlabel('NIST');
ylabel('MRF');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('log(T1)'), grid on
subplot(324)
plot(NISTplot(:, 2), NISTplot(:, 2)), hold on
plot(NISTplot(:, 2), t2, '.', MarkerSize=10);
xlabel('NIST');
ylabel('MRF');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('log(T2)'), grid on

subplot(325)
plot(flipud((t1'-NISTplot(:, 1))./NISTplot(:, 1)), '.', MarkerSize=10);
ylim([-1;1])
title('error T1 [%]'), grid on
subplot(326)
plot(flipud((t2'-NISTplot(:, 2))./NISTplot(:,2)), '.', MarkerSize=10);
ylim([-1;1])
title('error T2 [%]'), grid on

%% specific debug 

[t1,t2,b0] = patternmatchDEBUG(avg_fp(3,:),1,r,0,dict,dictnorm,1);

%%
figure, plot(unwrap(angle(avg_fp(3,:))))


%%

%b0 regularization
b0mat.name= 'meas_MID01030_FID112465_b0_mapping_rudy_OFFRES.mat';
load([RawData.folder b0mat.name]);
figure, imagesc(LPF), colorbar
figure, hist(LPF(:),100)

roundedLPF = LPF;
%%

b0values_in_dict = unique(r(:,3));
% Function to find the closest value in roundVector
for ix = 1:size(LPF,1)
    for iy = 1:size(LPF,2)
        [v,index] = min(abs(b0values_in_dict - LPF(ix,iy)));
        roundedLPF(ix,iy) = b0values_in_dict(index);
    end
end
imagine(cat(3, LPF, roundedLPF))

%%
figure, imagesc(LPF), colorbar 
figure, imagesc(roundedLPF), colorbar
figure, imagesc(LPF-roundedLPF), colorbar

%% sign correction for B0

roundedLPF = -roundedLPF;

%%
roundedLPF = roundedLPF*pi*2;
%%
% find b0 value in every vial
avg_b0 = zeros(size(rois,2),1);
for roiIDX = 1: size(rois,2)
    [ix,iy] = find(rois{roiIDX}==1);
    b0val = reshape(roundedLPF(ix,iy,:),[size(ix,1)*size(iy,1),1]);
    avg_b0(roiIDX) = mean(b0val);
end

%%
% round avg b0 to b0 in dictionary
for roiIDX = 1: size(rois,2)
    [v,index] = min(abs(b0values_in_dict - avg_b0(roiIDX)));
    round_avg_b0(roiIDX) = b0values_in_dict(index);
end
round_avg_b0 = round_avg_b0';
%%
% fit using subspace of dictionary
doPlot=0;
for roiIDX = 1: size(rois,2)
    i_dict_b0mask = find(r(:,3)==round_avg_b0(roiIDX));
    if ~lowrank_SVD
        sub_dict = dict(1:size(dict,1),1,i_dict_b0mask);
    else
        sub_dict = dict(1:size(dict,1),i_dict_b0mask);
    end
    sub_dict_norm = dictnorm(i_dict_b0mask);
    sub_r = r(i_dict_b0mask,:);

    signal = avg_fp(roiIDX,:);
    signal = real(signal) -j*imag(signal);
    [t1(roiIDX),t2(roiIDX),b0(roiIDX)] = patternmatchDEBUG(signal,1,sub_r,0,sub_dict,sub_dict_norm,1,doPlot);
end


%% save the NIST estimated phatom average fingerprints, t1 and t2 from NIST labels and b0 estimated in the area 

meas.avgB0 = avg_b0;
meas.signal = avg_fp;
nist.values =[NIST.T2array.T1.freemax', NIST.T2array.T2.freemax'];


%%
% image_combined = single(image_combined);
tempfilename = ['avgMEAS_' date '_' RawData.name(1:13) '_fp_b0_nist'];
RawData.folder = 'E:\scanner_data\twix_data\241011\';
save([RawData.folder tempfilename], 'meas', 'nist')
