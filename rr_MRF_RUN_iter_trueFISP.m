%
% MRF recon for 0.55T machine 
%%iterative reconstruction for trueFISP
% Rudy Rizzo
% Based on previous scripts by Jesus Fajardo (jesuserf@umich.edu)
% Based on previous scripts by Yun Jiang (yunjiang@med.umich.edu)

% beautiful conceptually but unfeasible: 30s dicitonary gen + 30s fitting


clear;clc;close all;
tic
%% setup reconstruction options
% Rudy: we need SVD on otherwise dictionary matching phase is too big
lowrank_SVD             = 1; % 0 :-> 'original' MRF reconstruction
% 1 :-> compress dictionary using SVD and perform the
%       reconstuction in low-rank coefficient images.
doCoilCompression       = 1;
doPrewhitenData         = 1;
PCA_denoising           = 0;
OffResonanceDeblurr     = 0; %TODO
B1correction            = 0; %TODO
SliceProfileCorrection  = 0; %TODO
PartialFourier          = 0; %TODO
doSave                  = 1;
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
% cd('./irt'); % setup IRT toolbox
% setup;
% cd('..');

% Path to IRT toolbox
irtPath = 'E:\POSTDOC_UoM\05_MATLAB\fessler-d2310\irt\';
addpath(irtPath); run setup;

addpath('E:\POSTDOC_UoM\05_MATLAB\rr_utilities\');
codePath = 'E:\POSTDOC_UoM\05_MATLAB\Imagine_old\';
addpath(genpath(codePath))


%% Setup data to be constructed
%rawdata file name
% path = './';
% measID = [159279];
% % spiral traj file name
% measID_spiral = [162867];
% 
% if B1correction
%     B1filename = 'meas_MID01029_FID14554_tfl_b1map_FOV400'; %Paste in 2D_MRF_Reconstruction-main
% end
% 
% if OffResonanceDeblurr
%     B0filename = 'meas_MID00512_FID92544_AdjGre';
% end
% 
% RawData = dir(fullfile(path,strcat('*',num2str(measID),'*.dat')));
NoiseDataFileName  = []; % measured noise data
% SpiralMeasFileName = dir(fullfile(path,strcat('*',num2str(measID_spiral),'*.dat')));


%{
% plot the dictionary
figure('name','Dictionary entries');
subplot(211);plot(abs(squeeze(dict(:,:,1000:6000:20000))),'LineWidth',2);
title('Dictionary entries'); xlabel('Time Points');ylabel('Signal Intensity (A.U.)');
%}


%% Open Siemens 3D Raw data
disp('open Siemens Raw data')
RawData.folder = 'E:\scanner_data\twix_data\240725\';


RawData.name = 'meas_MID00394_FID99994_rrMRFv2_hardSpecs_trueFISP.dat';

RawDataFileName = fullfile(RawData.folder,RawData.name);

[path,name,ext] = fileparts(RawDataFileName);
[raw,noise,log.rawinfo,~] = loadSiemensRawData(RawDataFileName);
size(raw)
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

    if ~doCoilCompression
        raw = permute(raw,[1,5,2,3,4]);
        raw = reshape(raw,[size(raw,1),size(raw,2),size(raw,3),size(raw,5)]);
    end
end
size(raw)
%% Coil Compression
% using SCC coil compression
if doCoilCompression
    fprintf('Coil Compression\n')
    raw = pcaCoilCompress(raw,0.90);

    raw = permute(raw,[1,5,2,3,4]);
    raw = reshape(raw,[size(raw,1),size(raw,2),size(raw,3),size(raw,5)]);
end
size(raw)

%% Setup acquisition parameter
fprintf('Step 0: Setup Acquisition Parameters\n');

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

%% load b0 map
b0mat.name= 'meas_MID00398_FID99998_b0_mapping_rudy_OFFRES.mat';
load([RawData.folder b0mat.name]);


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




G_gridding = cell(size(kx,2),1);

for s =1:size(kx,2), %Rudy: build up a G_gridding function for each spiral arm (= n. 48)
    k = double([kx(:,s) ky(:,s)]);
    G_gridding{s} = Gmri(k, mask, 'fov', log.rawinfo.fieldofview(1), 'basis', {'dirac'}, 'nufft', nufft_args);
    %nufft_st{ii,1} = nufft_init(2*pi*tempk,N,J,K,N/2,'minmax:kb');% Fessler's toolbox require the k in radian --YunJiang-01.27.16
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
[SpiralReadout,NumCoils,NumSlices,NumFrames] = size(raw);
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
%raw = circshift(raw,16,3);%This is just for some patients
raw_backup = raw;% JESUS
b1index = 0;
% for nn = 1:size(raw,3),
    nn = 11;
    disp('Slice number:')
    nn
    
   
    % Perform NUFFT
    image_uncombined = squeeze(single(zeros(log.rawinfo.matrix,log.rawinfo.matrix,NumCoils,NumFrames)));
    tempindex =  zeros(1, NumFrames);

    % f = waitbar(0, 'Gridding');
    for iproj= 1:NumFrames %
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
    % close(f)
   
    
    %%
    % CSM estimation
    fprintf('Step 4: Estimate coil sensitivities and perform the coil combination... \n');
        
        avg_image = sum(image_uncombined,4);
        %figure('name','Averaged Image along Time');
        % mr_imshow(abs(avg_image),[],[4 4]);
        csm = estimate_csm_walsh(avg_image);
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
    figure, 
    subplot(311), imagesc(abs(proxy_img)), colormap(gray), colorbar, title('abs'), axis square
    subplot(312), imagesc(real(proxy_img)), colormap(gray), colorbar, title('real'), axis square
    subplot(313), imagesc(imag(proxy_img)), colormap(gray), colorbar, title('imag'), axis square
    
    fprintf('Step 5: Iterative dictionary gen + pattern matching pixel-by-pixel... \n');
    %%
    h = waitbar(0, 'Recon...');
    j=0;

    startx = 95;
    stopx = 300;
    starty = 120;
    stopy = 310;
    dtot = size(image_combined,1)*size(image_uncombined,2);
    dtot = (stopx-startx+1)*(stopy-starty+1);
    for kxi = startx:stopx %1:size(image_combined,1)
        for kyi = starty:stopy %1:size(image_uncombined,2)
            waitbar(j/dtot, h, sprintf('Processing... (%d/%d)', j, dtot));
            b0_est = LPF_2D(kxi,kyi);

            % Calculate dictionary based on input text file
    
            log.delay = 3000;
            log.t1series = [10:10:2000 2020:20:3000];% 3050:50:3500 4000:500:5000];
            log.t2series = [2:2:100 105:5:300 310:10:500 520:20:800 850:50:1500 1600:100:2000];
            % log.offres = [-140:20:140];% [-200:10:50 -45:5:45 50:10:200]; %[-200:10:200];
    
    
            % EstMemSizeDict = size(log.t1series,2)*size(log.t2series,2)*size(log.offres,2)*log.rawinfo.Nex   *8/1024/1024/1024; %in GB
            % fprintf('Estimated dictionary size: %.2f GB \n', EstMemSizeDict)
        
    
    
            [dict,r] = Calculate_MRF_TRUEFISP_DictwithDelays(RawDataFileName,log.rawinfo,log.t1series,log.t2series,b0_est,...%[0.6:0.1:1.4] for B1
            tetxtfile,trtxtfile,fatxtfile,phtxtfile,1,log.delay,0);
            
            image_proxy(1,1,:) = squeeze(image_combined(kxi,kyi,:)); %it has to be NxNxNet 
            [t1map,t2map,b0map,m0map] = patternmatch(image_proxy,1,r,0,dict(1:NumFrames,:),1); % matching 1 pixel to the dictionary   
        
            t1map3d(kxi,kyi) = t1map;
            t2map3d(kxi,kyi) = t2map;
            m0map3d(kxi,kyi) = m0map;
            b0map3d(kxi,kyi) = b0map;
            
            j=j+1;
        end
        j=j+1;
    end
        % raw = raw_backup;%JESUS
    % cd('..');
% end % end

% cd(savedir);
%save("Blurred_Image1.mat","tt_1"); %JESUS
%save("Blurred_Image2.mat","tt_2");
%save("Blurred_Image3.mat","tt_3");
%save("Blurred_Image4.mat","tt_4");

%%

if lowrank_SVD
    tempfilename = strcat(tempfilename,'_SVDcompressed');
end
if SliceProfileCorrection
    tempfilename = strcat(tempfilename,'_SPcorrected');
end
if PCA_denoising
    tempfilename = strcat(tempfilename,'_PCAdenoising');
end
if OffResonanceDeblurr
    tempfilename = strcat(tempfilename,'_FSRdeblurr');
end
if B1correction
    tempfilename = strcat(tempfilename,'_B1Corr');
end
tempfilename = strcat(tempfilename,'_Nex',num2str(log.rawinfo.Nex));
tempfilename = strcat([saveFolder '_recon\' tempfilename 'slice_9_T2array'],'.mat');
if doSave
    save(tempfilename,'t1map3d','t2map3d','m0map3d','b0map3d','-v7.3'); %JESUS modified end
    disp('Total .mat file saved');
end
%dicomwrite(Matrix,'checkMatrix.dcm','MultiframeSingleFile',true) JESUS
%dicom save


%%
% t1map3d = rot90((reshape(t1map3d,[size(t1map3d,1) size(t1map3d,1) 1 size(t1map3d,3)])));
% t2map3d = rot90((reshape(t2map3d,[size(t2map3d,1) size(t2map3d,1) 1 size(t2map3d,3)])));
% m0map3d = rot90((abs(reshape(m0map3d,[size(m0map3d,1) size(m0map3d,1) 1 size(m0map3d,3)]))));

addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 3000];
climst2 = [0, 1500]; %uplimit:1200 for NIST T2 at 0.55T 
sliceSTART = 14;
sliceEND = 14;


toSliceViewer = t1map3d;
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst1, 'ScaleFactors', [1,1,1], 'Colormap', T1colormap), colorbar

toSliceViewer = t2map3d;
figure,
sliceViewer(toSliceViewer, "DisplayRange", climst2, 'ScaleFactors', [1,1,1], 'Colormap', T2colormap), colorbar

toSliceViewer = abs(m0map3d);
figure,
sliceViewer(mat2gray(toSliceViewer), "DisplayRange", [0, 1], 'ScaleFactors', [1,1,1]), colorbar

toSliceViewer = b0map3d;
figure,
sliceViewer(toSliceViewer, "DisplayRange", [-200, 200],'ScaleFactors', [1,1,1]), colorbar
%%