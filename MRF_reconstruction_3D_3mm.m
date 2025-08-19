% Simple script for demonstrating the MRF reconstruction to estimate T1 and
% T2 values.
%
% NUFFT is using the Michigan Imaging Reconstruction Toolbox
% https://web.eecs.umich.edu/~fessler/code/
%
%   Yun Jiang (yunjiang@med.umich.edu)
clear;clc;close all;
tic
%% setup reconstruction options
lowrank_SVD = 1; % 0 :-> 'original' MRF reconstruction
% 1 :-> compress dictionary using SVD and perform the
%       reconstuction in low-rank coefficient images.
PCA_denoising = 0;
CalculateDict = 1;
OffResonanceDeblurr = 0;
B1correction = 0;
SliceProfileCorrection = 0;
linearID = 0;
PartialFourier = 0;
%Nex = 3000;
BinaryGT = 0;
TotalSlices = 60; %Slices without partial fourier

%% setup path and irt
%activefilename = matlab.desktop.editor.getActiveFilename;
%[activepath,name,ext] = fileparts(activefilename);
%cd(activepath);
%cd('../MRF_recon_examples/');
addpath('./dictionary/DictSimulation_NoSP');
addpath('./dictionary/DictSimulation_SP');
addpath('./data');
addpath('./dependence');
addpath('./Colormap');
addpath('./Functions_Jesus');
addpath('./OpenSiemensRawData/');
%addpath('/data/Matlab_Github/MRF_spiral_3D');
%load('./Colormap/T1cm.mat');
%load('./Colormap/T2cm.mat');
cd('./irt'); % setup IRT toolbox
setup;
cd('..');
cd('..');
addpath('./data/');
cd("2D_MRF_Reconstruction-main/");

%% Setup data to be constructed
%rawdata file name
path = './data';
measID = [87654];
% spiral traj file name
measID_spiral = [162867];
B1filename = 'meas_MID05342_FID134887_tfl_b1map_FOV400'; %Paste in 2D_MRF_Reconstruction-main
B0filename = 'meas_MID01105_FID40569_AdjGre';

RawData = dir(fullfile(path,strcat('*',num2str(measID),'*.dat')));
NoiseDataFileName  = []; % measured noise data
SpiralMeasFileName = dir(fullfile(path,strcat('*',num2str(measID_spiral),'*.dat')));


%{
% plot the dictionary
figure('name','Dictionary entries');
subplot(211);plot(abs(squeeze(dict(:,:,1000:6000:20000))),'LineWidth',2);
title('Dictionary entries'); xlabel('Time Points');ylabel('Signal Intensity (A.U.)');
%}

disp('open Siemens Raw data')
%% Open Siemens 3D Raw data
RawDataFileName = fullfile(RawData.folder,RawData.name);
[path,name,ext] = fileparts(RawDataFileName);
[raw,noise,rawinfo,~] = loadSiemensRawData(RawDataFileName);

disp('Load raw noise')
if ~isempty(NoiseDataFileName)
    [noise_raw] = loadSiemensRawData(NoiseDataFileName);
    noise = noise_raw;
end

savedir = fullfile(path,num2str(measID));
addpath(path)
mkdir(savedir);
cd(savedir);

disp('calc noise decorrelation ')
if ~isempty(noise)
    [dmtx] = calculate_noise_decorrelation_mtx(permute(noise,[1,3,4,2]));

else
    dmtx = eye(rawinfo.ncoils);
end

disp('Prewhiten data ')
%% Prewhiten data
if ndims(raw) == 4,
    raw = reshape(raw,[size(raw,1),size(raw,2),size(raw,3),1,size(raw,4)]);
end
raw = permute(raw,[1,3,4,5,2]);
raw = apply_noise_decorrelation_mtx(raw,dmtx);

%% Coil Compression
% using SCC coil compression


if 1  
    fprintf('Coil Compression\n')
    raw = pcaCoilCompress(raw,0.92);
end

%% Setup acquisition parameter
fprintf('Step 0: Setup Acquisition Parameters\n');
%rawinfo.Nex = Nex(ii_Nex)
% rawinfo.Nex = 1000;% time point/ frame
% rawinfo.TR = 12100;% in us
% rawinfo.TE = 2200;% in us
% rawinfo.TI = 20640;% in us
% rawinfo.fieldofview = 400; % in unit of mm
% rawinfo.matrixsize = 400;
% rawinfo.dephase = 4; % 4*pi dephase within each TR

flipangle = importdata('FISP_FA.txt');
tr0 = rawinfo.TR(1)/1000/1000;
tr = importdata('FISP_TR.txt')*1e-6 + tr0;

% plot the parameters
figure('name','Acquisition Parameters');
subplot(211);plot(flipangle(1:rawinfo.Nex,1),'LineWidth',2);
title('Flip Angles'); xlabel('Time Points');ylabel('Flip Angles (degree)');
subplot(212);plot(tr(1:rawinfo.Nex,1)*1000,'LineWidth',2);
title('Repetition Time'); xlabel('Time Points');ylabel('TR (ms)');


%% Calculate dictionary based on input text file

if CalculateDict,

    delay = 3000;
    t1series = [10:10:2000 2020:20:3000 3050:50:3500 4000:500:5000];
    t2series = [2:2:100 105:5:300 310:10:500 520:20:800];% 850:50:1200];% 310:10:500 520:20:800 850:50:1500 1600:100:2000];
    if ~SliceProfileCorrection

        [dict,r] = Calculate_MRF_FISP_DictwithDelays(RawDataFileName,rawinfo,t1series,t2series,0,...%[0.6:0.1:1.4] for B1
            'FISP_TE.txt','FISP_TR.txt','FISP_FA.txt','FISP_PH.txt',1,delay,2);
        save('Dict_FISP_noSPcorr.mat','dict','r','-v7.3');
    else
        %%% --------------------------------------------------------------------------------------
        pulseduration = 2000; % in unit of us
        BandWidthTimeProduct = 8;
        sliceThickness = rawinfo.SliceThickness;% in unit of mm
        baseTR = rawinfo.repetitiontime(1); % in unit of us
        trtxtfile = 'FISP_TR.txt';
        tr = textread(trtxtfile,'%f\b');
        fatxtfile = 'FISP_FA.txt';
        theta = textread(fatxtfile,'%f\b');
        %PHtxtfile = 'FISP_PH.txt';
        %phase = textread(PHtxtfile,'%f\b');
        NumOfFrames = 1000;
        [dict,r] = calcFISP_dict_highres(pulseduration,BandWidthTimeProduct,sliceThickness,baseTR,theta,tr,NumOfFrames,1000);
        dictname = strcat('Dict_FISP_TBP',num2str(BandWidthTimeProduct),'_SliceTH',num2str(sliceThickness),'.mat');
        save(dictname,'dict','r','-v7.3');
        %dict = reshape(dict,[size(dict,1),1,size(dict,2)]);
        %%% --------------------------------------------------------------------------------------
    end

else

    fprintf('Step 1: Load Dicionary\n');
    if ~SliceProfileCorrection
        if B1correction
            cd('..');
            cd('..');
            load('Dict_FISP_noSPcorr_B1corr_0.6_1.4.mat','dict','r'); %Generate dict with script GenerateDict_B1_noSP.m
        else
            load('Dict_FISP_noSPcorr.mat','dict','r');
        end
    else
        BandWidthTimeProduct = 8;
        sliceThickness = rawinfo.SliceThickness;% in unit of mm
        dictname = strcat('Dict_FISP_TBP',num2str(BandWidthTimeProduct),'_SliceTH',num2str(sliceThickness),'.mat');
        load(dictname,'dict','r');

    end
end
%disp('dictionary saved!')
%   pause(9999999999999999999999999999999)
% dict = dict(1:rawinfo.Nex,:); JESUS commented
%{
% plot the dictionary
figure('name','Dictionary entries');
subplot(211);plot(abs(squeeze(dict(:,:,1000:6000:20000))),'LineWidth',2);
title('Dictionary entries'); xlabel('Time Points');ylabel('Signal Intensity (A.U.)');
%}

if lowrank_SVD,
    fprintf('Compressing dictionary\n')
    [dict,Vc,S] = svd_compress_dictionary(squeeze(dict),1e-2);

    subplot(212);plot(abs(squeeze(dict(:,1000:6000:end))),'LineWidth',2);
    title('Compressed Dictionary entries'); xlabel('Time Points');ylabel('Signal Intensity (A.U.)');
    plot(diag(S(1:20,1:20)),'o-','linewidth',2);xlabel('First 20 Singular Values');


end

%% load spiral trajectory
fprintf('Step 2: Prepare NUFFT \n');

%load('./data/Spiral_Traj.mat','kxall','kyall');

SpiralMeasFileName = fullfile(SpiralMeasFileName.folder,SpiralMeasFileName.name);
[kxall,kyall] = GradTrajMeas(SpiralMeasFileName);

%{
plot kspace
figure(2);
subplot(1,2,1);plot(kxall(:,1),kyall(:,1),'LineWidth',2);axis([-5,5,-5,5]);axis square;
title('1 spiral arm'); xlabel('kx (1/cm)');ylabel('ky (1/cm)');
subplot(1,2,2);
v = VideoWriter('spiral.avi');
open(v);
for ii = 1:size(kxall,2)
    plot(kxall(:,ii),kyall(:,ii),'LineWidth',2);xlim([-5,5]);ylim([-5,5]);axis square;
    title(strcat(num2str(ii),'/',num2str(size(kxall,2)), 'Spiral Trajectory'));
    xlabel('kx (1/cm)');ylabel('ky (1/cm)');
    pause(0.1);
    drawnow;
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
%}



%% Prepare NUFFT
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
resolution = rawinfo.fieldofview(1,1)/rawinfo.matrix; % spatial resolution
kmax = 1/(2*resolution);
resolution_measured = 1/(2*(kmax_measured/10));
% normalized spiral trajector to [-0.5 0.5]
kx = kx./(2*kmax_measured/10);
ky = ky./(2*kmax_measured/10);
k = [kx(:) ky(:)];

mask=true([rawinfo.matrix rawinfo.matrix]);
N=size(mask);

nufft_args = {N, [5 5], 2*N, N/2, 'table', 2^12, 'minmax:kb'};

% Calculate DCF
G_fid_dcf = Gmri(k, mask, 'fov', rawinfo.fieldofview(1), 'basis', {'dirac'}, 'nufft', nufft_args);
dcf_fid = reshape(abs( mri_density_comp_v2(k, 'pipe', 'G', G_fid_dcf.arg.Gnufft)),uplimit,size(kx,2));

% Calculate nufft
if lowrank_SVD,
    k= double([kx(:) ky(:)]);
    G_gridding = Gmri(k, mask, 'fov', rawinfo.fieldofview(1), 'basis', {'dirac'}, 'nufft', nufft_args);
else
    G_gridding = cell(size(kx,2),1);

    for s =1:size(kx,2),
        k = double([kx(:,s) ky(:,s)]);
        G_gridding{s} = Gmri(k, mask, 'fov', rawinfo.fieldofview(1), 'basis', {'dirac'}, 'nufft', nufft_args);
        %nufft_st{ii,1} = nufft_init(2*pi*tempk,N,J,K,N/2,'minmax:kb');% Fessler's toolbox require the k in radian --YunJiang-01.27.16
    end
end

%% load raw data

fprintf('Step 3: Gridding spiral data... \n');

%load('./data/MRF_rawdata.mat','raw');
raw = permute(raw,[1,5,2,3,4]);
raw = reshape(raw,[size(raw,1),size(raw,2),size(raw,3),size(raw,5)]);
raw = raw(:,:,:,1:rawinfo.Nex); % for saving time, we will only process first 1000 images
% Remove adc padding in the raw data
raw = raw(1+adcpad:end-adcpad,:,:,:); % for 1proj
raw = raw(1:uplimit,:,:,:);


disp(size(raw))
disp('Nonzero Slices:')
nsl = size(raw,3) %number of nonzero slices
%raw = flip(raw,3);
rawzeropad = zeros(size(raw,1),size(raw,2),TotalSlices,size(raw,4));
rawzeropad(:,:,1:size(raw,3),:) = raw;
raw = rawzeropad;

rawsvd = zeros(size(raw,1),48,size(raw,3),size(raw,2),size(dict,1)); %48 spiral rotations
disp(size(raw))

%rawbackup = raw;
disp('SVD-compressing raw data...')
for nn = 1:size(raw,3),
    disp('SVD-compressing slice number:')
    disp(nn)
    %raw = rawbackup(:,:,nn,:); %JESUS

    if linearID,
        projID = textread('FISP_ID.txt','%d\b');

    else
        projID = textread('FISP_ID.txt','%d\b');
    end

    if lowrank_SVD
        raw1s = compress_kspace_svd(permute(squeeze(raw(:,:,nn,:)),[1,3,2]),Vc,projID,(size(kx,2)));
        raw1s = permute(raw1s,[1,2,4,3]);
        [SpiralReadout,~,NumCoils,NumFrames] = size(raw1s);
        rawsvd(:,:,nn,:,:) = raw1s;
    else
        [SpiralReadout,NumCoils,NumSlices,NumFrames] = size(raw);
        % plot kspace from one coil
        %     figure('name','raw data');imagesc(log(abs(squeeze(raw(:,1,1,:)))));
        %     title('raw data from one coil');
        %     xlabel('Time Points'); ylabel('Data along Spiral');
    end
end
raw = rawsvd;
disp(size(raw))

if PartialFourier,
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Applying Partial Fourier
    disp('Using partial Fourier Fraction:')
    raw(:,:,nsl+1:end,:,:) = 0;
    symkspace = raw;
    symkspace(:,:,1:size(raw,3)-nsl,:,:) = 0; % center-symmetric k-space
    sym_im = transformKspaceToImage(symkspace,3);
    im_ph = angle(sym_im);



    %plot(abs(squeeze(symkspace(1,1,:,1))))


    %% Apply POCS method to retrieve phase
    disp('Applying POCS...')
    threshold_pocs = 0.001;
    %Zero padding for initial guess
    im_init = transformKspaceToImage(raw,3);
    %plot(abs(squeeze(im_init(1,1,:,1))))


    %take only magnitude term & apply phase term
    im_init = abs(im_init).*exp(1i*im_ph);
    disp((size(raw)))
    disp((round(size(raw,3)-nsl)))
    %FFT
    tmp_k = transformImageToKspace(im_init,3);
    diff_im = threshold_pocs +1;
    while (abs(diff_im) > threshold_pocs)
        bottomkspace = transformImageToKspace(im_init,3);
        tmp_k(:,:,round(size(raw,3)-nsl):end,:,:) = bottomkspace(:,:,round(size(raw,3)-nsl):end,:,:);
        %plot(abs(squeeze(tmp_k(1,1,:,1))))
        tmp_im = transformKspaceToImage(tmp_k,3); %inverse DFT
        % apply phase of the k-space center
        tmp_im = abs(tmp_im).*exp(1i*im_ph);
        tmp_k = transformImageToKspace(tmp_im,3);
        % Compare the reconstructed image
        diff_im = abs(tmp_im - im_init);
        diff_im = sum(diff_im(:).^2);
        fprintf('POCS Difference: %f\n',diff_im);

        im_init = tmp_im;
    end
    im_pocs = tmp_im;
    raw = im_pocs; %This is fourier-transformed just along z, that's why is in 'image' space
    raw = circshift(raw,-1,3);
    raw = flip(raw, 3);
else
    raw = transformImageToKspace(raw,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cd('data/');
if OffResonanceDeblurr
    LPF = Load_offresonance_multislice(B0filename);
    disp('Loading off-resonance map...')
end

if B1correction
    LPF_B1 = Load_B1map_multislice(B1filename);
    LPF_B1 = rot90(LPF_B1,1);
    LPF_B1 = flip(LPF_B1,3);
    disp('loading B1 map...')
end

% JESUS start
% for loop iterating over N slices JESUS
tt_1 = zeros(size(raw,3),size(dict,1)); %each of these variables is to save a through-time measurement in certain pixel JESUS
tt_2 = zeros(size(raw,3),size(dict,1));
tt_3 = zeros(size(raw,3),size(dict,1));
tt_4 = zeros(size(raw,3),size(dict,1));

raw = flip(raw, 3); % This is to match the B0 and B1 ordering
raw_backup = raw;% JESUS

for nn = 1:size(raw,3),
    disp('Slice number:')
    nn
    
    %JESUS
    if OffResonanceDeblurr
        %cd('..');
        %cd('..');
        %disp(size(raw))
        disp('Loading off-resonance map and performing deblurring...')
        L = 50; %FREQUENCIES

        
        %load('Raw_spiral_data.mat','raw'); %(raw spiral data)
        %InterpolatedField = rot90(InterpolatedField,3)
        %LPF = rot90(LPF,3);
        %LPF = flipud(LPF);
        % interpolatedDicom = rot90(InterpolatedDicom,3); %JESUS 
        %mr_imshow(abs((LPF(:,:,nn))));
        
        [demod_data, deltaw] = B0_CPR_deblurr_1(raw, LPF(:,:,nn), L);
        %save('c_i.mat','coeff_table')
        
    end

    if OffResonanceDeblurr
        %% Perform NUFFT
        image_uncombined = squeeze(single(zeros(rawinfo.matrix,rawinfo.matrix,NumCoils,NumFrames,1))); %(single(zeros(rawinfo.matrix,rawinfo.matrix,NumCoils,NumFrames,L+1)));
        %tempindex =  zeros(1, NumFrames);
        tic
        %f = waitbar(0, 'Gridding');
        parfor ll = 1:L+1
            tempindex =  zeros(1, NumFrames);
            for iproj= 1:NumFrames
                if lowrank_SVD
                    G = G_gridding;
                    dcf = dcf_fid;
                    for inc=1:NumCoils

                        k_temp = squeeze(demod_data(:,:,inc,iproj,ll));  % fully sampled data
                        image_temp = G'*(dcf(:).*k_temp(:));
                        image_uncombined(:,:,inc,iproj,ll) = embed(image_temp,mask);    
                        %index_progress = (iproj-1)*NumCoils+inc;
                        %waitbar(index_progress/(NumFrames*NumCoils*(L+1)), f,...
                            %sprintf('Progress: %d %%', floor(index_progress/(NumFrames*NumCoils*(L+1))*100)));
                    end
                else
                    tempindex(:,iproj) = mod(iproj - 1,size(kx,2))+1;
                    G = G_gridding{tempindex(:,iproj)};
                    dcf = dcf_fid(:,tempindex(:,iproj));
                    for inc=1:NumCoils

                        k_temp = squeeze(demod_data(:,:,inc,iproj,ll));  % single shot
                        image_temp = G'*(dcf(:).*k_temp(:));
                        image_uncombined(:,:,inc,iproj,ll) = embed(image_temp,mask);
                        %index_progress = (iproj-1)*NumCoils+inc;
                        %waitbar(index_progress/(NumFrames*NumCoils*(L+1)), f,...
                            %sprintf('Progress: %d %%', floor(index_progress/(NumFrames*NumCoils*(L+1))*100)));
                    end
                end
            end
        end
        ParNUFFT = toc
        %disp(image_uncombined)
        %save('Demod_Imgs.mat','image_uncombined') % Demodulated MFI images
        %pause(9999999999999999999)
        
        %image_uncombined = B0_CPR_deblurr_2(image_uncombined, coeff_table, InterpolatedDicom, L, num_of_interp);
        
    else
        %% Perform NUFFT
        image_uncombined = squeeze(single(zeros(rawinfo.matrix,rawinfo.matrix,NumCoils,NumFrames)));
        tempindex =  zeros(1, NumFrames);

        f = waitbar(0, 'Gridding');
        for iproj= 1:NumFrames
            if lowrank_SVD
                G = G_gridding;
                dcf = dcf_fid;
                for inc=1:NumCoils

                    k_temp = squeeze(raw_backup(:,:,nn,inc,iproj));  % fullly sampled data
                    image_temp = G'*(dcf(:).*k_temp(:));
                    image_uncombined(:,:,inc,iproj) = embed(image_temp,mask);
                    index_progress = (iproj-1)*NumCoils+inc;
                    waitbar(index_progress/(NumFrames*NumCoils), f,...
                        sprintf('Progress: %d %%', floor(index_progress/(NumFrames*NumCoils)*100)));
                end
            else
                tempindex(:,iproj) = mod(iproj - 1,size(kx,2))+1;
                G = G_gridding{tempindex(:,iproj)};
                dcf = dcf_fid(:,tempindex(:,iproj));
                for inc=1:NumCoils

                    k_temp = squeeze(raw_backup(:,:,nn,inc,iproj));  % single shot
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
        close(f)
    end
    
    if OffResonanceDeblurr
            %% CSM estimation
            fprintf('Step 4: Estimate coil sensitivities and perform the coil combination... \n');
            template_for_csm = squeeze(image_uncombined(:,:,:,:,1));
            if lowrank_SVD,
                %figure('name','1st coefficient image');
                %    mr_imshow(abs(squeeze(image_uncombined(:,:,:,1))),[],[4 4]);
                %    title('1st coefficient image')
                
                csm = estimate_csm_walsh(template_for_csm(:,:,:,1));
            else
                avg_image = sum(template_for_csm,4);
                %figure('name','Averaged Image along Time');
                % mr_imshow(abs(avg_image),[],[4 4]);
                csm = estimate_csm_walsh(avg_image);
            end
            %figure('name','Coil Sensitivity Maps');
            %mr_imshow(abs(csm),[],[4,4]);
            image_combined=single(zeros(rawinfo.matrix,rawinfo.matrix,NumFrames,L+1));
            for ll = 1:L+1
                %% Perform coil combination
                for iproj = 1:NumFrames
                    image_combined(:,:,iproj,ll)=squeeze(sum(conj(csm).*squeeze(image_uncombined(:,:,:,iproj,ll)),3));
                end
            end
    %for nn = 1:L+1 %JESUS start
    %    if lowrank_SVD,
    %        mr_imshow(abs(squeeze(image_combined(:,:,:,nn))));%abs(image_combined));%
    %    end
    %end %JESUS end
    image_combined_deblurred = B0_CPR_deblurr_2(image_combined, deltaw, LPF(:,:,nn));
    image_combined = zeros(size(image_combined_deblurred));
    image_combined = image_combined_deblurred;
    else
        %% CSM estimation
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

        %% Perform coil combination
        image_combined=single(zeros(rawinfo.matrix,rawinfo.matrix,NumFrames));

        for iproj = 1:NumFrames
            image_combined(:,:,iproj)=squeeze(sum(conj(csm).*squeeze(image_uncombined(:,:,:,iproj)),3));
        end
    end
    
    
    %imagenes a x,y,tp,1
    %
    % figure('name','Image of each time point');
    % set(gcf,'Position',[100 100 1500 600])
    % v = VideoWriter('MRFdata.avi');
    % open(v);
    % for ii = 1:NumFrames
    %
    %     subplot(121);
    %     imagesc(abs(image_combined(:,:,ii)));axis image;colormap gray;
    %     title(strcat(num2str(ii),'/',num2str(NumFrames), ' time point'));
    %
    %     subplot(122);
    %     imagesc(abs(mean(image_combined(:,:,1:ii),3)));axis image;colormap gray;
    %     title('Averaged Image along Time');
    %     frame = getframe(gcf);
    %     writeVideo(v,frame);
    %
    %     drawnow;
    % end
    % close(v);
    %disp(size(image_combined))
    %for nn = 1:L+1 %JESUS start
    %    if lowrank_SVD,
    %        mr_imshow(abs(squeeze(image_combined(:,:,:,nn))));%abs(image_combined));%
    %    end
    %end            %JESUS end
    
    %Plot timepoints images
    %mr_imshow(abs((image_combined(:,:,950))));
    %tt_1(nn,:) = image_combined(220,151,:);
    %tt_2(nn,:) = image_combined(96,123,:);
    %tt_3(nn,:) = image_combined(124,159,:);
    %tt_4(nn,:) = image_combined(159,124,:);
    
    
    %JESUS
    %save("Blurred_Image.mat","image_combined"); %This is to make the mask on
    %the LoadOffresonance function (probably the name is misleading ~~')
    %% Calculate the Inner product between MRF images and the dictionary to get T1 and T2 maps
    fprintf('Step 5: pattern match... \n');
    
    if B1correction
        disp('B1_index:')
        if mod(nn,2) == 1 % odd iteration
            b1index = (nn+1)/2 % increment index by 1
        else % even iteration
            b1index = nn/2 % leave index unchanged
        end
        b1val = zeros(size(LPF_B1(:,:,b1index),1),size(LPF_B1(:,:,b1index),2)); %B1 value to look for later on the dict matching (:,:,nn)

    if nn == 1,
        t1map3d = zeros(size(LPF_B1(:,:,b1index),1),size(LPF_B1(:,:,b1index),1),size(raw_backup,3));    
        t2map3d = zeros(size(LPF_B1(:,:,b1index),1),size(LPF_B1(:,:,b1index),1),size(raw_backup,3));
        m0map3d = zeros(size(LPF_B1(:,:,b1index),1),size(LPF_B1(:,:,b1index),1),size(raw_backup,3));
    end
    
    tic
    internalindex = size(LPF_B1(:,:,b1index),2);
    %poolsize = 20;
    %parpool(poolsize);
    parfor k = 1:size(LPF_B1(:,:,b1index),1)
        disp(k/size(LPF_B1(:,:,b1index),1)*100)
        for l = 1:internalindex
            if abs(LPF_B1(k,l,b1index) > 0.0001),
            % Calculate the absolute difference between the value and array elements
            diff_arr = abs(LPF_B1(k,l,b1index) - r(:,3));
            % Find the index or indexes of the minimum difference
            min_idx = find(diff_arr == min(diff_arr));
            % Store the row and column indices in the output matrix
            b1val(k,l) = r(min_idx(1),3);
            matchingIndexes = find(abs(r(:,3) - b1val(k,l)) < 0.0001);
            tempr =  r(matchingIndexes,:);
            tempdict = dict(:,matchingIndexes);
            %[t1val,t2val,~,m0val] = patternmatch(image_combined(k,l,:),mask,tempr,0,tempdict(1:NumFrames,:),4)
            
            trange = 1:size(image_combined(k,l,:),3);
            dictnorm = sqrt(sum(tempdict(trange,:).*conj(tempdict(trange,:))));
            xx = squeeze(image_combined(k,l,:)).';
            xxNorm = sqrt(sum(xx.*conj(xx),2));
            normAll = xxNorm*dictnorm; xxNorm = [];% clear xxNorm
            innerProduct = conj(xx)*tempdict(trange,:)./normAll; normAll = [];% clear normAll
            [value indexm] = max(abs(innerProduct),[],2); innerProduct = [];% clear innerProduct
            dictCol = tempdict(trange,indexm);
            [t1t2r dfc ] = ind2sub([size(tempr,1), 0],indexm);
            t1map3d(k,l,nn) = tempr(t1t2r(:),1);
            t2map3d(k,l,nn) = tempr(t1t2r(:),2);
            m0map3d(k,l,nn) = 0;
            %t1map3d(k,l,nn) = t1val;
            %t2map3d(k,l,nn) = t2val;
            
            end
        end
    end
    
    MatchTime = toc
    
    else 
    [t1map,t2map,~,m0map] = patternmatch(image_combined,mask,r,0,dict(1:NumFrames,:),4);
    if nn == 1,
        t1map3d = zeros(size(t1map,1),size(t1map,1),size(raw_backup,3));    
        t2map3d = zeros(size(t1map,1),size(t1map,1),size(raw_backup,3));
        m0map3d = zeros(size(m0map,1),size(m0map,1),size(raw_backup,3));
    end
    
    t1map3d(:,:,nn) = t1map;
    t2map3d(:,:,nn) = t2map;
    m0map3d(:,:,nn) = m0map;
    end
    %%mr_imshow(abs((image_combined(:,:,1))));
    %mr_imshow(abs((LPF_B1(:,:,b1index))));
    %if exist('Vc','var')
    %    uncompressed_image = decompress_images_svd(image_combined,Vc);
    %end

    %normeddict = dict./repmat(dictnorm,size(dict,1),1);


    %disp(size(t1map)) 
    %disp(size(t2map))
    %subplot()
    %imagesc(t2map);
    %title('T2')
    %colorbar
    %JESUS commented start
    %figure('name','T1 T2 and Proton Density');
    %set(gcf,'Position',[100 100 1500 600])
    %ax1 = subplot(121);imagesc(t1map(:,end:-1:1).');axis image;colorbar;colormap(ax1,T1colormap);title('T1');caxis([0 3000])
    %ax2 = subplot(122);imagesc(t2map(:,end:-1:1).');axis image;colorbar;colormap(ax2,T1colormap);title('T2');caxis([0 300])
    %ax3 = subplot(133);imagesc(mat2gray(abs(m0map)));axis
    %image;colorbar;colormap(ax3,gray);title('Proton Density'); JESUS
    %commented end

    if BinaryGT,

        dict_filename		= strcat('dictc_TBP',num2str(BandWidthTimeProduct),'TH',num2str(sliceThickness),'Nex',num2str(rawinfo.Nex),'.bin');
        dictnorm_filename	= strcat('dictcnorm_TBP',num2str(BandWidthTimeProduct),'TH',num2str(sliceThickness),'Nex',num2str(rawinfo.Nex),'.bin');
        V_filename		= strcat('Vmatc_TBP',num2str(BandWidthTimeProduct),'TH',num2str(sliceThickness),'Nex',num2str(rawinfo.Nex),'.bin');

        r_filename		= 'r.bin';
        DCW_filename		= 'DCW_MR3.bin';
        Spiral_filename 	= 'SpiralTraj_MR3.bin';
        projID_filename   = 'projID.bin'

        % Compressed dictionary
        %dict = dict.';%Yun Jiang --091321
        temp = single(zeros(size(dict,1)*size(dict,2)*2,1));
        fid = fopen(dict_filename, 'w');
        temp(1:2:end) = real(dict(:));
        temp(2:2:end) = imag(dict(:));
        fwrite(fid, temp, 'float32');
        fclose(fid);
        clear temp;

        fid = fopen(dictnorm_filename,'w');
        temp = single(zeros(size(normeddict,1)*size(normeddict,2)*2,1));
        temp(1:2:end) = real(normeddict(:));
        temp(2:2:end) = imag(normeddict(:));
        fwrite(fid,temp,'float32');
        fclose(fid);clear temp;

        temp = single(zeros(size(Vc,1)*size(Vc,2)*2,1));
        fid = fopen(V_filename, 'w');
        temp(1:2:end) = real(Vc(:));
        temp(2:2:end) = imag(Vc(:));
        fwrite(fid, single(Vc), 'float32');
        fclose(fid);clear temp;

        fid = fopen(r_filename, 'w');
        fwrite(fid,r(:),'int');
        fclose(fid);


        %temp = single(zeros(size(dcf_fid,1)*size(dcf_fid,2)*2,1));
        fid = fopen(DCW_filename,'w');

        %temp(1:2:end) = real(dcf_fid(:));
        %temp(2:2:end) = imag(dcf_fid(:));
        fwrite(fid,dcf_fid(:), 'float32');
        fclose(fid);

        fid = fopen(Spiral_filename,'w');
        SpiralTraj(1,:,:) = kx/max(abs(kx(:)))*0.5;
        SpiralTraj(2,:,:) = ky/max(abs(ky(:)))*0.5;
        fwrite(fid,SpiralTraj(:), 'float32');
        fclose(fid);


        %       fid = fopen(projID_filename,'w');
        % 	  fwrite(fid,projID(1:), 'float32');
        % 	  fclose(fid);

    end

    %dicomwrite(int16(t1map),'T1.dcm'); %JESUS commented
    %dicomwrite(int16(t2map),'T2.dcm'); %JESUS commented
        
    %disp('Size of T1, T2 map:')
    %disp(size(t1map))
    %disp('Size of M0 map:')
    %disp(size(m0map))
    raw = raw_backup;%JESUS
    cd('..');
end %JESUS end
%dicomwrite(reshape(int16(t1map3d),[size(t1map,1) size(t1map,1) 1 size(raw_backup,3)]),'T1.dcm'); % JESUS 
%dicomwrite(reshape(int16(t2map3d),[size(t1map,1) size(t1map,1) 1 size(raw_backup,3)]),'T2.dcm'); % JESUS
%dicomwrite(reshape(int16(m0map3d),[size(t1map,1) size(t1map,1) 1 size(raw_backup,3)]),'T1.dcm'); % JESUS
cd(savedir);
%save("Blurred_Image1.mat","tt_1"); %JESUS
%save("Blurred_Image2.mat","tt_2");
%save("Blurred_Image3.mat","tt_3");
%save("Blurred_Image4.mat","tt_4");


tempfilename = 'PM'; %JESUS modified start
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
tempfilename = strcat(tempfilename,'_Nex',num2str(rawinfo.Nex));
tempfilename = strcat(tempfilename,'.mat');

save(tempfilename,'t1map3d','t2map3d','m0map3d','-v7.3'); %JESUS modified end
disp('Total .mat file saved');
%dicomwrite(Matrix,'checkMatrix.dcm','MultiframeSingleFile',true) JESUS
%dicom save

%%
%% save dictionary to binary
if 0
    %cd('E:\Gadgetron\WIP881v23_DictionaryProstate');
    load('Dict_FISP_TBP8_SliceTH5.mat'); % [3084 48] for 2D prostate MRF
    nset = 3000;
    r0 = 6645;
    ts = 2;

    numpoints = nset*r0;
    interleaved = single(zeros(numpoints*2,1));

    % Alternate real and imaginary data
    interleaved(1:2:end) = real(dict(:));
    interleaved(2:2:end) = imag(dict(:));

    % Write the coordinate data into bin file
    fid1 = fopen('dict.bin','w');
    fwrite(fid1, interleaved, 'float32');
    fclose(fid1);

    %rmat = single(zeros(r0*ts,1));
    rmat = r(:);
    % Write the coordinate data into bin file
    fid1 = fopen('r.bin','w');
    fwrite(fid1, rmat, 'int');
    fclose(fid1);

    dict_permute = dict.'; % N x time
    [U,S,V] = svd(dict_permute,'econ');
    %nk = max(find(diag(S)/S(1)>1e-5));
    nk = 14
    V = V(:,1:nk);
    dictc = (dict_permute*V(:,1:nk)).';

    % compression dictionary
    numpoints = nk*r0;
    interleaved = single(zeros(numpoints*2,1));

    % Alternate real and imaginary data
    interleaved(1:2:end) = real(dictc(:));
    interleaved(2:2:end) = imag(dictc(:));

    % Write the coordinate data into bin file
    fid1 = fopen('dictc.bin','w');
    fwrite(fid1, interleaved, 'float32');
    fclose(fid1);

    dictnorm = sqrt(sum(dictc(:,:).*conj(dictc(:,:))));
    dictcnorm = dictc./repmat(dictnorm,nk,1);
    % save Vmat and dictc
    numpoints = nk*r0;
    interleaved = single(zeros(numpoints*2,1));

    % Alternate real and imaginary data
    interleaved(1:2:end) = real(dictcnorm(:));
    interleaved(2:2:end) = imag(dictcnorm(:));

    % Write the coordinate data into bin file
    fid1 = fopen('dictcnorm.bin','w');
    fwrite(fid1, interleaved, 'float32');
    fclose(fid1);

    numpoints = nset*nk;
    interleaved = single(zeros(numpoints*2,1));

    % Alternate real and imaginary data
    interleaved(1:2:end) = real(V(:));
    interleaved(2:2:end) = imag(V(:));

    % Write the coordinate data into bin file
    fid1 = fopen('Vmatc.bin','w');
    fwrite(fid1, interleaved, 'float32');
    fclose(fid1);
end
Totaltime = toc
