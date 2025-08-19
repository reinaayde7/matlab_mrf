%
% MRF recon Low Rank and Local low Rank giving the .mat files
%
% Based on previous scripts by Rudy Rizzo
% Based on previous scripts by Jesus Fajardo (jesuserf@umich.edu)
% Based on previous scripts by Yun Jiang (yunjiang@med.umich.edu)
clear;clc;close all;
tic
%% setup reconstruction options
% Rudy: we need SVD on otherwise dictionary matching phase is too big
lowrank_SVD             = 1; % 0 :-> 'original' MRF reconstruction
nproj = 48;
sliceToFit = 1;
%% setup path and irt

addpath('./rr_dictionary/DictSimulation_NoSP');
addpath('./dependence');
addpath('./Colormap');
addpath('./Functions_Jesus');
addpath('./OpenSiemensRawData/');

% Path to IRT toolbox
irtPath = 'C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\irt\irt';
addpath(irtPath); run setup;

codePath = 'C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\imagine\'; % ra: I commented momentarily
addpath(genpath(codePath)) % ra: I commented momentarily
%% if import .mat data
path = 'C:\Users\ayde\University of Michigan Dropbox\Reina Ayde\brain_1mmISO\Brain_HV_251217_1mmISO_s1_f300_Nex1000_FAbody15_TE2_8_TR17'
DATA = load(fullfile(path, 'DATA.mat')).DATA;
dict = load(fullfile(path, 'dictCompressed.mat'));
dictCompressed = dict.dictCompressed;
Phi = dict.Phi;
Vc = Phi;
r = dict.r;
coilmap = load(fullfile(path, 'coilmap.mat')).coilmap;
trajectory = load(fullfile(path, 'trajectory.mat'));
k = trajectory.k(:, 1:nproj);
w = trajectory.w;
raw = DATA;
%% Prepare Gridding operator
FOVx = 304;
FOVy = 304;
FOV = 300;
NumFrames = size(Phi, 2);
tempindex =  zeros(1, NumFrames);
NumCoils = size(coilmap, 3)
TotalSlices = 1;
uplimit = size(DATA, 1);
kx =squeeze(real(k));
ky =squeeze(imag(k));
k= double([kx(:) ky(:)]);
mask=true([FOVx FOVy]);
N=size(mask);
nufft_args = {N, [5 5], 2*N, N/2, 'table', 2^12, 'minmax:kb'};
G_gridding = Gmri(k, mask, 'fov', FOV, 'basis', {'dirac'}, 'nufft', nufft_args);
%%
if lowrank_SVD
    rawsvd = zeros(size(raw,1),nproj,size(raw,2),size(dictCompressed,1)); %48 spiral rotations
    disp(size(raw))
    disp('SVD-compressing raw data...')
    f = waitbar(0, 'SVD-compressing raw data...'); tic;
    % idtxtfile = [folderTXT '\FISP_ID.txt'];
    % projID = textread(idtxtfile,'%d\b');
    % idproj = projID(1:size(raw_orig,4))+1;

    projID = repmat(1:nproj, 1, ceil(size(raw,3)/nproj)); %Rudy 2504 - linear spiral sampling independent from txt file
    idproj = projID(1:size(raw,3))';

    %for nn = 1:size(raw,3)
        %waitbar(nn/size(raw,3), f, sprintf('Partition: %d/%d', nn, size(raw,3)));       

    raw_svd_1s = lowrank_2Dksp( squeeze(raw(:,:,:)), idproj, Vc, nproj);
    rawsvd(:,:,:,:) = raw_svd_1s;

    %end
    close(f); %time = toc; % ra I commented the following cause I don't have this function time2clock %SVD_comp_time = time2clock(time);
    raw = rawsvd;
    [SpiralReadout,~, NumCoils,NumFrames] = size(raw);
    disp('Data SVD-compressed!')
    disp(['data size (SVD-compressed): ' num2str(size(raw))])
else
    [SpiralReadout,NumCoils,NumFrames] = size(raw);
end

% % % %% saving for GRIDDING + PF on full Nex to be run on cluster
% % % tempfilename = ['data_GRIDDING_PF_cluster_' date '_' RawData.name(1:13)];
% % % save([RawData.folder tempfilename], 'NumFrames', 'TotalSlices', 'NumCoils', 'log', 'G_gridding', 'kx', 'dcf_fid', 'raw', 'mask', '-v7.3')

%% NUFFT gridding
disp('K-space GRIDDING...')
% Perform NUFFT

G_fid_dcf = Gmri(k, mask, 'fov', N(1), 'basis', {'dirac'}, 'nufft', nufft_args);
dcf_fid = reshape(abs(mri_density_comp_v2(k, 'pipe', 'G', G_fid_dcf.arg.Gnufft)),uplimit,size(kx,2));
dcf_fid = dcf_fid(:, 1:nproj);
adcpad = 20;


img_XYkz = squeeze(single(zeros(FOVx,FOVy,TotalSlices,NumCoils,NumFrames)));
tic
f = waitbar(0, 'K-space GRIDDING...'); tic;
for iFrame= 1:NumFrames %
    waitbar(iFrame/NumFrames, f, sprintf('Time-Point OR SVD-frame: %d/%d', iFrame, NumFrames));
    if lowrank_SVD %idproj becomes the SVD dimension (5) and you use only 1 gridding element instead of all 48
        G = G_gridding;
        dcf = dcf_fid;
        for inc=1:NumCoils

            k_temp = squeeze(raw(:,:,inc,iFrame));  % fullly sampled data
            image_temp = G'*(dcf(:).*k_temp(:)); %dcf(3088x1000) k_temp(3088x6)
            img_XYkz(:,:,inc,iFrame) = embed(image_temp,mask);

        end
    else
        tempindex(:,iFrame) = mod(iFrame - 1,size(kx,2))+1; %variable that matches Nex to spiral arm
        G = G_gridding{tempindex(:,iFrame)};
        dcf = dcf_fid(:,tempindex(:,iFrame));
        for inc=1:NumCoils
            k_temp = squeeze(DATA(:,inc,iFrame));  % single shot
            image_temp = G'*(dcf(:).*k_temp(:));
            img_XYkz(:,:,p,inc,iFrame) = embed(image_temp,mask);
        end
    end
end
ksp_cartesian = transformImageToKspace(transformImageToKspace(img_XYkz,1),2);

close(f); time = toc; %ra I commented the follwing % gridding_time = time2clock(time);
disp('K-space GRIDDED!')
disp(['data size (gridded): ' num2str(size(ksp_cartesian))])
% 60partition, 600 frames (NO SVD) = 48min
%% CSM estimation
% 
image_uncombined = transformKspaceToImage(transformKspaceToImage(ksp_cartesian,2),1);
%image_uncombined = flip(image_uncombined, 3); % This is to match the B0 and B1 ordering - flip in Kz domain

fprintf('3D CSM estimation... \n');
%average through Nex (coils do not depend on the eigenvalue/ex)
image_proxy_coil = sum(image_uncombined,4);

tic
smoothing = 200; %200 by rudy
%chunks = size(ksp_cartesian,3)/6;
csm = estimate_csm_walsh(single(image_proxy_coil), smoothing);
time = toc; %csm_estimate_time = time2clock(time);

% image_proxy_combined=squeeze(sum( conj(csm).*image_proxy_coil, 4 ));
for eig = 1:size(image_uncombined,4)
    image_combined_3D(:, :, eig) = squeeze(sum( conj(csm).*squeeze(image_uncombined(:,:,:,eig)), 3 ));
end

imagine(image_combined_3D) %visual check

% proxy_img = sum(image_combined_3D,3);

%% LLR settings
addpath(genpath('C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\matlab_scripts\CardiacMRF_2D\Reconstruction Code\Gridding'));
addpath(genpath('C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\matlab_scripts\CardiacMRF_2D\Reconstruction Code\Dictionary\LowRankRecon'));
addpath(genpath('C:\Users\ayde\OneDrive - Michigan Medicine\Documents\MATLAB\matlab_scripts\CardiacMRF_2D\Reconstruction Code\Miscellaneous'));
raw_orig_Z = DATA; % originally it's all the DATA.mat of all the slices and fliped in the 3rd axis
readOSFactor=1;
numSpiralArms=nproj;
nr = size(DATA,1);
nex = size(DATA,3);
numCoils = size(DATA,2); 
% i already have w
% w = repmat(dcf_fid, [1,ceil(nex/numSpiralArms)]);
% w = w(:,1:size(raw_orig,4));
wi = w(:, 1:nproj);
Phi = dict.Phi; %Vc;
use_gpu=0;

% I don't have the LowRankRecon
params = setupParameters_LowrankMRF2D( );
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
%% Pattern matching

%fprintf(['Single Slice Fit - Slice n:' num2str(sliceToFit) '\n']);
% -------------------------------------------------------------------------------------------------------------------------------------------------------
fprintf('LR recon\n');
%image_combined = squeeze(image_combined_3D(:,:,sliceToFit,:)); because I
%have only one slice
image_combined = squeeze(image_combined_3D(:,:,:));
[t1map,t2map,b0map,m0map] = patternmatch(image_combined,mask,r,0,dict.dictCompressed(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors

% -------------------------------------------------------------------------------------------------------------------------------------------------------
fprintf('LLR recon \n');
% See these papers for more details...
% Gastao Cruz, MRM 2019. "Sparsity and locally low rank regularization for MR fingerprinting".
% Jesse Hamilton, NMR Biomed 2019. "Simultaneous multislice cardiac magnetic resonance fingerprinting using low rank reconstruction".

%DATA = squeeze(raw_orig_Z(:,:,sliceToFit,:)); % [4004 4 38 600]
DATA = squeeze(raw_orig_Z(:,:,:));
DATA = DATA .* permute(repmat(sqrt(squeeze(w)),[1 1 numCoils]),[1 3 2]);
%coilmap = squeeze(dip.coilmap(:,:,sliceToFit,:));
%coilmap = squeeze(dip.coilmap(:,:,:));
coilmap = coilmap;

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
[t1mapLLR,t2mapLLR,b0mapLLR,m0mapLLR] = patternmatch(images_lowrank,mask,r,0,dict.dictCompressed(1:NumFrames,:),16); %when offres are in blocks 16 helps to avoid RAM outofbound errors

%% single slice plots
addpath('./Colormap');
addpath('./dependence');
load('./Colormap/T2cm.mat');
load('./Colormap/T1cm.mat');

climst1 = [0, 1500]; %3000/1500
climst2 = [0, 300];  %2000/250
s_vert = 40:260;
s_hor = 60:300;
idx1 = s_vert;
idx2 = s_hor;
figure, 
ax1 = subplot(311); imagesc(rot90(cat(1,t1map(idx1,idx2),t1mapLLR(idx1,idx2),t1map(idx1,idx2)-t1mapLLR(idx1,idx2))), climst1), colormap(ax1, T1colormap), colorbar
ax2 = subplot(312); imagesc(rot90(cat(1,t2map(idx1,idx2),t2mapLLR(idx1,idx2),t2map(idx1,idx2)-t2mapLLR(idx1,idx2))), climst2), colormap(ax2, T2colormap), colorbar

m00 = abs(m0map)/max(abs(m0map(:)));
m00LR = abs(m0mapLLR)/max(abs(m0mapLLR(:)));
ax3 = subplot(313); imagesc(rot90(cat(1,m00(idx1,idx2),m00LR(idx1,idx2),m00(idx1,idx2)-m00LR(idx1,idx2))),[0 0.6]), colormap(ax3, gray), colorbar

figure, imagesc(rot90(cat(1,b0map(idx1,idx2),b0mapLLR(idx1,idx2),b0map(idx1,idx2)-b0mapLLR(idx1,idx2))), [-100 100]), colormap(gray), colorbar
