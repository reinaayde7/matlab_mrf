% Rudy v1 240830
% optimized over a 2D recon from the 3D acquistiion
% requires the first part of the recon from rr_MRF_itSENSE until data are
% prepped.
addpath(genpath('E:\POSTDOC_UoM\02_Project_PE_lung_MRI_055T\Jesse-RealTimeCine_DeepImagePrior-main-d2310\Reconstruction_Code\3D\MoCom_recon\code'))
addpath(genpath('E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\PCAtensordenoiser'))

% parameters
% ITSENSE
itSense.maxit = 4; 
itSense.limit = 1E-8;

% ADMM (ITSENSE + Tensor MP-PCA denoiser)
itSenseTMPPCA.ADMM_maxit = 15; % number of iterations in outer loop of ADMM recon
itSenseTMPPCA.ADMM_CG_maxit = 2; % number of CG iterations (don't need many since it has warm start)
itSenseTMPPCA.ADMM_CG_limit = 1E-20; % can also stop iterations by residual, but using a fixed number is easier
%parameter looping
itSenseTMPPCA.ADMM_CG_lambda = [5e-2]; %regularizer weight
itSenseTMPPCA.TMPPCAsiz = 5; %kernel size for patch denoising



%%
% slice = 16;
% raw = transformImageToKspace(raw,3); %transform back to kspace
raw_kspace = fftshift(ifft(ifftshift(raw),[],3));
data = permute(raw_kspace,[1,2,3,5,4]); % [readout, spirals, slices, EIG, coils]
dataCS = data/max(abs(data(:)));
% dataCS2D = squeeze(dataCS(:,:,slice,:,:));
csmSENSE = permute(repmat(csm,[1 1 1 1 size(data,4)]), [1,2,3,5,4]);
% csmSENSE2D = squeeze(csmSENSE(:,:,slice,:,:));
siz = size(csmSENSE); siz = siz(1:4); %siz(4) = size(data,5);
Ksiz = size(dataCS);
% fully sampled data
masks = ones([Ksiz(2) Ksiz(3)]);
ksp_dcfXY = dcf_fid;
spirals = kx+1i*ky;
ksp_loc = spirals;
use_gpu=1;
SENSE_op = itSENSE_SoS_GPU_MRF(masks,csmSENSE,ksp_loc,ksp_dcfXY,siz,Ksiz,use_gpu);

%% test 
tt = SENSE_op'*(dataCS.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4) Ksiz(5)])));
%%
% itSENSE recon
tic
[itSENSE_recon,residuals_sense,sense_its] = Conjugate_Gradient((dataCS.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4) Ksiz(5)]))),...
    SENSE_op,itSense.maxit,itSense.limit);
time2clock(toc);
out.itSense.recon_its = sense_its;
out.itSense.maxit = itSense.maxit;
out.itSense.limit = itSense.limit;
out.itSense.res = residuals_sense;
out.mrfStandard_recon = squeeze(image_combined);

%%
% ADMM (ITSENSE + Tensor MP-PCA denoiser) recon 
tic
dataADMM = dataCS;
itSENSE_recon0 = sense_its(:,:,:,:,2); %pick a resonable starting point (no to much noise enhancement + ok undersmp artifact)
% itSENSE_recon0 = itSENSE_recon0/ max(abs(itSENSE_recon0(:)));

for i1 = 1:size(itSenseTMPPCA.ADMM_CG_lambda,2)
    disp(['i1: ' num2str(i1)])
    admmCGl = itSenseTMPPCA.ADMM_CG_lambda(i1);

    
    [ADMM_recon,tMPPCA_its,noise_its,ADMM_its,res_it] = ADMM_itSENSE_TMPPCA_MRF_solver_3D((dataADMM.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4) Ksiz(5)]))),...
                                            SENSE_op,itSenseTMPPCA.ADMM_maxit,admmCGl,itSenseTMPPCA.ADMM_CG_maxit,itSenseTMPPCA.ADMM_CG_limit,itSENSE_recon0,itSenseTMPPCA.TMPPCAsiz);


    out.admm.log.admmCGl{i1} = admmCGl;
    out.admm.recon_its{i1} = ADMM_its;
    out.admm.log.res{i1} = res_it;
    out.admm.log.noise{i1} = noise_its;
end
time2clock(toc);


%% straight denoising
denois_straight = denoise_recursive_tensor(tt,[5 5]);
out.denois_straight = denois_straight;

%% savings
saveName = 'FISP_volunteer_FreeMax_3D_CS_TMPPCAdenoising_v2.mat';
save([saveFolder saveName], 'out','-v7.3')

%%
size(out.admm.recon_its{1})
%%
imagine(squeeze(out.admm.log.noise{1}(:,:,16,:,:)))