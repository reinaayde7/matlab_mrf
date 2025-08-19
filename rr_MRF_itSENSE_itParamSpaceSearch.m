% Rudy v1 240830
% optimized over a 2D recon from the 3D acquistiion
% requires the first part of the recon from rr_MRF_itSENSE until data are
% prepped.
addpath(genpath('E:\POSTDOC_UoM\02_Project_PE_lung_MRI_055T\Jesse-RealTimeCine_DeepImagePrior-main-d2310\Reconstruction_Code\3D\MoCom_recon\code'))
addpath(genpath('E:\POSTDOC_UoM\08_Project_MRF\3D_MRF_FISP_Prostate-main_d240531\PCAtensordenoiser'))

% parameters
% ITSENSE
itSense.maxit = 10; 
itSense.limit = 1E-8;

% ADMM (ITSENSE + BM3D)
% itSenseBM3D.ADMM_maxit = 10; % number of iterations in outer loop of ADMM recon
% itSenseBM3D.ADMM_CG_maxit = 2; % number of CG iterations (don't need many since it has warm start)
% itSenseBM3D.ADMM_CG_limit = 1E-20; % can also stop iterations by residual, but using a fixed number is easier
% %parameter looping
% itSenseBM3D.ADMM_CG_lambda = [1e-3 1e-2 1e-1];
% itSenseBM3D.bm4d_lambda = [0.1 0.5 1 2 5]; % 5e-3 1e-2 2e-2]; %lambda denoiser regularization (the higher the smoother)

% ADMM (ITSENSE + Tensor MP-PCA denoiser)
itSenseTMPPCA.ADMM_maxit = 10; % number of iterations in outer loop of ADMM recon
itSenseTMPPCA.ADMM_CG_maxit = 2; % number of CG iterations (don't need many since it has warm start)
itSenseTMPPCA.ADMM_CG_limit = 1E-20; % can also stop iterations by residual, but using a fixed number is easier
%parameter looping
itSenseTMPPCA.ADMM_CG_lambda = [1e-3 5e-2 1e-2 5e-2 1e-1]; %regularizer weight
itSenseTMPPCA.TMPPCAsiz = 5; %kernel size for patch denoising

slice = 16;
data = permute(raw_backup,[1,2,3,5,4]); % [readout, spirals, slices, EIG, coils]
dataCS = data/max(abs(data(:)));
dataCS2D = squeeze(dataCS(:,:,slice,:,:));
csmSENSE = permute(repmat(csm,[1 1 1 1 size(data,4)]), [1,2,3,5,4]);
csmSENSE2D = squeeze(csmSENSE(:,:,slice,:,:));
siz = size(csmSENSE2D); siz = siz(1:3); %siz(4) = size(data,5);
Ksiz = size(dataCS2D);
% fully sampled data
masks = ones([Ksiz(2) Ksiz(3)]);
ksp_dcfXY = dcf_fid;
spirals = kx+1i*ky;
ksp_loc = spirals;
use_gpu=1;
SENSE_op2D = itSENSE_SoS_GPU_MRF_2D(masks,csmSENSE2D,ksp_loc,ksp_dcfXY,siz,Ksiz,use_gpu);

%%
% itSENSE recon
tic
[itSENSE_recon,residuals_sense,sense_its] = Conjugate_Gradient((dataCS2D.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4)]))),...
    SENSE_op2D,itSense.maxit,itSense.limit);
time2clock(toc);
out.itSense.recon_its = sense_its;
out.itSense.maxit = itSense.maxit;
out.itSense.limit = itSense.limit;
out.itSense.res = residuals_sense;
out.mrfStandard_recon = squeeze(image_combined(:,:,slice,:));

%%
% ADMM (ITSENSE + BM3D) recon 
tic
dataADMM = dataCS2D/max(abs(dataCS2D(:)));
itSENSE_recon0 = sense_its(:,:,:,2); %pick a resonable starting point (no to much noise enhancement + ok undersmp artifact)
itSENSE_recon = itSENSE_recon/ max(abs(itSENSE_recon(:)));

for i1 = 1:size(itSenseBM3D.ADMM_CG_lambda,2)
    for i2 = 1:size(itSenseBM3D.bm4d_lambda,2)
        disp(['i1: ' num2str(i1) ' - i2: ' num2str(i2)])
        admmCGl = itSenseBM3D.ADMM_CG_lambda(i1);
        bm4dl = itSenseBM3D.bm4d_lambda(i2);

        
        [ADMM_recon,bm3d_its,noise_its,ADMM_its,res_it] = ADMM_itSENSE_BM4D_MRF_solver_2D((dataADMM.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4)]))),...
                                                SENSE_op2D,itSenseBM3D.ADMM_maxit,admmCGl,itSenseBM3D.ADMM_CG_maxit,itSenseBM3D.ADMM_CG_limit,itSENSE_recon0,bm4dl);


        out.admm.log.admmCGl{i1,i2} = admmCGl;
        out.admm.log.bm4dl{i1,i2} = bm4dl;
        out.admm.recon_its{i1,i2} = ADMM_its;
        out.admm.log.res{i1,i2} = res_it;
    end
end
time2clock(toc);

%%
% ADMM (ITSENSE + Tensor MP-PCA denoiser) recon 
tic
dataADMM = dataCS2D/max(abs(dataCS2D(:)));
itSENSE_recon0 = sense_its(:,:,:,2); %pick a resonable starting point (no to much noise enhancement + ok undersmp artifact)
itSENSE_recon = itSENSE_recon/ max(abs(itSENSE_recon(:)));

for i1 = 1:size(itSenseTMPPCA.ADMM_CG_lambda,2)
    disp(['i1: ' num2str(i1)])
    admmCGl = itSenseTMPPCA.ADMM_CG_lambda(i1);
    
    [ADMM_recon,tMPPCA_its,noise_its,ADMM_its,res_it] = ADMM_itSENSE_TMPPCA_MRF_solver_2D((dataADMM.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4)]))),...
                                            SENSE_op2D,itSenseTMPPCA.ADMM_maxit,admmCGl,itSenseTMPPCA.ADMM_CG_maxit,itSenseTMPPCA.ADMM_CG_limit,itSENSE_recon0,itSenseTMPPCA.TMPPCAsiz);


    out.admm.log.admmCGl{i1} = admmCGl;
    out.admm.recon_its{i1} = ADMM_its;
    out.admm.log.res{i1} = res_it;
end
time2clock(toc);



%% savings
saveName = 'FISP_patient_Recon_v2_slice19.mat';
save([saveFolder saveName], 'out')
%%
imagine(out.admm.recon_its{1,4})



%%
imagine(permute(sense_its,[2,1,4,3]))




%%
% plain BM3D denoising of direct recon
dir_recon = SENSE_op2D'*(dataCS2D.*sqrt(repmat(ksp_dcfXY,[1 1 Ksiz(3) Ksiz(4)])));

%%
bm3d_power = [0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5];
for dp = 1:length(bm3d_power) %denoising power
    slc = dir_recon; %[Nx, Ny, EIG]
    norm = max(abs(slc(:)))/100; %Rudy: crucial scaling for BM3D
    slc = slc / norm;
     for e = 1:size(dir_recon,3)
        dd(:,:,e) = (BM3D(real(squeeze(slc(:,:,e))), bm3d_power(dp)) + 1i*BM3D(imag(squeeze(slc(:,:,e))), bm3d_power(dp)))*norm;
     end
     dir_recon_denoised(:,:,:,dp) = dd;
end

%%
ddd = zeros([400,400,6,8]);
ddd(:,:,:,1) = dir_recon;
ddd(:,:,:,2:end) = dir_recon_denoised;
imagine(permute(ddd,[1,2,4,3]))


%%
out.recon = ddd;
out.denoising_power = [0,bm3d_power];
%%
saveName = 'FISP_healthyVolunteer_direct_BM3D_denoiser_slice16_v3.mat';
save([saveFolder saveName], 'out')