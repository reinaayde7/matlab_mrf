
%%
tempfilename = 'data_GRIDDING_PF_cluster_07-Feb-2025_meas_MID00029';
RawData.folder = 'E:\scanner_data\twix_data\250131\';
load([RawData.folder tempfilename]);

disp('K-space GRIDDING - noLOWRANK...')
% Perform NUFFT

tempindex =  zeros(1, NumFrames);
img_XYkz = squeeze(single(zeros(log.rawinfo.matrix,log.rawinfo.matrix,TotalSlices,NumCoils,NumFrames)));
tic
f = waitbar(0, 'K-space GRIDDING...'); tic;
for iFrame= 1:NumFrames %
    waitbar(iFrame/NumFrames, f, sprintf('Time-Point OR SVD-frame: %d/%d', iFrame, NumFrames));
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
ksp_cartesian = transformImageToKspace(transformImageToKspace(img_XYkz,1),2);
close(f); time = toc; gridding_time = time2clock(time);
clearvars img_XYkz


%% TODO -> set it up on the cluster

disp('Partial Fourier reconstruction...')
    % ksp_cartesian_full = ksp_cartesian;
    f = waitbar(0, 'run POCS...'); tic;
    tic
    for svd_i = 1:size(raw,5)
        waitbar(svd_i/size(raw,5), f, sprintf('Time-Point OR SVD-frame: %d/%d', svd_i, size(raw,5)));
        ksp = permute(squeeze(ksp_cartesian(:,:,:,:,svd_i)), [4,1,2,3]);
        iter = 20;
        watchProgress = 0;
        [~, kspFull] = pocs( ksp, iter, watchProgress );
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
