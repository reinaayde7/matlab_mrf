% data import and dictionary import:
% run rr_MRF_RUN_v2 up to if lowrank_SVD
dataFolder = 'E:\scanner_data\twix_data\250131\';
save([dataFolder 'demo_SVDcompres_decompress.mat'], 'raw_orig', 'Vc', 'idproj','-v7.3')

%%
dataFolder = 'E:\scanner_data\twix_data\250131\';
load([dataFolder 'demo_SVDcompres_decompress.mat'], 'raw_orig', 'Vc', 'idproj')
% raw_orig: data 
% projID: spiral index
% Vc: compression matrix

size(raw_orig) %raw data size: [Nreadout, Ncoils, Npartitions, Nexcitations]
%%
par = 31; %selected the partition to process
f_o = squeeze(raw_orig(1,1,par,:));
figure, plot(abs(f_o)), title('center kspace'), xlabel('Nex')

%%

% rr = transformImageToKspace(raw_orig,3);
ksp_orig = squeeze(raw_orig(:,:,31,:));

np=48; %number of projections
% ksp_orig = ksp_orig./max(ksp_orig(:));

%%
[ksp_comp, ksvd] = lowrank_2Dksp(ksp_orig,idproj,Vc,np);
[ksp_decompr, ksvd_compr] = lowrank_2Dksp_inv(ksp_comp,Vc,idproj);
[ksp_comp2, ksvd2] = lowrank_2Dksp(ksp_decompr,idproj,Vc,np);

%% plot center of kspace
figure, 
subplot(131), plot(abs(squeeze(ksp_orig(10,1,:)))), title('original')
subplot(132), plot(abs(squeeze(ksp_decompr(10,1,:)))), title('decompressed')
subplot(133), plot(abs(squeeze(ksp_orig(10,1,:))-squeeze(ksp_decompr(10,1,:)) )), title('difference')

%%
figure, 
subplot(131), imagesc(abs(squeeze(ksp_comp(1,1,1,:)))), colorbar
subplot(132), imagesc(abs(squeeze(ksp_comp2(1,1,1,:)))), colorbar
subplot(133), imagesc(abs(squeeze(ksp_comp(1,1,1,:))-squeeze(ksp_comp2(1,1,1,:)) )), colorbar
%%




%%
k =zeros(np,size(ksp_orig,1), size(ksp_orig,2), size(ksp_orig,3));
mask_proj_t = zeros(np,size(ksp_orig,3));
for coil = 1:size(ksp_orig,2) %coils
    for t=1:size(ksp_orig,3) %Nex
        p = idproj(t); % spiral interleaf
        k(p,:,coil,t) = ksp_orig(:,coil,t);
        mask_proj_t(p,t) = 1;
    end
end

%%
size(k) %[proj, Nro, ncoils, 600]
kv = reshape(k, [size(k,1)*size(k,2)*size(k,3), size(k,4)] ); % [everything x 600]
kcv = kv*Vc;% [everything x 7]
kc = reshape(kcv, [size(k,1), size(k,2), size(k,3), size(kcv,2)] );%[proj, Nro, ncoils, 7]

%% decompression test
kcv1 =  reshape(kc, [size(k,1)*size(k,2)*size(k,3), size(kc,4)] );
kv1 = kcv1*Vc'; % [everything x 600]
imtx = Vc*Vc';
k1=reshape(kv1, [size(k,1), size(k,2), size(k,3), size(kv1,2)] );

k1filtered =zeros(np,size(ksp_orig,1), size(ksp_orig,2), size(ksp_orig,3));
for coil = 1:size(ksp_orig,2) %coils
    for t=1:size(ksp_orig,3) %Nex
        p = idproj(t); % spiral interleaf
        k1filtered(p,:,:,t) = k1(p,:,:,t) *mask_proj_t(p,t);
    end
end

%% recompression test
k1fv = reshape(k1filtered, [size(k,1)*size(k,2)*size(k,3), size(k,4)] ); % [everything x 600]
k1fcv = k1fv*Vc;% [everything x 7]
k1fc = reshape(k1fcv, [size(k,1), size(k,2), size(k,3), size(kcv,2)] );%[proj, Nro, ncoils, 7]

%%
figure, 
subplot(121), plot(abs(squeeze(k(1,1,1,:)))), title('original')
subplot(122), plot(abs(squeeze(k1filtered(1,1,1,:)))), title('decompressed')

%%
figure, 
subplot(121), plot(abs(squeeze(kc(1,1,1,:)))), title('original compressed')
subplot(122), plot(abs(squeeze(k1fc(1,1,1,:)))), title('decompressed recompressed')
% subplot(133), plot(abs(squeeze(ksp_orig(10,1,:))-squeeze(ksp_decompr(10,1,:)) )), title('difference')

%% NUFFT gridding
tempindex =  zeros(1, 600);
ksp_kxkyZ = transformImageToKspace(raw_orig,3);
p=31; %slice to process
ksp2D = squeeze(ksp_kxkyZ(:,:,p,:)); %
img_XY = squeeze(single(zeros(log.rawinfo.matrix,log.rawinfo.matrix,NumCoils,NumFrames)));

tic
f = waitbar(0, '2D K-space GRIDDING...'); tic;
for iFrame= 1:NumFrames %
    waitbar(iFrame/NumFrames, f, sprintf('Time-Point OR SVD-frame: %d/%d', iFrame, NumFrames));
    tempindex(:,iFrame) = mod(iFrame - 1,size(kx,2))+1; %variable that matches Nex to spiral arm
    G = G_gridding{tempindex(:,iFrame)};
    dcf = dcf_fid(:,tempindex(:,iFrame));
    for inc=1:NumCoils       
        k_temp = squeeze(ksp2D(:,inc,iFrame));  % single shot
        image_temp = G'*(dcf(:).*k_temp(:));
        img_XY(:,:,inc,iFrame) = embed(image_temp,mask);        
    end
end
close(f); time = toc; gridding_time = time2clock(time);
disp('2D K-space GRIDDED!')
disp(['data size (gridded): ' num2str(size(img_XY))])
% 60partition, 600 frames (NO SVD) = 48min

%%

for i=1:304
   img_XY_comp(i,:,:,:,:) = lowrank_2Dksp(squeeze(img_XY(i,:,:,:)),idproj,Vc,np);
end

%%
size(img_XY_comp)

%%
imagine(sum(img_XY_comp,3))
