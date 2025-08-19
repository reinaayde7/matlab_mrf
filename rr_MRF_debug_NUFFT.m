
size(raw)
% nRO, nCoil, nPar, Nex

% rr = raw(1+adcpad:end-adcpad,:,:,:);
% rr = rr(1:uplimit,:,:,:);
rr = raw;
size(rr)

% DATA_FT = ifft1c(rr, 3); % FFT along partitions (kz --> z) no density compensation since fully sampled
DATA_FT = rr;
[nr,nc,nz,Nex] = size(rr);
np = 48;
% Nex=1000;
idproj = repmat([1:48],1,100);
%%
parfor z = 1:nz % loop over all partitions
    ktot = zeros(nr,np,nc); % sum the k-space data over time
    dcf = zeros(1,np);               % account for uneven sampling of spiral projections

    for t = 1:Nex  % loop over measurement frames
        p = idproj(t); % spiral arm index
        ktot(:,p,:) = squeeze(ktot(:,p,:)) + DATA_FT(:,:,z,t); % add to running sum of k-space
        dcf(p) = dcf(p) + 1;
    end

    dcf(dcf==0) = 1;
    ktot = ktot./repmat(dcf,[nr 1 nc]); % account for uneven sampling of 48 spiral arms
    ktot = ktot.*repmat(sqrt(densityComp/max(densityComp(:))),[1 1 nc]); % density compensation function
    
    FT = NUFFT(kxall/N(1)+1i*kyall/N(1),densityComp/max(densityComp(:)),[0 0],[N(1)*readoutOS N(1)*readoutOS]);
    image_coil(:,:,z,:) = FT'*ktot;
end

%%
size(image_coil)
temp = makesos(image_coil,4);
figure, sliceViewer(temp)

%%
figure, 
subplot(131),imagesc(squeeze(temp(:,:,4))), colorbar
subplot(132),imagesc(squeeze(temp(:,:,14))), colorbar
subplot(133),imagesc(squeeze(temp(:,:,24))), colorbar