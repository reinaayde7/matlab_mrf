% % % %%
% % % raw_orig = transformImageToKspace(raw_orig,3); %Rudy: FFT on kz dimension (from kspace to image space)
% % % raw_orig = flip(raw_orig, 3);
% % % size(raw_orig)
% % % dip.kspace_data = squeeze(raw_orig(:,:,11,:));
% % % 
% % % %%
% % % 
% % % % JESUS start
% % % % for loop iterating over N slices JESUS
% % % % tt_1 = zeros(size(raw,3),size(dict,1)); %each of these variables is to save a through-time measurement in certain pixel JESUS
% % % % tt_2 = zeros(size(raw,3),size(dict,1));
% % % % tt_3 = zeros(size(raw,3),size(dict,1));
% % % % tt_4 = zeros(size(raw,3),size(dict,1));
% % % 
% % % raw = flip(raw, 3); % This is to match the B0 and B1 ordering - flip in Kz domain


destFolder = 'E:\POSTDOC_UoM\08_Project_MRF\DIP\';
Folder = 'TRUEFISP\';
what = 'NIST_T2L';
extraN = '250109_f300_Nex400_TE1_9_TR13_FAbody_shifted';

%%
%raw = circshift(raw,16,3);%This is just for some patients
raw_backup = raw;% JESUS
% b1index = 0;
% for nn = 6:13%1:size(raw,3),
    nn = 14;
    disp('Slice number:')
    nn
    
    
    dip.kspace_data = squeeze(raw_orig(:,:,nn,:));
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
    dip.coilmap = csm;
    
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

% end
    
% % image_combined = real(image_combined) - j*imag(image_combined);
   %
    % if isTRUEFISP
    %     image_combined = image_combined*exp(j*pi); %to compensate for error in simulation wehre dictionary is conj
    % end

    % save for DIP
    % subFolderName =['NIST_T2L_s' num2str(nn) '_241218_f300_ISO1mm_Nex1000_TE2_8_TR14_FAbody15'];
    subFolderName =[what '_s' num2str(nn) '_' extraN];
    % subFolderName =['carotid_volume2_s' num2str(nn) '_241218_f300_Nex1000_FAbody15'];
    mkdir(fullfile(destFolder,Folder,subFolderName))
    
    coilmap = dip.coilmap;
    save([destFolder Folder subFolderName '/coilmap.mat'], 'coilmap', '-v7.3'); %[Nx, Ny, coils]
    
    DATA = dip.kspace_data;  
    save([destFolder Folder subFolderName '/DATA.mat'],'DATA','-v7.3'); %[Nro, coils, Nex];
    
    dictCompressed = dict;
    Phi = Vc;
    save([destFolder Folder subFolderName '/dictCompressed.mat'], 'Phi', 'dictCompressed', 'r', '-v7.3'); %Phi: [Nex,EIG], dictCompressed: [EIG,DicSIZE], r: [DicSIZE, DicVAR]
    
    clearvars k w
    k(:,1,:) = dip.k;
    w(:,1,:) = dip.w;
    save([destFolder Folder subFolderName '/trajectory.mat'], 'k', 'w', '-v7.3'); %k: [Nro,1,Nex], w:[Nro,1,Nex]
    
    disp('DIP savings DONE!')
% end