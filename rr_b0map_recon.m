clear all, close all, clc
dirName = 'E:\scanner_data\twix_data\241018\';
fileName = 'meas_MID00343_FID113798_b0_mapping_rudy_200x200_st3_5';

dirName = 'E:\scanner_data\twix_data\240731\';
fileName = 'meas_MID00073_FID100464_b0_mapping_rudy_0_55T';
toImport = [dirName fileName '.dat'];

% recon based on 

% function Offres = OffResMap( filename )x  
% Script to load Off-Resonance map using
% 2 TE echoes.
% filename is the .dat raw Siemens data
% Maskbase is the through-time MRF average image
% loaded in order of doing a mask.
% This 1.0 version doesn't handle any B0 data shape
% for multislice B0 maps
% Jesus Fajardo (jesuserf@med.umich.edu)

% set(gca,'DefaultTextFontSize',28);
%clear;clc;close all;
%cd('./data/');
addpath('./dependence');
addpath('./OpenSiemensRawData');
codePath = 'E:\POSTDOC_UoM\05_MATLAB\Imagine_old\';
addpath(genpath(codePath))
addpath('./Functions_Jesus/');
%addpath('./data/');
%addpath('./data/896/');
%cd('..');


%% Open Siemens 2D Raw data
raw = mapVBVD(toImport);
raw = raw{2}; %not using noise data

%%
raw = squeeze(raw.image( :, :, :, :, :, :, :, :, :, :, :)); % Select slice in 5th pos
disp(size(raw))
%rudy, with one slice: Kx, Ncoils, Ky, NAve, NTE

FOV = 400; %Field of view
LPF = zeros(FOV,FOV);

if numel(size(raw)) < 5 %handling acquisition of 1slice only
    rawn(:,:,:,1,:) = raw;
    clearvars raw
    raw = rawn;
end
disp(size(raw))

%%
for nn = 1:size(raw,4), %num of slices
% for nn = 20:30
    disp(['Slice: ' num2str(nn)])
    % data1 = squeeze(raw( :, :, :, 1)); % TE 1 data
    % data2 = squeeze(raw( :, :, :, 2)); % TE 2 data
    % 
    %when you acquire multiple slices
    data1 = squeeze(raw( :, :, :, nn, 1)); % TE 1 data
    data2 = squeeze(raw( :, :, :, nn, 2)); % TE 2 data

    %averaging
    % data1 = squeeze(sum(raw(:,:,:,:,1),4)/size(raw,4));
    % data2 = squeeze(sum(raw(:,:,:,:,2),4)/size(raw,4));

    %data2 = squeeze(raw( :, :, :, :, 4, :, :, :, :, :, 2)); % TE 2 data
    %data1 = squeeze(raw.image( :, :, :, :, :, :, :, 1, :, :, :)); % TE 1 data
    %data2 = squeeze(raw.image( :, :, :, :, :, :, :, 2, :, :, :)); % TE 2 data


    %Nc =size(data1,2); %Number of coils
    [~,Nc,~] = size(data1);

    % Classic FFT
    fprintf('FFT...\n')
    image_coil_te1 = zeros(FOV*2,FOV,Nc);
    image_coil_te2 = zeros(FOV*2,FOV,Nc);

    for i = 1:Nc
        image_coil_te1(:,:,i) = fftshift(ifft2(ChangeArraySize(squeeze(data1(:,i,:)),[FOV*2,FOV])));
        image_coil_te2(:,:,i) = fftshift(ifft2(ChangeArraySize(squeeze(data2(:,i,:)),[FOV*2,FOV])));
    end

    % Coil Sensitivity Map estimation
    fprintf('Step 4: Estimate coil sensitivities and perform the coil combination... \n');

    avg_image1 = sum(image_coil_te1,3); %sum over slices I assume.
    avg_image2 = sum(image_coil_te2,3);
    % csm1 = estimate_csm_walsh(image_coil_te1);
    % csm2 = estimate_csm_walsh(image_coil_te2);
    
    csm1 = estimate_csm_walsh(avg_image1);
    csm2 = estimate_csm_walsh(avg_image2);

%
    % Perform coil combination
    image1 = single(zeros(size(data1,1),size(data1,3),Nc));
    image2 = single(zeros(size(data2,1),size(data2,3),Nc));

    image1 = sum(image_coil_te1.*conj(csm1),3)./sum(csm1.*conj(csm1),3);
    image2 = sum(image_coil_te2.*conj(csm2),3)./sum(csm2.*conj(csm2),3);

    image1 = image1((2*FOV)/4+1:3*(2*FOV)/4,:); %
    image2 = image2((2*FOV)/4+1:3*(2*FOV)/4,:);
    %image_combined1 = image_combined1(size(image_combined1,1)/4:3*size(image_combined1,1)/4-1,:); % center image and (JESUS added -1)
    %image_combined2 = image_combined2(size(image_combined2,1)/4:3*size(image_combined2,1)/4-1,:); % remove oversampling

    % Calculate Off-resonance map
    Delta_Phi = angle(image2)-angle(image1);
    % figure, imagesc(Delta_Phi)

    %Delta_Phi = unwrap_phase(Delta_Phi); %unwrap phase
    %Delta_Phi = sunwrap(Delta_Phi); %unwrap phase
    Delta_Phi = Unwrap_TIE_DCT_Iter(Delta_Phi); %unwrap phase
    %Delta_Phi = qualityGuidedUnwrapping(Delta_Phi); %unwrap phase
    %Delta_T = (22.94e-3 - 10.00e-3); % 0.55 T
    %Delta_T = (8.61e-3 - 3.85e-3); % 1.5 T
    
    %wf = 220; %Hz (fat - water difference) 1.5T

    % 3T
    % Delta_T = (11 - 6.00)*1e-3;
    % Delta_T = (8.68 - 6.00)*1e-3;
    % wf = 430; %Hz 

    % 1.5T
    % Delta_T = (11 - 6.00)*1e-3;
    % wf = 220; %Hz 

    % 0.55T
    % Delta_T = (18.24 - 6.00)*1e-3;
    Delta_T = (12.94)*1e-3;
    wf = 80; %Hz 

    
    OffRes = Delta_Phi/(2*pi*Delta_T);

    % Flip the map to match the reconstruction orientation
    OffRes = rot90(OffRes,1);

    LPF_2D = Lowpass(OffRes, 15); %smoothing
    
    %Cropping off-res magnitude in +/- Thres (aliasing!)
    Thres = wf/2; 
    OffRes(OffRes>Thres)=Thres;
    OffRes(OffRes<-Thres)=-Thres;
    LPF_2D(LPF_2D>Thres)=Thres;
    LPF_2D(LPF_2D<-Thres)=-Thres;

    LPF_2D = rot90(LPF_2D,3);
    % figure, imagesc(OffRes)
    % figure, imagesc(LPF_2D), colorbar
    %hereafter to consider only when acquiring multiple slices.
    LPF(:,:,nn) = LPF_2D;
    % LPF(:,:,1) = LPF_2D;
end

%% to run this outside the slice for loop
A = LPF(:,:,1:int16(size(LPF,3)/2));
B = LPF(:,:,int16(size(LPF,3)/2)+1:size(LPF,3));

disp('sizes')
disp(size(A))
disp(size(B))
disp(size(LPF))

idxeven = 0;
idxodd = 0;
for i=1:size(LPF,3),
    if bitget(i,1) %true if isODD
        if idxeven < size(B,3)
        LPF(:,:,i) = A(:,:,idxeven+1);
        idxeven = idxeven+1
        end
    else
        if idxodd < size(A,3)
        LPF(:,:,i) = B(:,:,idxodd+1);
        idxodd = idxodd+1
        end
    end
end
% save('off_res.mat','LPF')

%%
imagine(LPF)
%%  
saveName = [dirName fileName '_OFFRES.mat'];
save(saveName,'LPF')



