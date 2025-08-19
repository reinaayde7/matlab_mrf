function DeltaB = Load_offresonance(filename)%, maskbase)

% function Offres = OffResMap( filename )
% Script to load Off-Resonance map using
% 2 TE echoes.
% filename is the .dat raw Siemens data
% Maskbase is the through-time MRF average image
% loaded in order of doing a mask.
% This 1.0 version doesn't handle any B0 data shape
% For now just 2D of the Siemens gre_field_mapping sequence
% Jesus Fajardo (jesuserf@med.umich.edu)

set(gca,'DefaultTextFontSize',28);
%clear;clc;close all;
cd('..');
addpath('./dependence');
addpath('./OpenSiemensRawData');
addpath('./Functions_Jesus/');
addpath('./data/');
%addpath('./data/896/');
%cd('..');


%% Open Siemens 2D Raw data
filename = 'data/meas_MID01158_FID75125_field_mapping_FOV400'



raw = mapVBVD(filename);
raw = raw{2}; %not using noise data


raw = squeeze(raw.image( :, :, :, :, 4, :, :, :, :, :, :)); % Select slice in 5th pos
disp(size(raw))
pause(99999999999999)
%disp(size(raw))
data1 = squeeze(raw( :, :, :, 1)); % TE 1 data
data2 = squeeze(raw( :, :, :, 2)); % TE 2 data
disp(size(data1))
data1 = data1(size(data1,1)/4:3*size(data1,1)/4-1,:,:); %crop before zero pad
data2 = data2(size(data2,1)/4:3*size(data2,1)/4-1,:,:);
FOV = 400;
data1pad = zeros(FOV,size(data1,2),FOV);
data2pad = zeros(FOV,size(data2,2),FOV);

disp(size(data1))

for i = 1:size(data1pad,2)
    data1pad(:,i,:) = padarray(squeeze(data1(:,i,:)),[FOV/2-size(data1,1)/2,FOV/2-size(data1,1)/2],0,'both');
    data2pad(:,i,:) = padarray(squeeze(data2(:,i,:)),[FOV/2-size(data2,1)/2,FOV/2-size(data2,1)/2],0,'both');
end

%{
subplot()
imagesc( squeeze(abs(data1pad(:,4,:))));
title('PD Image from MRF')
colorbar
pause(99999999)
%}

%data2 = squeeze(raw( :, :, :, :, 4, :, :, :, :, :, 2)); % TE 2 data
%data1 = squeeze(raw.image( :, :, :, :, :, :, :, 1, :, :, :)); % TE 1 data
%data2 = squeeze(raw.image( :, :, :, :, :, :, :, 2, :, :, :)); % TE 2 data


%Nc =size(data1,2); %Number of coils
[~,Nc,~] = size(data1);

%% Classic FFT
fprintf('FFT...\n')
image1 = zeros(size(data1pad,1),size(data1pad,3),Nc);
image2 = zeros(size(data2pad,1),size(data2pad,3),Nc);

for i = 1:Nc
    image1(:,:,i) = fftshift(ifft2(squeeze(data1pad(:,i,:))));
    image2(:,:,i) = fftshift(ifft2(squeeze(data2pad(:,i,:))));
end

%image1 = fftshift(ifft2(squeeze(data1(:,1,:))));
%disp(size(image1))
%subplot();
%imagesc(abs(squeeze(image2(:,:,1))))


%% Coil Sensitivity Map estimation
fprintf('Step 4: Estimate coil sensitivities and perform the coil combination... \n');

avg_image1 = sum(image1,3);
avg_image2 = sum(image2,3);
csm1 = estimate_csm_walsh(avg_image1);
csm2 = estimate_csm_walsh(avg_image2);

%subplot();
%imagesc(abs(squeeze(csm1(:,:,8))));


%% Perform coil combination
image_combined1 = single(zeros(size(data1pad,1),size(data1pad,3),Nc));
image_combined2 = single(zeros(size(data2pad,1),size(data2pad,3),Nc));

image_combined1 = sum(image1.*conj(csm1),3)./sum(csm1.*conj(csm1),3); 
image_combined2 = sum(image2.*conj(csm2),3)./sum(csm2.*conj(csm2),3); 

image_combined1 = image_combined1(size(image_combined1,1)/4:3*size(image_combined1,1)/4-1,:); % center image and (JESUS added -1)
image_combined2 = image_combined2(size(image_combined2,1)/4:3*size(image_combined2,1)/4-1,:); % remove oversampling

%% Calculate Off-resonance map
Delta_Phi = angle(image_combined2)-angle(image_combined1);


%Delta_Phi = unwrap_phase(Delta_Phi); %unwrap phase
%Delta_Phi = sunwrap(Delta_Phi); %unwrap phase 
Delta_Phi = Unwrap_TIE_DCT_Iter(Delta_Phi); %unwrap phase
%Delta_Phi = qualityGuidedUnwrapping(Delta_Phi); %unwrap phase
Delta_T = (8.61e-3 - 3.85e-3);
wf = 220; %Hz (fat - water difference)
OffRes = Delta_Phi/(2*pi*Delta_T);




%% Load Dicom to compare
[info] = dicominfo('Jesus_research.MR._.8001.1.2022.06.30.11.31.16.152.71364596.dcm');
[B0map] = dicomread('Jesus_research.MR._.8001.1.2022.06.30.11.31.16.152.71364596.dcm');
%B0map = squeeze(B0map(:,:,:,4))%Multislice only!!!!

B0 = rescale(B0map,-1,1)*((2*pi-4096)/2048)*((1/(Delta_T))/2);% (2*p -4096)/2048*range (double(B0map)*(max(B0map(:)) - min(B0map(:)))/4095) *(1/(Delta_T)/2);

%(info.RescaleSlope)%*double(B0map))% +info.RescaleIntercept)/4095*((1/(Delta_T)/2);
I = B0; %dicomread('74689485');
I = double(I);
%I = unwrap_phase(I);
%pause(99999999999)



%% Load through-time MRF average
%  to create the mask
%load("Blurred_Image.mat","image_combined") %in folder 2D_MRF_Reconstruction_main

%C = imquantize(abs(mean(image_combined(:,:,1),3)),500); %mask from MRF through-time average THIS FILE IS IN THE 2D_MRF_ReconstructionMain folder
A = OffRes;
% get some 2D matrix, and plot as surface
%subplot(121), imagesc(C);

%{
disp(size(C))
subplot()
imagesc(squeeze(C(:,:,1)));
title('PD Image from MRF')
colorbar
pause(99999999)
%}

%% Resize Off-resonance map
%# create interpolant
[X,Y] = meshgrid(1:size(A,2), 1:size(A,1));
[XX,YY] = meshgrid(1:size(I,2), 1:size(I,1));
F = scatteredInterpolant(X(:), Y(:), A(:), 'linear');
FF = scatteredInterpolant(XX(:), YY(:), I(:), 'linear');

%# interpolate over a finer grid

[U,V] = meshgrid(linspace(1,size(A,2),FOV), linspace(1,size(A,1),FOV));
[UU,VV] = meshgrid(linspace(1,size(I,2),FOV), linspace(1,size(I,1),FOV));
%subplot(122), imagesc(F(U,V));

InterpolatedField = F(U,V);
LPF = Lowpass(InterpolatedField, 10);
InterpolatedDicom = rot90(FF(UU,VV),3); %This is for the DICOM map to be in the same direction of the acquired map (single slice)
InterpolatedDicom = FF(UU,VV); %This is for the DICOM map to be in the same direction of the acquired map (single slice)
Thres = 70; %Cropping off-res magnitude in +/- Thres
InterpolatedField(InterpolatedField>Thres)=Thres;
InterpolatedField(InterpolatedField<-Thres)=-Thres;
LPF(LPF>Thres)=Thres;
LPF(LPF<-Thres)=-Thres;
InterpolatedDicom(InterpolatedDicom>170)=170;
InterpolatedDicom(InterpolatedDicom<-170)=-170;
%InterpolatedDicom(1:126 , :)= 0;
%InterpolatedDicom(300:end , :)= 0;
%InterpolatedDicom(: , 1:120)= 0;
%InterpolatedDicom(: , 250:end)= 0;

%c = [150 150 300 300];
%r = [150 300 150 300];
%InterpolatedDicom = roipoly(InterpolatedDicom,c,r);

%{
imshowpair(abs(InterpolatedField),abs(squeeze(image_combined(:,:,1))),'diff')


maskidx = C == 2;
disp(size(maskidx));
disp(size(InterpolatedField));
disp(size(InterpolatedField));
masked1 = maskidx .* (squeeze(image_combined(:,:,1)));
masked2 = maskidx .* (squeeze(InterpolatedField(:,:))); %masked fieldmap
masked3 = maskidx .* (squeeze(InterpolatedDicom(:,:))); %masked dicom fieldmap
%}


figure('name','T1 T2 and Proton Density');
subplot(1,2,1)
imagesc(InterpolatedField);
colormap('gray');
title('Calculated map')
colorbar

subplot(1,2,2)
imagesc(LPF);
colormap('gray');
title('Low Pass Filtered Field Map')
colorbar
%{



subplot(1,4,1)
imagesc(masked1);
colorbar
title('Subplot 1')

subplot(1,4,2)
imagesc(masked2);
title('Subplot 2')
colorbar

subplot(1,4,3)
imagesc(maskidx);
title('Subplot 3')
colorbar

subplot(1,4,4)
imagesc(masked3);
title('Subplot 4')
colorbar
%}

%save('off_res.mat','LPF')


end

