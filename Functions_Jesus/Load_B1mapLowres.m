function DeltaB1 = Load_B1map(filename)%, maskbase)

% Script to load B1 map using TFL
% filename is the .dat raw Siemens data
% Maskbase is the through-time MRF average image
% loaded in order of doing a mask.
% This 1.0 version doesn't handle any B1 data shape
% For now just 2D of the Siemens tfl_b1map sequence
% Jesus Fajardo (jesuserf@med.umich.edu)

set(gca,'DefaultTextFontSize',28);
%clear;clc;close all;
cd('..');
addpath('./dependence');
addpath('./OpenSiemensRawData');
addpath('./Functions_Jesus/');
addpath('./data/');
%addpath('./data/896/');Sli
%cd('..');


%% Open Siemens 2D Raw data
filename = 'data/meas_MID05161_FID52101_tfl_b1map_FOV400'
    

raw = mapVBVD(filename);
raw = raw{2}; %not using noise data
raw.image
raw = squeeze(raw.image( :, :, :, :,:, :, :, :, :, :, :)); % Select slice in 4th pos
%raw = circshift(raw,int16(size(raw,4)/2),4);

FOV = 400; %Field of view
LPF = zeros(FOV,FOV,size(raw,4));

for nn = 1:size(raw,4),
%raw = squeeze(raw( :, :, :, nn,:, :)); % When is not a single slice

disp(size(raw))
FOV = 64; %Field of view
dataFA1 = squeeze(raw(:,:,:,:,1));
dataFA2 = squeeze(raw(:,:,:,:,2));
dataFA3 = squeeze(raw(:,:,:,:,3));
data1 = squeeze(dataFA2( :, :, :, 1)); % SS Pre data
data2 = squeeze(dataFA3( :, :, :, 2)); % PD data CHEQUEAR ESTO Jesus
%Nc =size(data1,2); %Number of coils


[~,Nc,~] = size(data1);
disp(size(data1))
%data1 = data1(1+size(data1,1)/4:3*size(data1,1)/4-1,:,:);
%data2 = data2(1+size(data2,1)/4:3*size(data2,1)/4-1,:,:);
disp(size(data1))


%% Classic FFT
fprintf('FFT...\n')
image1 = zeros(FOV*2,FOV,Nc);
image2 = zeros(FOV*2,FOV,Nc);

for i = 1:Nc
    image1(:,:,i) = fftshift(ifft2(ChangeArraySize(squeeze(data1(:,i,:)),[FOV*2,FOV])));
    image2(:,:,i) = fftshift(ifft2(ChangeArraySize(squeeze(data2(:,i,:)),[FOV*2,FOV])));
end


%% Coil Sensitivity Map estimation
fprintf('Step 4: Estimate coil sensitivities and perform the coil combination... \n');

avg_image1 = sum(image1,3);
avg_image2 = sum(image2,3);
csm1 = estimate_csm_walsh(avg_image1);
csm2 = estimate_csm_walsh(avg_image2);

%% Perform coil combination
image_combined1 = single(zeros(size(data1,1),size(data1,3),Nc));
image_combined1 = sum(image1.*conj(csm1),3)./sum(csm1.*conj(csm1),3); 
%image_combined1 = image_combined1(18:81,:); % center image THESE DIMENSIONS ARE FOR 100x64 IMGS
image_combined2 = single(zeros(size(data2,1),size(data2,3),Nc));
image_combined2 = sum(image2.*conj(csm1),3)./sum(csm2.*conj(csm2),3); 
image_combined1 = image_combined1((2*FOV)/4+1:3*(2*FOV)/4,:); %
image_combined2 = image_combined2((2*FOV)/4+1:3*(2*FOV)/4,:);

%image_combined1 = image_combined1(size(image_combined1,1)/4:3*size(image_combined1,1)/4,1:size(image_combined1,2)); % center image and (JESUS added -1)
%image_combined2 = image_combined2(size(image_combined2,1)/4:3*size(image_combined2,1)/4,1:size(image_combined2,2)); % remove oversampling
%image_combined2 = image_combined2(18:81,:); % center image THESE DIMENSIONS ARE FOR 100x64 IMGS

%% Calculate kmap (Correction coefficient map)
ratio = abs(image_combined2)./abs(image_combined1);
Thres = 1.; %Cropping off-res magnitude in +/- Thres
ratio(ratio>Thres)=Thres;
ratio(ratio<-Thres)=-Thres;
%this threshold above is to prevent the cosine to be complex
FA=80.; % CHECK FLIP-ANGLE!
b1map = acosd(ratio)/FA;

%% Unwrap
%b1map = Unwrap_TIE_DCT_Iter(b1map); %unwrap phase
%% Flip the map to match the reconstruction orientation
%b1map = rot90(b1map,1);
%{
%# interpolate over a finer grid
A = b1map;
[X,Y] = meshgrid(1:size(A,2), 1:size(A,1));
%[XX,YY] = meshgrid(1:size(I,2), 1:size(I,1));
F = scatteredInterpolant(X(:), Y(:), A(:), 'linear');
FOV = 400;
[U,V] = meshgrid(linspace(1,size(A,2),FOV), linspace(1,size(A,1),FOV));
%[UU,VV] = meshgrid(linspace(1,size(I,2),FOV), linspace(1,size(I,1),FOV));
%subplot(122), imagesc(F(U,V));
%}
LPF_B1 = b1map;%F(U,V);

%{
%% Crop handrawn field
quant_map = InterpolatedField;
size_image_x = size(quant_map,1);
size_image_y = size(quant_map,2);
figure, imagesc(1:size_image_x, 1:size_image_y, quant_map);
h = drawfreehand;
ROI = createMask(h);
MaskedImage = ROI.*quant_map;
data = nonzeros(MaskedImage);
%histogram(data,15)
%imshow(MaskedImage,[]);
InterpolatedField = MaskedImage;


LPF_B1_2D = Lowpass(InterpolatedField, 15);
ThresHigh = 1.5; %Cropping off-res magnitude in +/- Thres
ThresLow = 0.;
InterpolatedField(InterpolatedField>ThresHigh)=ThresHigh;
InterpolatedField(InterpolatedField<ThresLow)=ThresLow;
LPF_B1_2D(LPF_B1_2D>ThresHigh)=ThresHigh;
LPF_B1_2D(LPF_B1_2D<ThresLow)=ThresLow;
cmlist = [ThresLow, ThresHigh];
%}
%{
figure('name','T1 T2 and Proton Density');
subplot(1,2,1)
imagesc(rot90(((InterpolatedField(:,:))),1),cmlist);
%colormap('gray');
title('Calculated map')
colorbar
daspect([1 1 1]) 
%}
%LPF_B1_2D = rot90(LPF_B1_2D,3); %Rotate to match the coefficients images orientation
%mr_imshow(LPF_B1_2D(:,:));

subplot()


subplot()
imagesc((((LPF_B1(:,:)))),cmlist);
%colormap('gray');
title('Low Pass Filtered Field Map')
colorbar
daspect([1 1 1]) 
end

save('B1_map.mat','LPF_B1') 


end