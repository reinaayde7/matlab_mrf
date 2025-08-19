function DeltaB1 = Load_B1map(filename)%, maskbase)

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
filename = 'data/meas_MID02690_FID34521_MOLLI_T1'
    

raw = mapVBVD(filename);
raw = raw{2}; %not using noise data

disp(size(raw.image))
raw = squeeze(raw.image( :, :, :, :, :, :, :, :, :, :, :)); % Select slice in 5th pos
%disp(size(raw))


data1 = squeeze(raw( :, :, :, 1)); % 
data2 = squeeze(raw( :, :, :, 2)); % 
data3 = squeeze(raw( :, :, :, 3));
data4 = squeeze(raw( :, :, :, 4));
data5 = squeeze(raw( :, :, :, 5));
data6 = squeeze(raw( :, :, :, 6));
data7 = squeeze(raw( :, :, :, 7));
data8 = squeeze(raw( :, :, :, 8));

%data2 = squeeze(raw( :, :, :, :, 4, :, :, :, :, :, 2)); % TE 2 data
%data1 = squeeze(raw.image( :, :, :, :, :, :, :, 1, :, :, :)); % TE 1 data
%data2 = squeeze(raw.image( :, :, :, :, :, :, :, 2, :, :, :)); % TE 2 data



%Nc =size(data1,2); %Number of coils
[~,Nc,~] = size(data1)

%% Classic FFT
fprintf('FFT...\n')
image1 = zeros(size(data1,1),size(data1,3),Nc);
image2 = zeros(size(data2,1),size(data2,3),Nc);
image3 = zeros(size(data3,1),size(data3,3),Nc);
image4 = zeros(size(data4,1),size(data4,3),Nc);
image5 = zeros(size(data5,1),size(data5,3),Nc);
image6 = zeros(size(data6,1),size(data6,3),Nc);
image7 = zeros(size(data7,1),size(data7,3),Nc);
image8 = zeros(size(data8,1),size(data8,3),Nc);

for i = 1:Nc
    image1(:,:,i) = fftshift(ifft2(squeeze(data1(:,i,:))));
    image2(:,:,i) = fftshift(ifft2(squeeze(data2(:,i,:))));
    image3(:,:,i) = fftshift(ifft2(squeeze(data3(:,i,:))));
    image4(:,:,i) = fftshift(ifft2(squeeze(data4(:,i,:))));
    image5(:,:,i) = fftshift(ifft2(squeeze(data5(:,i,:))));
    image6(:,:,i) = fftshift(ifft2(squeeze(data6(:,i,:))));
    image7(:,:,i) = fftshift(ifft2(squeeze(data7(:,i,:))));
    image8(:,:,i) = fftshift(ifft2(squeeze(data8(:,i,:))));
end

%image1 = fftshift(ifft2(squeeze(data1(:,1,:))));
%disp(size(image1))
%subplot();
%imagesc(abs(squeeze(image2(:,:,1))))


%% Coil Sensitivity Map estimation
fprintf('Step 4: Estimate coil sensitivities and perform the coil combination... \n');

avg_image1 = sum(image1,3);
avg_image2 = sum(image2,3);
avg_image3 = sum(image3,3);
avg_image4 = sum(image4,3);
avg_image5 = sum(image5,3);
avg_image6 = sum(image6,3);
avg_image7 = sum(image7,3);
avg_image8 = sum(image8,3);


csm1 = estimate_csm_walsh(avg_image1);
csm2 = estimate_csm_walsh(avg_image2);
csm3 = estimate_csm_walsh(avg_image3);
csm4 = estimate_csm_walsh(avg_image4);
csm5 = estimate_csm_walsh(avg_image5);
csm6 = estimate_csm_walsh(avg_image6);
csm7 = estimate_csm_walsh(avg_image7);
csm8 = estimate_csm_walsh(avg_image8);

%subplot();
%imagesc(abs(squeeze(csm1(:,:,8))));


%% Perform coil combination
image_combined1 = single(zeros(size(data1,1),size(data1,3),Nc));
image_combined2 = single(zeros(size(data2,1),size(data2,3),Nc));
image_combined3 = single(zeros(size(data3,1),size(data3,3),Nc));
image_combined4 = single(zeros(size(data4,1),size(data4,3),Nc));
image_combined5 = single(zeros(size(data5,1),size(data5,3),Nc));
image_combined6 = single(zeros(size(data6,1),size(data6,3),Nc));
image_combined7 = single(zeros(size(data7,1),size(data7,3),Nc));
image_combined8 = single(zeros(size(data8,1),size(data8,3),Nc));


image_combined1 = sum(image1.*conj(csm1),3)./sum(csm1.*conj(csm1),3); 
image_combined2 = sum(image2.*conj(csm2),3)./sum(csm2.*conj(csm2),3); 
image_combined3 = sum(image3.*conj(csm3),3)./sum(csm3.*conj(csm3),3); 
image_combined4 = sum(image4.*conj(csm4),3)./sum(csm4.*conj(csm4),3); 
image_combined5 = sum(image5.*conj(csm5),3)./sum(csm5.*conj(csm5),3); 
image_combined6 = sum(image6.*conj(csm6),3)./sum(csm6.*conj(csm6),3); 
image_combined7 = sum(image7.*conj(csm7),3)./sum(csm7.*conj(csm7),3); 
image_combined8 = sum(image8.*conj(csm8),3)./sum(csm8.*conj(csm8),3); 

image_combined1 = image_combined1(size(image_combined1,1)/4:3*size(image_combined1,1)/4-1,size(image_combined1,2)/4:3*size(image_combined1,2)/4-1); % center image and (JESUS added -1)
image_combined2 = image_combined2(size(image_combined2,1)/4:3*size(image_combined2,1)/4-1,size(image_combined2,2)/4:3*size(image_combined2,2)/4-1); % center image and (JESUS added -1)
image_combined3 = image_combined3(size(image_combined3,1)/4:3*size(image_combined3,1)/4-1,size(image_combined3,2)/4:3*size(image_combined3,2)/4-1); % center image and (JESUS added -1)
image_combined4 = image_combined4(size(image_combined4,1)/4:3*size(image_combined4,1)/4-1,size(image_combined4,2)/4:3*size(image_combined4,2)/4-1); % center image and (JESUS added -1)
image_combined5 = image_combined5(size(image_combined5,1)/4:3*size(image_combined5,1)/4-1,size(image_combined5,2)/4:3*size(image_combined5,2)/4-1); % center image and (JESUS added -1)
image_combined6 = image_combined6(size(image_combined6,1)/4:3*size(image_combined6,1)/4-1,size(image_combined6,2)/4:3*size(image_combined6,2)/4-1); % center image and (JESUS added -1)
image_combined7 = image_combined7(size(image_combined7,1)/4:3*size(image_combined7,1)/4-1,size(image_combined7,2)/4:3*size(image_combined7,2)/4-1); % center image and (JESUS added -1)
image_combined8 = image_combined8(size(image_combined8,1)/4:3*size(image_combined8,1)/4-1,size(image_combined8,2)/4:3*size(image_combined8,2)/4-1); % center image and (JESUS added -1)

DataCube = zeros(size(image_combined1,1),size(image_combined1,2),8);
T1map = zeros(size(image_combined1,1),size(image_combined1,2));
DataCube(:,:,1) = image_combined1;
DataCube(:,:,2) = image_combined2;
DataCube(:,:,3) = image_combined3;
DataCube(:,:,4) = image_combined4;
DataCube(:,:,5) = image_combined5;
DataCube(:,:,6) = image_combined6;
DataCube(:,:,7) = image_combined7;
DataCube(:,:,8) = image_combined8;

x = [1: 1: 8]; %Put the TI's here!
y = abs(squeeze(DataCube(33,33,:)));
%disp(size(y(:)))
%scatter(x,y)

subplot()
imagesc(abs(image_combined8));
colormap('gray');
daspect([1 1 1])
colorbar

pause(99999999999999999999999)
%% Flip the map to match the reconstruction orientation
%OffRes = rot90(OffRes,1);

%{

disp(size(Gx))
subplot()
%imagesc(squeeze(Gx(:,11)));
P = polyfit(linspace(0,size(Gx,1),size(Gx,1)),Gx_thres(:,33,1),1)
yfit = polyval(P,linspace(0,size(Gx,1),size(Gx,1)));
plot(1:64,Gx(:,33,1),1:64,yfit, 'r-.')
%plot(1:64,Gx(:,33,1),1:64,Gx_thres(:,33,1), 'o')
title('G');
pause(9999999999999999999)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DERIVATIVE REGULARIZATION ENDS   %%%%%%%
%}

%{
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
%}


%% Load through-time MRF average
%  to create the mask
%load("Blurred_Image.mat","image_combined") %in folder 2D_MRF_Reconstruction_main

%C = imquantize(abs(mean(image_combined(:,:,1),3)),300); %mask from MRF through-time average THIS FILE IS IN THE 2D_MRF_ReconstructionMain folder

% get some 2D matrix, and plot as surface
%subplot(121), imagesc(C);

%{
disp(size(A))
subplot()
imagesc(squeeze(C(:,:,1)));
title('Off res map')
colorbar
pause(99999999)
%}

%% Resize Off-resonance map
%# create interpolant
A = OffRes;
[X,Y] = meshgrid(1:size(A,2), 1:size(A,1));
%[XX,YY] = meshgrid(1:size(I,2), 1:size(I,1));
F = scatteredInterpolant(X(:), Y(:), A(:), 'linear');
%FF = scatteredInterpolant(XX(:), YY(:), I(:), 'linear');

%# interpolate over a finer grid
FOV = 400;
[U,V] = meshgrid(linspace(1,size(A,2),FOV), linspace(1,size(A,1),FOV));
%[UU,VV] = meshgrid(linspace(1,size(I,2),FOV), linspace(1,size(I,1),FOV));
%subplot(122), imagesc(F(U,V));

InterpolatedField = F(U,V);
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
%}
LPF = Lowpass(InterpolatedField, 15);
%InterpolatedDicom = rot90(FF(UU,VV),3); %This is for the DICOM map to be in the same direction of the acquired map (single slice)
%InterpolatedDicom = FF(UU,VV); %This is for the DICOM map to be in the same direction of the acquired map (single slice)
Thres = 220; %Cropping off-res magnitude in +/- Thres
InterpolatedField(InterpolatedField>Thres)=Thres;
InterpolatedField(InterpolatedField<-Thres)=-Thres;
LPF(LPF>Thres)=Thres;
LPF(LPF<-Thres)=-Thres;
%InterpolatedDicom(InterpolatedDicom>170)=170;
%InterpolatedDicom(InterpolatedDicom<-170)=-170;
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

%% Plot the map to compare orientation
%mr_imshow(abs((LPF)));
%{
cmlist = [-Thres, Thres];
subplot()
imagesc(rot90(((LPF(:,:))),1),cmlist);
colormap('grey');
title('Off-resonance (Hz)')
colorbar
daspect([1 1 1]) 
%}

cmlist = [-Thres, Thres];
figure('name','T1 T2 and Proton Density');
subplot(1,2,1)
imagesc(rot90(((InterpolatedField(:,:))),1),cmlist);
colormap('gray');
title('Calculated map')
colorbar

subplot(1,2,2)
imagesc(rot90(((LPF(:,:))),1),cmlist);
colormap('gray');
title('Low Pass Filtered Field Map')
colorbar

%{

subplot(1,4,1)
imagesc(squeeze(A(:,:,1)));
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

