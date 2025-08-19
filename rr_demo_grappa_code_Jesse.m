% Parameters
N = 8; % Number of coils
matrixSize = 256; % Image matrix size

phantomImg = phantom('Modified Shepp-Logan', matrixSize);
imagesc(phantomImg);

% Create coil sensitivity maps
[x, y] = meshgrid(linspace(-1,1,matrixSize), linspace(-1,1,matrixSize));
coilSensitivity = zeros(matrixSize, matrixSize, N);

for c = 1:N
    % Define coil center positions in a circular arrangement
    angle = (c-1) * (2*pi / N);
    x0 = 0.6 * cos(angle);
    y0 = 0.6 * sin(angle);
    
    % Generate Gaussian-like sensitivity profile
    coilSensitivity(:,:,c) = exp(-((x-x0).^2 + (y-y0).^2) / (0.8^2)) .* exp(1i * angle);
end

% Normalize sensitivities so the sum of squared magnitudes is ~1
sumSensitivity = sqrt(sum(abs(coilSensitivity).^2, 3));
coilSensitivity = coilSensitivity ./ max(sumSensitivity(:));

% Generate coil images
coilImages = zeros(matrixSize, matrixSize, N);
for c = 1:N
    coilImages(:,:,c) = phantomImg .* abs(coilSensitivity(:,:,c));
end

% Display results
figure;
for c = 1:N
    subplot(2, ceil(N/2), c);
    imagesc(abs(coilImages(:,:,c))); colormap gray; axis off;
    title(['Coil ', num2str(c)]);
end

ksp = transformImageToKspace(transformImageToKspace(coilImages,1),2);
% Display results
figure;
for c = 1:N
    subplot(2, ceil(N/2), c);
    imagesc(abs(ksp(matrixSize/2-30:matrixSize/2+30,matrixSize/2-30:matrixSize/2+30,c))); colormap gray; axis off;
    title(['Coil ', num2str(c)]);
end

%%
load('data_test_grappa.mat');

%%
% ksp =  transformImageToKspace(transformImageToKspace(data.image_coils,1),2);
% mask = [1:2:304];

ksp = permute(squeeze(data.kspace_spiral(:,1,:,:,1)),[2,1,3]);
mask = [1:2:48];

%%
ksp_u_plot = zeros(size(ksp));
ksp_u_plot(mask,:,:) = ksp(mask,:,:);

ksp_u = ksp(mask,:,:);
sig = permute(ksp_u,[3,1,2]);

acs = permute(ksp(24-5:24+5,:,:),[3,1,2]);

af=2;
% ws = [3];

[sigrecon,ws,src,targ]=opengrappa_no_comments(sig,acs,af); 

ksp_G = permute(sigrecon,[2,3,1]);

% ksp_u_cc = squeeze(sum( conj(squeeze(data.image_csm)).*squeeze(ksp_u_plot), 3 ));
% ksp_G_cc = squeeze(sum( conj(squeeze(data.image_csm)).*squeeze(ksp_G), 3 ));
% ksp_cc = squeeze(sum( conj(squeeze(data.image_csm)).*squeeze(ksp), 3 ));
% 
% i_u = transformKspaceToImage(transformKspaceToImage(ksp_u_cc,1),2);
% i_G = transformKspaceToImage(transformKspaceToImage(ksp_G_cc,1),2);
% i = transformKspaceToImage(transformKspaceToImage(ksp_cc,1),2);

% figure, 
% subplot(231), imagesc(abs(i))
% subplot(232), imagesc(abs(i_G))
% subplot(233), imagesc(abs(i_u))
% subplot(234), imagesc(abs(i-i))
% subplot(235), imagesc(abs(i_G-i))
% subplot(236), imagesc(abs(i_u-i))


% figure, 
% subplot(231), imagesc(abs(ksp_cc(152-20:152+20,152-20:152+20))), colorbar
% subplot(232), imagesc(abs(ksp_G_cc(152-20:152+20,152-20:152+20))), colorbar
% subplot(233), imagesc(abs(ksp_u_cc(152-20:152+20,152-20:152+20))), colorbar
% 
% subplot(234), imagesc(abs(ksp_cc(152-20:152+20,152-20:152+20)-ksp_cc(152-20:152+20,152-20:152+20))), colorbar
% subplot(235), imagesc(abs(ksp_G_cc(152-20:152+20,152-20:152+20)-ksp_cc(152-20:152+20,152-20:152+20))), colorbar
% subplot(236), imagesc(abs(ksp_u_cc(152-20:152+20,152-20:152+20)-ksp_cc(152-20:152+20,152-20:152+20))), colorbar

% figure, 
% subplot(231), imagesc(abs(ksp_cc)), colorbar
% subplot(232), imagesc(abs(ksp_G_cc)), colorbar
% subplot(233), imagesc(abs(ksp_u_cc)), colorbar
% 
% subplot(234), imagesc(abs(ksp_cc-ksp_cc)), colorbar
% subplot(235), imagesc(abs(ksp_G_cc-ksp_cc)), colorbar
% subplot(236), imagesc(abs(ksp_u_cc-ksp_cc)), colorbar


figure, 
subplot(131), imagesc(squeeze(abs(ksp_G(:,1:100,1))))
subplot(132), imagesc(squeeze(abs(ksp_u_plot(:,1:100,1))))
subplot(133), imagesc(squeeze(abs(ksp(:,1:100,1))))
%%
maskZ(abs(squeeze(ksp_u_plot(:,1,1)))>eps) = 1;
figure,imagesc(maskZ); axis square;

%%
varargin{1} = 0;%-1:1; %makes it a grappa1D on kz (work around to weights on DCF through out spiral readout 
varargin{2} = 2;
% varargin{3} = permute(ksp(152-50:152+50,152-50:152+50,:),[2,1,3]);
varargin{3} = permute(ksp(24-5:24+5,:,:),[2,1,3]);

ksp_G2 = grappa2(permute(ksp_u_plot,[2,1,3]),maskZ','idx',varargin{1}, 'width', varargin{2}, 'cal', varargin{3});  

ksp_G2 = permute(ksp_G2,[2,1,3]);
figure, 
subplot(131), imagesc(squeeze(abs(ksp_G2(:,1:100,1))))
subplot(132), imagesc(squeeze(abs(ksp_u_plot(:,1:100,1))))
subplot(133), imagesc(squeeze(abs(ksp(:,1:100,1))))

% ksp_u_cc = squeeze(sum( conj(squeeze(data.image_csm)).*squeeze(ksp_u_plot), 3 ));
% ksp_G_cc = squeeze(sum( conj(squeeze(data.image_csm)).*squeeze(permute(ksp_G2,[2, 1, 3])), 3 ));
% ksp_cc = squeeze(sum( conj(squeeze(data.image_csm)).*squeeze(ksp), 3 ));
% 
% i_u = transformKspaceToImage(transformKspaceToImage(ksp_u_cc,1),2);
% i_G = transformKspaceToImage(transformKspaceToImage(ksp_G_cc,1),2);
% i = transformKspaceToImage(transformKspaceToImage(ksp_cc,1),2);
% 
% figure, 
% subplot(231), imagesc(abs(i))
% subplot(232), imagesc(abs(i_G))
% subplot(233), imagesc(abs(i_u))
% subplot(234), imagesc(abs(i-i))
% subplot(235), imagesc(abs(i_G-i))
% subplot(236), imagesc(abs(i_u-i))
% 
% 
% figure, 
% subplot(231), imagesc(abs(ksp_cc(152-20:152+20,152-20:152+20))), colorbar
% subplot(232), imagesc(abs(ksp_G_cc(152-20:152+20,152-20:152+20))), colorbar
% subplot(233), imagesc(abs(ksp_u_cc(152-20:152+20,152-20:152+20))), colorbar
% 
% subplot(234), imagesc(abs(ksp_cc(152-20:152+20,152-20:152+20)-ksp_cc(152-20:152+20,152-20:152+20))), colorbar
% subplot(235), imagesc(abs(ksp_G_cc(152-20:152+20,152-20:152+20)-ksp_cc(152-20:152+20,152-20:152+20))), colorbar
% subplot(236), imagesc(abs(ksp_u_cc(152-20:152+20,152-20:152+20)-ksp_cc(152-20:152+20,152-20:152+20))), colorbar

%%
function [sigrecon,ws,src,targ]=opengrappa_no_comments(sig,acs,af)

    %   Please read the license text at the bottom of this program. By using this program, you 
    %   implicity agree with the license. 
    %
    %   The main points of the license:
    %
    %   1) This code is strictly for non-commercial applications. The code is
    %   protected by
    %      multiple patents.
    %   2) This code is strictly for research purposes, and should not be used in any
    %      diagnostic setting.
    %
    %   22.10.2004  Mark Griswold (mark@physik.uni-wuerzburg.de)
    
    
    
    [nc,ny,nx]=size(sig);
    [nc,nyacs,nxacs]=size(acs);    
    
    
    tic;         
    src=zeros((nyacs-af)*(nxacs-2),nc*6);
    targ=zeros((nyacs-af)*(nxacs-2),nc*(af-1));
    cnt=0;                          
    
    for xind=2:nxacs-1,  
        for yind=1:nyacs-af,        
            cnt=cnt+1;      
            src(cnt,:)=reshape(acs(:,yind:af:yind+af,xind-1:xind+1),1,nc*6);    
            % src(cnt,:)=reshape(acs(:,yind:af:yind+af,xind),1,nc*6);
            targ(cnt,:)=reshape(acs(:,yind+1:yind+af-1,xind),1,nc*(af-1));        
        end
    end
    
    
    ws=pinv(src)*targ;
    sigrecon=zeros(nc,ny*af,nx);
    sigrecon(:,1:af:end,:)=sig;  
    
    for xind=2:nx-1,    
        for yind=1:af:ny*af-af,             
            srcr=reshape(sigrecon(:,yind:af:yind+af,xind-1:xind+1),1,nc*6);    
            % srcr=reshape(sigrecon(:,yind:af:yind+af,xind),1,nc*6); 
            sigrecon(:,yind+1:yind+af-1,xind)=reshape(srcr*squeeze(ws),[nc af-1]);            
        end
    end
    
    % recon=circshift(ifft(ifft(circshift(sigrecon,[0 round(af*ny./2) round(nx./2)]),[],2),[],3),[0 round(af*ny./2) round(nx./2)]);
end