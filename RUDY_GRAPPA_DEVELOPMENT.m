

cal = permute(raw_noG(:,:,15:35,[1:48:end]),[1,3,2,4]); %13:37
data_u = permute(raw(:,:,:,1),[1,3,2]);
data = permute(raw_noG(:,:,:,1),[1,3,2]);

[nro, npar, ncoils] = size(data_u);
%%
% clearvars src trg
pad=2; % [-pad 0 +pad] calibration indexes
kx = 3;
ky=2;

cal_padded = zeros(size(cal,1)+pad*2, size(cal,2), size(cal,3), size(cal,4));
cal_padded(pad+1:end-pad,:,:,:) = cal;
data_u_padded = zeros(size(data_u,1)+pad*2, size(data_u,2)+1, size(data_u,3));
data_u_padded(pad+1:end-pad,1:end-1,:,:) = data_u;
%%
j=1;
f = waitbar(0, 'GRAPPA recon...'); tic;
for ro=pad+1:nro %ro: readout point to recon
    waitbar(ro/(nro-pad+1), f, sprintf('RO: %d/%d', ro, nro-pad+1));
    cal_i = cal_padded(ro-pad:ro+pad,:,:,:); %calibration [-pad 0 +pad] area around ro point
    
    % calibration block around RO point
    k=1;
    for r = ceil(kx/2):pad*2+1-floor(kx/2) %indexes in cal_i to be explored with a given kx kernel
        for p =ky:size(cal_i,2)-ky+1 %indexes in cal_i to be recon with a given ky kernel
            if ky == 2
                p_idx = [-1,1];
            elseif ky == 4
                p_idx = [-3,-1,1,3];
            else %ky=6;
                p_idx = [-5,-3,-1,1,3,5];
            end

            src(:,:,k) = reshape(cal_i(r-floor(kx/2):r+floor(kx/2),p+p_idx,:,:),kx*ky*ncoils, size(cal,4)); %source points
            trg(:,:,k) = reshape(cal_i(r,p,:,:),1*ncoils,size(cal,4)); %target points
            k=k+1;
        end
    end

    src = reshape(src,kx*ky*ncoils, size(src,2)*size(src,3));
    trg = reshape(trg,1*ncoils, size(trg,2)*size(trg,3));
    ws = pinv(src.')*trg.';
    
    if j<2
        fprintf('Linear system: #of Eq. %d ; #of Var. %d ; Ratio of Overdet: %.3f \n',size(src,2), size(src,1), size(src,2) / size(src,1));
    end
    clearvars src trg

    % reconstruction: RO point all partitions (NB: accounts only for accel
    % of 2 now
    for p =2:2:npar
        srcr = reshape(data_u_padded(ro-floor(kx/2):ro+floor(kx/2),p+p_idx,:),kx*ky*ncoils,1);
        trgr = (srcr.'*ws).'; % target points
        data_u_padded(ro,p,:) = reshape(trgr,[1,ncoils]);
    end
    j=j+1;
end
close(f); time = toc; GRAPPA_time = time2clock(time);
%%
data_r = data_u_padded(pad+1:end-pad,1:end-1,:,:);
a(:,1) = squeeze(data(5,:,1));
a(:,2) = squeeze(data_u(5,:,1));

figure, plot(abs(a(:,1))), hold on 
plot(abs(a(:,2)))
%%
 % idx = 1:50;
 % idx = 50:150;
% idx = 150:300;
% idx = 500:1000;
% idx = 1000:2000;
% idx = 2000:3000;

coil = 3;

figure, 
subplot(222),imagesc(angle(squeeze(data_r(idx,:,coil)))), title('grappa'), colorbar
subplot(221),imagesc(angle(squeeze(data(idx,:,coil)))), title('orig'), colorbar
subplot(223),imagesc(angle(squeeze(data_u(idx,:,coil)))), title('undersmp'), colorbar
subplot(224),imagesc( angle( squeeze( data_r( idx,:,coil ) ) - squeeze( data(idx,:,coil) )) ), title('diff'), colorbar

figure, 
subplot(222),imagesc(abs(squeeze(data_r(idx,:,coil)))), title('grappa'), colorbar
subplot(221),imagesc(abs(squeeze(data(idx,:,coil)))), title('orig'), colorbar
subplot(223),imagesc(abs(squeeze(data_u(idx,:,coil)))), title('undersmp'), colorbar
subplot(224),imagesc( abs( squeeze( data_r( idx,:,coil ) ) - squeeze( data(idx,:,coil) )) ), title('diff'), colorbar

%%
error = abs( squeeze( data_r( idx,:,coil ) ) - squeeze( data(idx,:,coil) ));
error = sum(sum(error(:,2:2:10),1),2);

%% aliasing check
data_u(:,2:2:end,:) = 0;
i_d = fftshift(fft(data_u,[],2),2);
figure, imagesc(squeeze(abs(i_d(1:100,:,1))))

i_d = fftshift(fft(data,[],2),2);
figure, imagesc(squeeze(abs(i_d(1:100,:,1))))

i_d = fftshift(fft(data_r,[],2),2);
figure, imagesc(squeeze(abs(i_d(1:100,:,1))))








%% function testout
% nb: HERE WE ARE NOT COPYING BACK THE CENTER OF KSPACE TO THE GRAPPA
% RESULTS
cal = permute(raw_noG(:,:,15:35,[1:48:end]),[1,3,2,4]); %13:37
data = permute(raw_noG(:,:,:,49),[1,3,2]);

data_u = data;
data_u(:,2:2:end,:) = 0;
% data_u = permute(raw(:,:,:,1),[1,3,2]);

pad=2;
kx = 3;
ky=2;
% 
% tic
% [data_r, ws1] = RR_spiralGRAPPA_MRF_parundersmp(data_u,pad,kx,ky,1,cal,[]);
% toc

tic
[data_r, ~] = RR_spiralGRAPPA_MRF_parundersmp(data_u,pad,kx,ky,0,cal,ws1);
toc

%% test
% nb: HERE WE ARE NOT COPYING BACK THE CENTER OF KSPACE TO THE GRAPPA
% RESULTS

coil = 10;
figure,
subplot(231), imagesc(abs(squeeze(data_u(50:150,:,coil)))), title('data_u','Interpreter','none')
subplot(232), imagesc(abs(squeeze(data_r(50:150,:,coil)))), title('data_r','Interpreter','none')
subplot(233), imagesc(abs(squeeze(data(50:150,:,coil)))), title('data','Interpreter','none')


i_data_u = fftshift(fft(data_u,[],2),2);
i_data = fftshift(fft(data,[],2),2);
i_data_r = fftshift(fft(data_r,[],2),2);

subplot(234), imagesc(abs(squeeze(i_data_u(1:150,:,coil)))), title('i_data_u','Interpreter','none')
subplot(235), imagesc(abs(squeeze(i_data_r(1:150,:,coil)))), title('i_data_r','Interpreter','none')
subplot(236), imagesc(abs(squeeze(i_data(1:150,:,coil)))), title('i_data','Interpreter','none')





%% iterative over tp testaout
% nb: HERE WE ARE NOT COPYING BACK THE CENTER OF KSPACE TO THE GRAPPA
% RESULTS
tic

cal_idx = 1; %takes care we keep taking all indexes correpondent to same spiral arm for cal
for tp = 1:50%600

    data = permute(raw_noG(:,:,:,tp),[1,3,2]);    
    data_u = data;
    data_u(:,2:2:end,:) = 0;
    % data_u = permute(raw(:,:,:,1),[1,3,2]);

    pad=2;
    kx = 3;
    ky=2;
    
    if tp<48+1
        doCal=1;
        wsc=[];  
        v=cal_idx:48:600
        cal = permute(raw_noG(:,:,15:35,[cal_idx:48:end]),[1,3,2,4]); %13:37     
        [data_r, ws{cal_idx}] = RR_spiralGRAPPA_MRF_parundersmp(data_u,pad,kx,ky,doCal,cal,wsc);
        cal_idx=cal_idx+1;
    else %calibration already done!
        doCal=0;
        [data_r, ~] = RR_spiralGRAPPA_MRF_parundersmp(data_u,pad,kx,ky,doCal,cal,ws{cal_idx});
        cal_idx=cal_idx+1;
    end

    if cal_idx>48
        cal_idx=1;
    end    
    
    raw_r(:,:,:,tp) = data_r;
   
end
toc
raw_r = permute(raw_r,[1,3,2,4]);
%% test
% nb: HERE WE ARE NOT COPYING BACK THE CENTER OF KSPACE TO THE GRAPPA
% RESULTS
coil = 10;
tp = 50;

raw_u = raw_noG;
raw_u(:,:,2:2:end,:) = 0;

figure,
subplot(231), imagesc(abs(squeeze(raw_u(50:150,coil,:,tp)))), title('data_u','Interpreter','none')
subplot(232), imagesc(abs(squeeze(raw_r(50:150,coil,:,tp)))), title('data_r','Interpreter','none')
subplot(233), imagesc(abs(squeeze(raw_noG(50:150,coil,:,tp)))), title('data','Interpreter','none')

i_raw_u = fftshift(fft(raw_u,[],3),3);
i_raw_noG = fftshift(fft(raw_noG,[],3),3);
i_raw_r = fftshift(fft(raw_r,[],3),3);

subplot(234), imagesc(abs(squeeze(i_raw_u(1:150,coil,:,tp)))), title('i_data_u','Interpreter','none')
subplot(235), imagesc(abs(squeeze(i_raw_r(1:150,coil,:,tp)))), title('i_data_r','Interpreter','none')
subplot(236), imagesc(abs(squeeze(i_raw_noG(1:150,coil,:,tp)))), title('i_data','Interpreter','none')


%%
