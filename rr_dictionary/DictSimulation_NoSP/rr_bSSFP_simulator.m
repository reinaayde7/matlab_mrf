
% Rudy v1 - 240822


clear all, close all, clc

% RF properties
nt = 200; %RF pulse train length
flipA = 90/180*pi;
flipA_v = [10:10:150]/180*pi;

flip = repmat(flipA,nt,1); %flip angle train


phase = zeros(nt,1);
for i = 1:2:nt
    phase(i) = pi;
end

dirRF = zeros(nt,1);

%spin system properties
spins = 20;

%fat
% T1 = 187/1000; %s
% T2 = 94/1000; %s

%blood
T1 = 1122/1000; %s
T2 = 263/1000; %s

df = 0; %off resonance
df_v = [-100:1:100];

% sequence properties
te = ones(nt,1) * 1/1000; %s
tr = ones(nt,1) * 10/1000; %s

tic
for fa = 1: length(flipA_v)
    for offr = 1:length(df_v)
    flipA = flipA_v(fa);
    flip = repmat(flipA,nt,1); %flip angle train  
    [Mx(:,offr,fa), My(:,offr,fa), Mz(:,offr,fa), Mss(:,offr,fa)] = runsim(nt, te, tr, spins, flip, phase, T1, T2, df_v(offr), dirRF);
    end
end
toc

% frequency response of bSSFP at steady state
figure, 
for fa= 15: length(flipA_v)
    plot(df_v,Mss(end,:,fa)), hold on
    pause()
end
xlabel('off resonances')
ylabel('Mss (% to 1)')

%
figure, imagesc(df_v,flipA_v,squeeze(Mss(end,:,:))'), colorbar



function [Mx, My, Mz, Mss] = runsim(nt, te, tr, spins, flip, phase, T1, T2, df, dirRF)
    Mx = zeros(nt,1,'single');
    My = zeros(nt,1,'single');
    Mz = zeros(nt,1,'single');
    M = repmat([0 0 1]',1,spins);

    onemat = ones(1,spins);

    for k = 1:nt %Rudy debug
        
        % alpha pulse
        M = throtXY(flip(k),phase(k),dirRF(k))*M;
        [Ate,Bte] = freeprecess(te(k),T1,T2,df);
        Bte = Bte*onemat;
        M = Ate*M+Bte;

        %Rudy: rephase signal (ADC phase == RF phase)
        % to speed up, this is taken care outside spin simulation
        M = zrot(-phase(k))*M;

        Mx(k) = mean(M(1,:),2);
        My(k) = mean(M(2,:),2);
        Mz(k) = mean(M(3,:),2);
        Mss(k) = abs(Mx(k) + 1i*My(k));
        
        %Rudy: dephase signal (RF oscillation)
        % to speed up, this is taken care outside spin simulation
        M = zrot(phase(k))*M;

        [A,B]=freeprecess(tr(k)-te(k),T1,T2,df);
        B = B*onemat;
        M=A*M+B;
              
        
    end % timepoints
end