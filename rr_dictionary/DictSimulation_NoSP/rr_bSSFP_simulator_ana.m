
% Rudy v1 - 240822


clear all, close all, clc

% RF properties
nt = 400; %RF pulse train length
flipA(1) = 50/180*pi;
flipA(2) = 50/180*pi;
flip(:,1) = ones(nt,1)*flipA(1);
flip(:,2) = ones(nt,1)*flipA(2);
% flip(1,2) = flipA(2)*0.5;
flip(1:8,2) = [57, 90, 29, 7, 148, 84, 5, 180]' /180*pi;

phase = zeros(nt,1);
for i = 1:2:nt
    phase(i) = pi;
end

phase2 = phase;
phase2(1:8) = [0, pi, pi, 0, pi, pi, pi, 0];
dirRF = zeros(nt,1);
dirRF2 = dirRF;
dirRF2(1:8) = [0, 1, 0, 1, 0, 1, 0 ,0];

%spin system properties
spins = 5;
T1 = 1000/1000; %s
T2 = 52/1000; %s
df = 0; %off resonance

% sequence properties
te = ones(nt,1) * 1/1000;
tr = ones(nt,1) * 5/1000;
tr2 = tr;
tr2(1:8) = tr(1)/2;
% tr(1) = tr(1)/2;

[Mx(:,1), My(:,1), Mz(:,1), Mss(:,1)] = runsim(nt, te, tr, spins, flip(:,1), phase, T1, T2, df, dirRF);
[Mx(:,2), My(:,2), Mz(:,2), Mss(:,2)] = runsim(nt, te, tr2, spins, flip(:,2), phase2, T1, T2, df, dirRF2);

figure,
plot(Mss), grid on


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