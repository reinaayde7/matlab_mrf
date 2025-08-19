% Calcuate the MRF-FISP dictionary
% The function sets up the slice selective gradient, and the RF pulses
% call bloch.c to calcuate the signal
% Yun Jiang --- yun.jiang@case.edu  07/26/15
function [dict,r] = calcFISP_dict_highres(pulseduration,BandWidthTimeProduct,sliceThickness,baseTR,theta,tr,NumOfFrames,InvPerFrames)
if nargin < 6,
    trtxtfile = 'TR_FISP.txt';
    tr = textread(trtxtfile,'%f\b');
end;

if nargin < 5,
    fatxtfile = 'FA_FISP.txt';
    theta = textread(fatxtfile,'%f\b');
end;

if nargin < 4, baseTR = 12000;% % in unit of us
end;

if nargin < 3, sliceThickness = 5; % in unit of mm
end;

if nargin < 2, BandWidthTimeProduct = 8;
end;

if nargin < 1, pulseduration = 2000; % in unit of us
end;
% Fixed parameter
gamma = 42.58; % in unit of MHz/T, 4258 Hz/G
pulseNumber = 200; % This can be adjusted based on RF pulse in the sequence.
% pulseduration/pulseNumber is the duration per step in the RF pulse
sliceThickFactor = 2;
NumOfSpins = 200; % Number of spins across the (sliceThickness * sliceThickFactor)

TI = 20640-10240-800;% in unit of us...based on the event block in Siemens IDEA.
%TI = 12700 -10240 - 800;
%The gap between the end of adiabetic pulse to the slice seletive gradient


%% Calculate the sinc pulse and slice selective gradient
[SincPulse,rf,SliceSelGrad,t,dT] = CalcRFPulseandGrad(pulseduration,pulseNumber,BandWidthTimeProduct,sliceThickness);

%% T1 T2 series
%t1series = 1000;%[10:10:90 100:20:1000 1040:40:2000 2050:100:3000 3200:200:5000];
%t2series = 100;%[2:2:8 10:5:100 100:10:300 350:50:800];% 900:100:2000];

%t1series = 1000;%[500:200:1000];%[10:10:90 100:20:1000 1040:40:2000 2050:100:3000 3100:200:4500];
%t2series = [50:20:300];%[2:2:98 100:5:150 160:10:300 320:20:500 550:50:1000];% 1100:100:1600 1800:200:3000];
% t1 = zeros(1000,1);
% t1(1,1) = 10;
% for tt = 2:1000,
%     t1(tt,1) = t1(tt-1,1)*1.02;
% end
% t1series = t1(find(t1<3000));
% t2 = zeros(1000,1);
% t2(1,1) = 2;
% for tt = 2:1000,
%     t2(tt,1) = t2(tt-1,1)*1.02;
% end
% t2series = t2(find(t2<500));
t1series = [10:10:2000 2020:20:3000 3050:50:3500 4000:500:5000];
t2series = [2:2:100 105:5:300 310:10:500 520:20:800];% 850:50:1200];% 310:10:500 520:20:800 850:50:1500 1600:100:2000];

cnt=0;
for it1 = 1:length(t1series)
    for it2 = 1:length(t2series)
        if (t1series(it1)>=t2series(it2))
            cnt=cnt+1;
            r(cnt,1)=t1series(it1);
            r(cnt,2)=t2series(it2);
        end
    end
end






%% Create the Dictionary
T1 = r(:,1);
T2 = r(:,2);

dict = zeros(NumOfFrames,cnt);
tic
parfor jj = 1:cnt,
    t1 = T1(jj);
    t2 = T2(jj);
    mx = zeros(NumOfSpins,1);
    my = zeros(NumOfSpins,1);
    mz = ones(NumOfSpins,1);% from -1 to 0.95 to account the inversion efficiency
    
    x = [-bitshift(NumOfSpins,-1):bitshift(NumOfSpins,-1)-1]*((sliceThickness)*sliceThickFactor)/(NumOfSpins);
    %[mx,my,mz] = bloch(0,0,TI/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
    for ii = 1:NumOfFrames,
        if (mod(ii,InvPerFrames) == 1)
            mz = (-0.95)*mz;
            mx = zeros(NumOfSpins,1);
            my = zeros(NumOfSpins,1);
            [mx,my,mz] = bloch(0,0,TI/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
        end;
        theta1 = theta(ii)/180*pi;       
        F1 = sqrt((sum(real(SincPulse),2)).^2 + (sum(imag(SincPulse),2)).^2);
        
        A1 = (theta1./((gamma*100)*2*pi*F1*(dT/(1000*1000))));%T->gauss
        rf_scaled = (A1*rf);
        
        grad = SliceSelGrad;
        
        RepTime = baseTR + tr(ii);
        TE = 2180;%in us
        
        [mx,my,mz] = bloch(rf_scaled,grad/10,(t(2)-t(1))/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
        
        [mx,my,mz] = bloch(0,0,(TE-t(end))/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
        msig=mx+1i*my;
        [mx,my,mz] = bloch(0,0,(RepTime - TE)/1000/1000,t1/1000,t2/1000,0,x/10,0,mx,my,mz);
    
        dict(ii,jj) = mean(msig);

    
    end
       
    
end
toc

function [SincPulse,rf,SliceSelGrad,t,dT] = CalcRFPulseandGrad(pulseduration,pulseNumber,BandWidthTimeProduct,sliceThickness)
gamma = 42.58; % in unit of MHz/T, 4258 Hz/G

%% Setting parameters and timings
% From the beginning of RF pulse to the beginning the ADC. The ramptime is
% ignored. The rewind gradient can be freedomly distributed in the
% allowed time that is after the slice selective gradient and
% before the acquisition.

% Here pulseduration + pulseduration*rewindgradientratio
rewindgradient_ratio = 0.1;
Duration_SliGrad = (1 + rewindgradient_ratio)*pulseduration; % in unit of us


dT = pulseduration/pulseNumber; % in unit of us
t = (1:dT:Duration_SliGrad); % in unit of us
rf = zeros(size(t));

tpulse = [-bitshift(pulseNumber,-1):bitshift(pulseNumber,-1)-1]*dT;

%% Construct the Hanning filtered Sinc pulse
% Hanning filtered Sinc Pulse  -- based on Siemens Education Materials
% on RealTimeEvents p23
Gs = BandWidthTimeProduct*1000*1000/(gamma*pulseduration*sliceThickness);%mT/m
tempT = tpulse;%[-bitshift(pulseNumber,-1):bitshift(pulseNumber,-1)-1].*(pulseduration)/pulseNumber;%us
tempH = (1+cos((2*pi.*tempT/(pulseduration))))/2;
Emp = 1;% What is EMP??-- default Empirical factor is 1 that is the default setting in the pulse sequence.
tempS = (pulseduration/BandWidthTimeProduct).*sin(Emp*pi*gamma*Gs*sliceThickness/(1000*1000).*tempT)./(pi.*tempT);
indx = find(tempT==0);
tempS(1,indx) = 1;
SincPulse = tempH.*tempS;

rf(1,1:length(SincPulse)) = SincPulse;

%% Set Slice Gradient Strength
SliceSelGradAmp = Gs; %mT/m

tmpSliceSelGrad(1,1:pulseNumber) = SliceSelGradAmp;
tmpSliceRewGrad(1,1:pulseNumber*rewindgradient_ratio) = -(SliceSelGradAmp*length(1:pulseNumber))/2/...
    length(pulseNumber+1:pulseNumber*(1+rewindgradient_ratio));
SliceSelGrad = cat(2,tmpSliceSelGrad,tmpSliceRewGrad);

