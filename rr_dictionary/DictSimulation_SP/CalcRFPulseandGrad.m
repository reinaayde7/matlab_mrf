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