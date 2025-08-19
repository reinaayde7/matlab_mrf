pulseduration = 2000; % in unit of us
BandWidthTimeProduct = 8;
sliceThickness = 5;% in unit of mm
baseTR = 11.2*1000; % in unit of us
trtxtfile = 'FISP_TR_Body.txt';
tr = textread(trtxtfile,'%f\b');
fatxtfile = 'FISP_FA_Body.txt';
theta = textread(fatxtfile,'%f\b');
PHtxtfile = 'PH_FISP.txt';
phase = textread(PHtxtfile,'%f\b');
%theta = theta/74*69.38;
NumOfFrames = 3000;
[dict,r] = calcFISP_dict_highres_msft(pulseduration,BandWidthTimeProduct,sliceThickness,baseTR,theta,tr,phase,NumOfFrames,3000);
%save('FISP_DICT_TR10ms_msft.mat','dict','r','-v7.3');