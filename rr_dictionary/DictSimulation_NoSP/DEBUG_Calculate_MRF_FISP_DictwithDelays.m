% Calculate Dictionary based on the raw data

% Yun Jiang - yun.jiang@case.edu
% Case Western Reserve University
% -- All copyright reserved




function [dict,r] = DEBUG_Calculate_MRF_FISP_DictwithDelays(rawinfo,t1series,t2series,offres,TEtxtfile,TRtxtfile,FAtxtfile,PHtxtfile,Inv,Delays,dephase)
%if nargin < 10, repeats = [];  end;
if nargin < 9, Delays  = [];  end;
if nargin < 8, PHtxtfile = [];end;
if nargin < 7, FAtxtfile = [];end;
if nargin < 6, TRtxtfile = [];end;
if nargin < 5, TEtxtfile = [];end;
if nargin < 4, offres = [];end;
if nargin < 3, t2series = []; end;
if nargin < 2, t1series = []; end;
if nargin < 1, rawdatafilename = [];end;
if isempty(FAtxtfile),
    [tempfilename,temppathname] = uigetfile('*.txt','Select a flip angle text file');
    FAtxtfile = fullfile(temppathname,tempfilename);clear tempfilename temppathname;
end;
if isempty(PHtxtfile),
    [tempfilename,temppathname] = uigetfile('*.txt','Select a flip phase text file');
    PHtxtfile = fullfile(temppathname,tempfilename);clear tempfilename temppathname;
end;
if isempty(TEtxtfile),
    [tempfilename,temppathname] = uigetfile('*.txt','Select a TE text file');
    TEtxtfile = fullfile(temppathname,tempfilename);clear tempfilename temppathname;
end;
if isempty(TRtxtfile),
    [tempfilename,temppathname] = uigetfile('*.txt','Select a repetition text file');
    TRtxtfile = fullfile(temppathname,tempfilename);clear tempfilename temppathname;
end
if isempty(t1series),
    t1series = [60:10:2000 2020:20:3000 3050:50:3500 4000:500:5000];
end
if isempty(t2series),
    t2series = [10:5:100 105:5:200 210:10:300 350:50:500];
end
if isempty(offres),
    offres = [-100:5:100];
end
% if isempty(TE),
%     TE = 2.2; % ms
% end;

if isempty(Delays),
    Delays = 3000; % ms
end

% if isempty(repeats),
%     repeats = 3;
% end
r = paramTable(t1series,t2series,offres); % stores T1 in first column, T2 in second column
cnt = size(r,1);
% [~,~,rawinfo,~] = rawinfo;


seqtiming = diff(rawinfo.pmutime)*2.5;% in unit of 2.5 ms
% read the basic sequence parameters from the raw data
Nex = rawinfo.Nex;
tr0 = rawinfo.TR(1)/1000/1000;
disp('TR')
disp(tr0)
tr = importdata(TRtxtfile)*1e-6 + tr0;
tr = reshape(tr(1:Nex),[],1);

flip = importdata(FAtxtfile)*pi/180;% temporary need to remove
%flip = round(flip*10)/10;
%flip = flip/70.0*50; %comment out the scaling on the flip angles Yun Jiang
%08/16/2022
flip = reshape(flip(1:Nex,1),[],1);

phase = importdata(PHtxtfile)*pi/180;
phase = reshape(phase(1:Nex),[],1);

ti = rawinfo.TI(1)/1000/1000;
if isempty(ti)
    ti = 0;
end;

te0 = rawinfo.TE(1)/1000/1000;
disp('TE')
disp(te0)
te = importdata(TEtxtfile)*1e-6 + te0;

% te = [te; te(1,1)]; % append the first TE to the end, use the first TE signal for signal to reach pseudo -steady state.
% Rudy: 240716 we do not need to simulate twice Nex to have both non steady
% state and steady state. it is RAM consuming.
%te = te(1,1);
waitingdelay = Delays/1000;

%------------------------------------------------------------------------
% Rduy: GPU based speed up
%------------------------------------------------------------------------
p = gcp('nocreate');
if isempty(p),
    p = gcp;
end;
poolsize = p.NumWorkers;  %rudy debug
blocksize = ceil(cnt/poolsize);  %rudy debug
% poolsize = 1; blocksize = 1;  %rudy debug
dictr = zeros(Nex,blocksize,poolsize);
dicti = dictr;
dictz = dicti;
% % % % tic
parfor ipool = 1:poolsize %rudy debug
    r_index = (ipool-1)*blocksize+1:ipool*blocksize;
    r_index(r_index > cnt) = [];
    dr = zeros(Nex*length(te),blocksize);
    di = dr;
    dz = di;


    %Rudy: this is the one for speed 
    [dr(:,1:length(r_index)),di(:,1:length(r_index)),dz(:,1:length(r_index))] =dictionary_FISP_withRelaxDecay_mex(...
        r(r_index,:),flip,tr,te,phase,ti,dephase*pi,250,Nex,Inv,waitingdelay);

    dictr(:,:,ipool) = dr;
    dicti(:,:,ipool) = di;
    dictz(:,:,ipool) = dz;
    % disp('done3')
end
TotalTime_Dictionary = toc
dict = dictr + 1i*dicti;
dict = reshape(dict,[Nex,length(te) blocksize*poolsize]);
dict = dict(:,:,1:cnt);
dict = squeeze(dict);



%------------------------------------------------------------------------
% Rduy: for debug out of MEX
%------------------------------------------------------------------------
% [dictr,dicti,dictz] =dictionary_FISP_withRelaxDecay(...
%         r,flip,tr,te,phase,ti,dephase*pi,250,Nex,Inv,waitingdelay);
% TotalTime_Dictionary = toc
% dict = dictr + 1i*dicti;
%------------------------------------------------------------------------




clearvars dictr dicti dictz

% if 0
% dict = reshape(dict,[Nex,length(te) blocksize*poolsize]);
% dict = dict(:,:,1:cnt);
% end;
% 
% if 1  % throw away the previous repeats, and only take the last one for the pseudo steady-state signal
%     dict = dict(:,2:end,:); % throw away all the previous repeats, and the first echo signal
%     dict = circshift(dict,1,2); % shift the last echo to the first
% end

return
