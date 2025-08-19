% Calculate Dictionary based on the raw data

% Yun Jiang - yun.jiang@med.umich.edu
% University of Michigan
% -- All copyright reserved




function [dict,r] = Calculate_MRF_FISP_Dict(rawinfo,t1series,t2series,offres,TEtxtfile,TRtxtfile,FAtxtfile,PHtxtfile,Delays,dephase)
if nargin < 8, Delays  = [];  end;
if nargin < 7, PHtxtfile = [];end;
if nargin < 6, FAtxtfile = [];end;
if nargin < 5, TRtxtfile = [];end;
if nargin < 4, TEtxtfile = [];end;
if nargin < 3, offres = []; end;
if nargin < 2, t2series = []; end;
if nargin < 1, t1series = []; end;

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

if isempty(Delays),
    Delays = 1500; % ms
end


r = paramTable(t1series,t2series,offres); % stores T1 in first column, T2 in second column
cnt = size(r,1);

% read the basic sequence parameters from the raw data
Nex = rawinfo.Nex;
tr0 = rawinfo.TR(1)/1000/1000;
tr = importdata(TRtxtfile)*1e-6 + tr0;
tr = reshape(tr(1:Nex),[],1);

flip = importdata(FAtxtfile)*pi/180;% temporary need to remove
flip = reshape(flip(1:Nex,1),[],1);

phase = importdata(PHtxtfile)*pi/180;
phase = reshape(phase(1:Nex),[],1);

ti = rawinfo.TI(1)/1000/1000;
if isempty(ti)
    ti = 0;
end;

te0 = rawinfo.TE(1)/1000/1000;
te = importdata(TEtxtfile)*1e-6 + te0;

waitingdelay = Delays/1000/1000;
p = gcp('nocreate');
if isempty(p),
    p = gcp;
end;
poolsize = p.NumWorkers;
blocksize = ceil(cnt/poolsize);
dictr = zeros(Nex*length(te),blocksize,poolsize);
dicti = dictr;
dictz = dicti;
tic
if (1)
    parfor ipool = 1:poolsize
        r_index = (ipool-1)*blocksize+1:ipool*blocksize;
        r_index(r_index > cnt) = [];
        dr = zeros(Nex*length(te),blocksize);
        di = dr;
        dz = di;
        [dr(:,1:length(r_index)),di(:,1:length(r_index)),dz(:,1:length(r_index))] =dictionary_FISP_withRelaxDecay_mex(...
            r(r_index,:),flip,tr,te,phase,ti,dephase*pi,250,Nex,1,waitingdelay);
        dictr(:,:,ipool) = dr;
        dicti(:,:,ipool) = di;
        dictz(:,:,ipool) = dz;
    end;
else
    [dictr,dicti,dictz] =dictionary_FISP_withRelaxDecay(...
        r,flip,tr,te,phase,ti,dephase*pi,250,Nex,1,waitingdelay);
end;


TotalTime_Dictionary = toc
dict = dictr + 1i*dicti;

if 1
    dict = reshape(dict,[Nex,length(te) blocksize*poolsize]);
    dict = dict(:,:,1:cnt);
end;



return
