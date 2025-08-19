% PATTERNMATCH Performs dot product pattern matching between MRF images and
% the dictionary. Returns maps of T1, T2, M0, off-resonance, and B1. Also
% returns a map showing the inner product value of the best-matching
% dictionary atom at each pixel location.
% 
% Inputs:
%   imadapt: (N x N x timepoints) MRF images
%   roimask: (NxN, or empty) ROI defining the object support. Pixels
%               outside this region will not be matched. You can also pass
%               in an emptry matrix, and then it will match all the pixels 
%               in the image.
%   r: (parameters x 4) list of parameter combinations
%               column 1: T1 (ms)
%               column 2: T2 (ms)
%               column 3: off-resonance (Hz)
%               column 4: B1 (scaling factor greater than 0)
%          If you aren't fitting df or B1, you can pass in a matrix with
%          only 2 columns
%   dict: (timepoints x parameters) the MRF dictionary
%   image_blocks: To speed up computation time, divide the image into this
%           many blocks and run pattern matching in parallel.
%               
% Outputs
%       t1big: (N x N) T1 map
%       t2big: (N x N) T2 map
%       m0big: (N x N) M0 map
%       dfbig: (N x N) off-resonance map
%       b1big: (N x N) B1 map
%       scores: (N x N) inner product value at each pixel with best matching
%               dictionary atom
% 
% Based on Jesse Hamilton's code
% 
% Added additional adc map, but remove B1 map, and off-resonance map
% Yun Jiang (yunjiang@med.umich.edu)
% University of Michigan

function [t1big,t2big,adcbig,m0big] = patternmatch(imadapt,roimask,r,d,dict,nblocks)

% donorm = true;
% if ~isempty(varargin)
%     donorm=varargin{1};
% end

[N,~,nt] = size(imadapt);

t1big = zeros(N,N);
t2big = zeros(N,N);
m0big = zeros(N,N);
adcbig = zeros(N,N);

% cnt = size(r,1);

if size(imadapt,1) ~= nt
    tmp=imadapt;
    imadapt=zeros(nt,N,N,'single');
    for t=1:nt
        imadapt(t,:,:)=squeeze(tmp(:,:,t));
    end
    clear tmp
end

trange = 1:nt;
%dfrange = 0;

dictnorm = sqrt(sum(dict(trange,:).*conj(dict(trange,:))));
% if ~donorm
%     dictnorm=ones(size(dictnorm));
% end


%prevlen=0;
hw = waitbar(0,'Pattern matching: Please wait');
for iblock=1:nblocks
        waitbar(iblock/nblocks,hw,sprintf('Pattern Matching - wait... (%d/%d)', iblock, nblocks));
    %     fprintf('Pattern recognition: %d\n',iblock);
    
%     msg = sprintf('Pattern recognition: %d/%d\n',iblock,nblocks);
%     fprintf([repmat('\b',1,prevlen) '%s'],msg);
%     prevlen = numel(msg);
    
    idx = (iblock-1)*N*N/nblocks+1:iblock*N*N/nblocks;
    index = idx(roimask(idx)==1);
    Nindex = length(index);
    if isempty(index), continue; end
    
    xx = squeeze(imadapt(trange,index)).';
    xxNorm = sqrt(sum(xx.*conj(xx),2));
    normAll = xxNorm*dictnorm; clear xxNorm

    f = exp(-1i*angle(xx));
    xx = xx.*f;
    dd = dict(trange,:).*f;
    innerProduct = conj(xx)*dd./normAll; clear normAll
    % innerProduct = conj(xx)*dict(trange,:)./normAll; clear normAll
    % innerProduct = xx*dict(trange,:)./normAll; clear normAll 
    
    [value indexm] = max(abs(innerProduct),[],2); clear innerProduct
    % [value indexm] = max(imag(innerProduct),[],2); clear innerProduct
    dictCol = dict(trange,indexm);
    % indexm: has size cnt*Ndf
    %     r: size cnt x 2
    
    coefSave = single(zeros(1,Nindex)); % M0 image
    for iindex = 1:Nindex
        coefSave(iindex) = pinv(dictCol(:,iindex))*xx(iindex,:).';
    end
    
 
    [t1t2r dfc ] = ind2sub([size(r,1), length(d)],indexm);
    t1Save = r(t1t2r(:),1);
    t2Save = r(t1t2r(:),2);
    adcSave = r(t1t2r(:),3); % Rudy: saving B0 map (offresonance map)
    % adcSave = d(dfc(:));
    
    t1big(index) = t1Save;
    t2big(index) = t2Save;
    m0big(index) = coefSave;
    adcbig(index) = adcSave;
    
end

close(hw)