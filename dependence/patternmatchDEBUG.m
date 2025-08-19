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

function [t1,t2,adc,innerProduct] = patternmatchDEBUG(imadapt,roimask,r,d,dict,dictnorm,nblocks, doPlot)

    % donorm = true;
    % if ~isempty(varargin)
    %     donorm=varargin{1};
    % end
    
    dict = single(dict);
    dictnorm = single(dictnorm);

    [N,nt] = size(imadapt);
    
    
    % cnt = size(r,1);
    
    % if size(imadapt,1) ~= nt
    %     tmp=imadapt;
    %     imadapt=zeros(nt,N,N,'single');
    %     for t=1:nt
    %         imadapt(t,:,:)=squeeze(tmp(:,:,t));
    %     end
    %     clear tmp
    % end
    
    trange = 1:nt;
    %dfrange = 0;
    
    % dictnorm = sqrt(sum(dict(trange,:).*conj(dict(trange,:))));
    % if ~donorm
    %     dictnorm=ones(size(dictnorm));
    % end
    
    
    %prevlen=0;
    
        
    % idx = (iblock-1)*N*N/nblocks+1:iblock*N*N/nblocks;
    % index = idx(roimask(idx)==1);
    % Nindex = length(index);
    % if isempty(index), continue; end
    
    % xx = squeeze(imadapt(trange,index)).';
    xx = conj(imadapt);
     % xx = imadapt;

    xxNorm = sqrt(sum(xx.*conj(xx),2));
    normAll = xxNorm*dictnorm; %clear xxNorm
    

    dd = squeeze(dict);
    % dd = abs(imag(dd));
    % xx= -real(xx);
    
    % dd = conj(dd);
    % dd=dd1;
    % dd1 = -real(dd) + 1i*imag(dd);
    % dd= dd1;
    % xx = imag(xx) + 1i*real(xx);
    % xx = real(xx);

    %RR 240823 these three approaches are all equal!
    % RR 240823: for FISP is has to be conj(xx) and not conj(dd)!!!!!!
    % innerProduct = conj(xx)*dd./normAll; %clear normAll
    % RR 240823: for trueFISP is has to be xx and not conj(xx)!!!!!!
    % otherwise we loose the sign on the real part given hint on where B0
    % is leaning (+ or -)
    
    innerProduct = conj(xx)*dd./normAll; %conj run before, up here

    [value indexm] = max(abs(innerProduct),[],2); %clear innerProduct
    % [value indexm] = max(imag(innerProduct),[],2); %clear innerProduct
    % [value indexm] = max(real(innerProduct),[],2); %clear innerProduct

    % [value indexm] = max(real(innerProduct)+imag(innerProduct),[],2); %clear innerProduct
    
    %Rudy: this evaluation of inner product + looking for abs max ==
    %maximising for absolute value

    % dictCol = dict(trange,indexm);
    % indexm: has size cnt*Ndf
    %     r: size cnt x 2
    
    % indexm = indexmA;
    if doPlot
        figure, 
        subplot(311)        
        plot(abs(xx/xxNorm),'b'), hold on
        % plot(abs(squeeze(dict(:,1,indexm)))/dictnorm(indexm), 'LineWidth',2)
         plot(abs(squeeze(dict(:,indexm)))/dictnorm(indexm), 'LineWidth',2)
        ylabel('abs')
        title(['pm strength: ' num2str(value)])

        subplot(312)
        plot(real(conj(xx)/xxNorm),'b'), hold on
        plot(real(dd(:,indexm))/dictnorm(indexm), 'LineWidth',2)
        ylabel('real')

        subplot(313)
        plot(imag(conj(xx)/xxNorm),'b'), hold on
        plot(imag(dd(:,indexm))/dictnorm(indexm), 'LineWidth',2)
        xlabel('Nex'), ylabel('imag')


        % % % % % subplot(212)
        % % % % % plot(angle(xx/xxNorm),'b'), hold on
        % % % % % plot(angle(dd(:,indexm))/dictnorm(indexm), 'LineWidth',2)
        % % % % % ylabel('phase')
    end
    % coefSave = single(zeros(1,Nindex)); % M0 image
    % for iindex = 1:Nindex
    %     coefSave(iindex) = pinv(dictCol(:,iindex))*xx(iindex,:).';
    % end
    
    
    [t1t2r dfc ] = ind2sub([size(r,1), length(d)],indexm);
    t1 = r(t1t2r(:),1);
    t2 = r(t1t2r(:),2);
    adc = r(t1t2r(:),3); % Rudy: saving B0 map (offresonance map)
    % s=1;
    % adcSave = d(dfc(:));
    
    
end



% % %%
% % figure, 
% % subplot(211), plot(imag(dd(:,indexm)))
% % subplot(212), plot(real(dd(:,indexm)))