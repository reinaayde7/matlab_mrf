%Script to calculate 2D k-space from measurements as in Duyn et al.  JMR 132, p.150-153 (1998)
%Averages across k-space measurements to reduce noise in the measurement.
%Combines data across multiple coils to improve SNR
%
%Assumes there is only a single Siemens meas.out or .dat file in the directory
%
%The raw data must have the following shape
% [COLUMNS,COILS,SHOTS,AVERAGES,1,1,1,1,REPETITIONS] -- VB data structure
% [COLUMNS,COILS,SHOTS,1,1,AVERAGES,1,1,REPETITIONS] -- VD data structure

%
% Repetitions are assumed to be ordered as follows
% 1.) Slice at +y, Y GRAD ON
% 2.) Slice at -y, Y GRAD ON
% 3.) Slice at -x, X GRAD ON
% 4.) Slice at +x, X GRAD ON
% 5.) Slice at +y, Y GRAD OFF
% 6.) Slice at -y, Y GRAD OFF
% 7.) Slice at -x, X GRAD OFF

% 8.) Slice at +x, X GRAD OFF
% For resolutions below 2mm or so, the above 8 repetitions should be sufficient for the measurement
%
%
% For higher resolution cases, additional measurements with shifted k-spaces should be performed
% as described in Beaumont et al. MRM 58, p.200-205 (2007)
%
% e.g.
%
% 9.) Slice at +y, Y GRAD ON,  SHIFTED IN Y
% 10.) Slice at -y, Y GRAD ON,  SHIFTED IN Y
% 11.) Slice at -x, X GRAD ON,  SHIFTED IN X
% 12.) Slice at +x, X GRAD ON,  SHIFTED IN X
% 13.) Slice at +y, Y GRAD OFF,  SHIFTED IN Y
% 14.) Slice at -y, Y GRAD OFF,  SHIFTED IN Y
% 15.) Slice at -x, X GRAD OFF,  SHIFTED IN X
% 16.) Slice at +x, X GRAD OFF,  SHIFTED IN X
% 17.) Slice at +y, Y GRAD ON,  SHIFTED IN -Y
% 18.) Slice at -y, Y GRAD ON,  SHIFTED IN -Y
% 19.) Slice at -x, X GRAD ON,  SHIFTED IN -X
% 20.) Slice at +x, X GRAD ON,  SHIFTED IN -X
% 21.) Slice at +y, Y GRAD OFF,  SHIFTED IN -Y
% 22.) Slice at -y, Y GRAD OFF,  SHIFTED IN -Y
% 23.) Slice at -x, X GRAD OFF,  SHIFTED IN -X
% 24.) Slice at +x, X GRAD OFF,  SHIFTED IN -X

% Modify from previous script
% v0.1 -- turn the script into a function
% v0.2 -- save as an binary file for Siemens C2P package --not critical
% v0.3 -- set uplimit as an output based the kmax
% v0.4 -- set the b0 eddy current measurement as outputs -- Yun Jiang
%         09/06/2016

% Yun Jiang -- yun.jiang@case.edu 
function [kxall,kyall,uplimit,kxall_b0,kyall_b0] = GradTrajMeas(filename,saveBinary)

% I probably should read this from the header instead of hard-coding!
loc=3;  %location of slice from isocenter (cm)

doPLOTS=1;  %Set doPLOTS to 0 to suppress plots during trajectory calculation

nshifts=1  %number of k-space shifts done in the pulse sequence

%Read in the raw data and save to rawdata.mat
if nargin < 2, saveBinary=[];end;
if isempty(saveBinary),saveBinary = 0;end; % default is not saving as the binary
if nargin < 1, filename=[];end;
if isempty(filename),
    p = pwd;
    %     while 1,
    %         datapath = fullfile(BASEPATH,'data');
    %         if exist(datapath,'dir'), cd(datapath); break; end
    %         datapath = fullfile(BASEPATH,'test','data');
    %         if exist(datapath,'dir'), cd(datapath); break; end
    %         break;
    %     end; clear datapath;
    [tmpfile, tmppath] = uigetfile('*.dat', 'Pick a dataset file to load');
    if isequal(tmpfile,0) | isequal(tmppath,0), cd(p); clear p; return; end
    filename = fullfile(tmppath,tmpfile); %clear tmpfile tmppath;
    %cd(p); clear p;
end

% [raw noise ref phasecore centerlines header]=read_meas_vb13(filename,'STD',0);

%[raw,info] = loadSiemensRawData(filename);
datastruct = mapVBVD(filename);
if (length(datastruct) == 1)
    raw = (datastruct.image());
else
    raw = (datastruct{1,2}.image());
end
%raw = raw(1:4000,:,:,:,:,:,:,:,:);
%which averages to use
%navgs=1:size(raw,4);  %Yun Jiang VB - average on the 4th

navgs=1:size(raw,6);  %Yun Jiang VD - average on the 6th

%n3D=size(raw,6);
ncoils = size(raw,2);

for proj=1:size(raw,3)  %interleaves (shots) loop
    
    %preallocate
    posx=zeros(size(raw,1),ncoils); posx0=posx; negx=posx; negx0=posx;
    posy=zeros(size(raw,1),ncoils); posy0=posy; negy=posy; negy0=posy;
    posz=zeros(size(raw,1),ncoils); posz0=posz; negz=posz; negz0=posz;
    
    for coil=1:ncoils  %COIL loop
        
        for shifts=1:nshifts  %k-space SHIFTS loop
            %shifts = shifts + 2;
            %if(shifts==1 || shifts==2 || shifts==3)
            %tmp=mean(raw(:,coil,proj,navgs,1,1,1,1,1+8*(shifts-1)),4)./nshifts;
            tmp=mean(raw(:,coil,proj,1,1,navgs,1,1+8*(shifts-1)),6)./nshifts;
            posy(:,coil)=posy(:,coil)+tmp.*exp(-1i*angle(tmp(1)));  %5th argument is the slice index, 8th argument is rep

            %tmp=mean(raw(:,coil,proj,navgs,1,1,1,1,5+8*(shifts-1)),4)./nshifts;
            tmp=mean(raw(:,coil,proj,1,1,navgs,1,1,5+8*(shifts-1)),6)./nshifts;
            posy0(:,coil)=posy0(:,coil)+tmp.*exp(-1i*angle(tmp(1)));
            
            %tmp=mean(raw(:,coil,proj,navgs,1,1,1,1,2+8*(shifts-1)),4)./nshifts;
            tmp=mean(raw(:,coil,proj,1,1,navgs,1,1,2+8*(shifts-1)),6)./nshifts;
            negy(:,coil)=negy(:,coil)+tmp.*exp(-1i*angle(tmp(1)));
            
            %tmp=mean(raw(:,coil,proj,navgs,1,1,1,1,6+8*(shifts-1)),4)./nshifts;
            tmp=mean(raw(:,coil,proj,1,1,navgs,1,1,6+8*(shifts-1)),6)./nshifts;
            negy0(:,coil)=negy0(:,coil)+tmp.*exp(-1i*angle(tmp(1)));
            
            %tmp=mean(raw(:,coil,proj,navgs,1,1,1,1,4+8*(shifts-1)),4)./nshifts;
            tmp=mean(raw(:,coil,proj,1,1,navgs,1,1,4+8*(shifts-1)),6)./nshifts;
            negx(:,coil)=negx(:,coil)+tmp.*exp(-1i*angle(tmp(1)));
            
            %tmp=mean(raw(:,coil,proj,navgs,1,1,1,1,8+8*(shifts-1)),4)./nshifts;
            tmp=mean(raw(:,coil,proj,1,1,navgs,1,1,8+8*(shifts-1)),6)./nshifts;
            negx0(:,coil)=negx0(:,coil)+tmp.*exp(-1i*angle(tmp(1)));
            
            %tmp=mean(raw(:,coil,proj,navgs,1,1,1,1,3+8*(shifts-1)),4)./nshifts;
            tmp=mean(raw(:,coil,proj,1,1,navgs,1,1,3+8*(shifts-1)),6)./nshifts;
            posx(:,coil)=posx(:,coil)+tmp.*exp(-1i*angle(tmp(1)));
            
            %tmp=mean(raw(:,coil,proj,navgs,1,1,1,1,7+8*(shifts-1)),4)./nshifts;
            tmp=mean(raw(:,coil,proj,1,1,navgs,1,1,7+8*(shifts-1)),6)./nshifts;
            posx0(:,coil)=posx0(:,coil)+tmp.*exp(-1i*angle(tmp(1)));
            
        end
        if(coil>1) %force consistent phase accross coils
            posx(:,coil)=posx(:,coil).*exp(-1i*angle(posx(1,1)));
            posx0(:,coil)=posx0(:,coil).*exp(-1i*angle(posx0(1,1)));
            posy(:,coil)=posy(:,coil).*exp(-1i*angle(posy(1,1)));
            posy0(:,coil)=posy0(:,coil).*exp(-1i*angle(posy0(1,1)));
            
            negx(:,coil)=negx(:,coil).*exp(-1i*angle(negx(1,1)));
            negx0(:,coil)=negx0(:,coil).*exp(-1i*angle(negx0(1,1)));
            negy(:,coil)=negy(:,coil).*exp(-1i*angle(negy(1,1)));
            negy0(:,coil)=negy0(:,coil).*exp(-1i*angle(negy0(1,1)));
        end
    end  %coils
    
    
    
    method=1;
    %     %calculate trajectory using the positive shifts only
%          kx4=-1/(2*pi*loc)*(unwrap(angle(posx))-unwrap(angle(posx0)));  %need a 2*pi in the denominator to get to cm^-1
%          ky4=-1/(2*pi*loc)*(unwrap(angle(posy))-unwrap(angle(posy0)));
%     %
%     %     %calculate trajectory using the negative shifts only
%          kx4=1/(2*pi*loc)*(unwrap(angle(negx))-unwrap(angle(negx0)));  %need a 2*pi in the denominator to get to cm^-1
%          ky4=1/(2*pi*loc)*(unwrap(angle(negy))-unwrap(angle(negy0)));
    
    %calculate trajectory using both
    kx4=-1/(4*pi*loc)*unwrap(angle((posx.*negx0)./(posx0.*negx)));
    if 1
    ky4=-1/(4*pi*loc)*unwrap(angle((posy.*negy0)./(posy0.*negy)));% temporary comment out as Emerging has weird 
    else
    ky4=-1/(2*pi*loc)*unwrap(angle((posy./posy0)));
    end
    %kx4=1/(2*pi*loc)*(unwrap(angle(posx))-unwrap(angle(negx)));  %need a 2*pi in the denominator to get to cm^-1
    %ky4=1/(2*pi*loc)*(unwrap(angle(posy))-unwrap(angle(negy)));
    if (0),
        kx4 = kx4 - repmat(kx4(bitshift(size(kx4,1),-1)+1,:),[1080,1]);
        ky4 = ky4 - repmat(ky4(bitshift(size(ky4,1),-1)+1,:),[1080,1]);
    end;
    
    %    kx4=-1/(4*pi*loc)*unwrap(angle((posx./negx)));
    %    ky4=-1/(4*pi*loc)*unwrap(angle((posy./negy)));
    
    %b0 calculations unused...
    b0_kx=-1/(4*pi*loc)*unwrap(angle((posx.*negx)./(posx0.*negx0)));
    b0_ky=-1/(4*pi*loc)*unwrap(angle((posy.*negy)./(posy0.*negy0)));
    b0_kx=-1/(4)*unwrap(angle((posx.*negx)./(posx0.*negx0)));
    b0_ky=-1/(4)*unwrap(angle((posy.*negy)./(posy0.*negy0)));
    
    
    if(1) %combine coils
        %Now we have a kx,ky and kz measurement for each coil
        %Weight each coil by its received signal magnitude (or square of this magnitude?)
        posweights = mean(abs(posx),1).'; posweights=posweights./sum(posweights);
        negweights = mean(abs(negx),1).'; negweights=negweights./sum(negweights);
        
        kx4v2=kx4*posweights;
        posweights2 = posweights.^2;  posweights2=posweights2./sum(posweights2); %square the relative weights
        kx4v3=kx4*posweights2;
        b0_kx_combined = b0_kx*posweights2;
        
        posweights = mean(abs(posy),1).'; posweights=posweights./sum(posweights);
        negweights = mean(abs(negy),1).'; negweights=negweights./sum(negweights);
        
        ky4v2=ky4*posweights;
        posweights2 = posweights.^2;  posweights2=posweights2./sum(posweights2); %square the relative weights
        ky4v3=ky4*posweights2;
        b0_ky_combined = b0_ky*posweights2;
        
    end
    
    %only plot data from the coil with the highest signal
    [~,max_idx_posx]=max(mean(abs(posx),1));
    [~,max_idx_negx]=max(mean(abs(negx),1));
    [~,max_idx_posy]=max(mean(abs(posy),1));
    [~,max_idx_negy]=max(mean(abs(negy),1));
    
    %average the measurements from the positive and negative shifts
    kx4 = kx4(:,max_idx_posx);
    ky4 = ky4(:,max_idx_posy);
    
    
    
    if(doPLOTS)
        kmax=max(sqrt(kx4v3(:).^2+ky4v3(:).^2));
        
        figure(100),subplot(221),plot(abs([posx(:,max_idx_posx) posx0(:,max_idx_posx) posy(:,max_idx_posy) posy0(:,max_idx_posy)]))
        subplot(223),plot(abs([negx(:,max_idx_negx) negx0(:,max_idx_negx) negy(:,max_idx_negy) negy0(:,max_idx_negy)]))
        
        %figure(100),subplot(222),plot(kx4v3,ky4v3),axis (kmax*[-1 1 -1 1 -1 1]),axis square,drawnow
        figure(100),subplot(222),plot(kx4v3,ky4v3),axis (kmax*[-1 1 -1 1]),axis square,drawnow
    end
    drawnow
    % keyboard
    
    kxall(:,proj)=kx4v3;
    kyall(:,proj)=ky4v3;
    kxall_b0(:,proj) = b0_kx_combined;
    kyall_b0(:,proj) = b0_ky_combined ;
end

[kmax,uplimit]=max(sqrt(kxall(:,1).^2+kyall(:,1).^2));

if saveBinary, % save to binary file
    adc_pad_len = numADCpad*2;
    kx_norm = kxall ./ max(max(kxall));
    kx_norm = kyall ./ max(max(kyall));
    kx_norm = kxall ./ max(max(kxall));
    ky_norm = kyall ./ max(max(kyall));
    kx_norm_clean =kx_norm( 1+adc_pad_len : end - adc_pad_len, :);
    ky_norm_clean =ky_norm( 1+adc_pad_len : end - adc_pad_len, :);
    
    spiral_coordinates = zeros(2, 1600, 48);
    spiral_coordinates(1, :,:) = kx_norm_clean(1:1600,:);
    spiral_coordinates(2, :,:) = ky_norm_clean(1:1600,:);
    
    tmpfile1 = 'Traj_148_FOV300_res128_1600_FISP.bin';
    tmpfile2 = 'Traj_148_FOV300_res128.bin';
    filename1 = fullfile(tmppath,tmpfile1);
    filename2 = fullfile(tmppath,tmpfile2);

    
    fid=fopen(filename1,'wb');
    fwrite(fid, single(spiral_coordinates), 'float32');
    fclose(fid);
    
    fid=fopen(filename2,'wb');
    fwrite(fid, single(spiral_coordinates), 'float32');
    fclose(fid);
end;


