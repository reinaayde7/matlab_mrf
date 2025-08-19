% Functions: load Siemens Raw data for Siemens VD System
% The kspace data arrange as
% [rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.nshots,rawinfo.nframes]
%
% Version 0.1 -- The first version --12.06.2015
% Version 0.2 -- Add (Number of Average) and (Number of Repetition) in
% outputs -- 12.08.2015
% Version 0.3 -- update mapVBVD (Sept 2015) version
%             -- read the noise from the multiraid data
% Version 0.4 -- read the noise from single raid.
% Version 0.5 -- export PMUtime,and timestamp from MDH  --12.29.15
%             -- 2.5 ms per size for PMUtime
%             -- triggerdelay, and triggerwindow
% Version 0.6 -- add TR and TI into rawinfo struct -- 1/4/2016
% Version 0.7 -- add TE into rawinfo struct --03/04/2018

% Yun Jiang -- yun.jiang@case.edu -- All copyright reserved
%

%[COLUMNS,COILS,SHOTS,1,1,AVERAGES,1,1,REPETITIONS]
function [raw,noise,rawinfo,noiseinfo] = loadSiemensRawData(filename)
%if nargin < 2, doCoilCompression = [];end;
if nargin < 1, filename = []; end;
%if isempty(doCoilCompression), doCoilCompression = 0; end;
if isempty(filename),
    [tempfilename,temppathname] = uigetfile('*.dat','Select the Siemens raw data file');
    filename = fullfile(temppathname,tempfilename);
end;

datastruct = mapVBVD(filename);
switch length(datastruct),
    case 1
        raw = squeeze(datastruct.image());
        % Acquision parameters
        try
            rawinfo.ncoils         = length(datastruct.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
        catch
            rawinfo.ncoils         = length(datastruct.hdr.MeasYaps.asCoilSelectMeas{1,1}.asList);
        end
        rawinfo.nshots         = datastruct.hdr.Config.NLinMeas;
        rawinfo.nframes        = datastruct.hdr.Config.NSetMeas;
        rawinfo.matrix         = datastruct.hdr.Config.BaseResolution;
        rawinfo.fieldofview    = datastruct.hdr.Config.RoFOV;
        rawinfo.repetitiontime = datastruct.hdr.Config.TR;
        rawinfo.inversiontime  = datastruct.hdr.Config.TI;
        rawinfo.nslices        = datastruct.hdr.Config.NSlcMeas;
        rawinfo.naverages      = datastruct.hdr.Config.NAveMeas;
        rawinfo.nrepetitions   = datastruct.hdr.Config.NRepMeas;
        rawinfo.nADCpoints     = size(raw,1);
        rawinfo.TI             = datastruct.hdr.Config.TI;% Yun Jiang -- 1/4/2016
        rawinfo.TR             = datastruct.hdr.Config.TR;% Yun Jiang -- 1/4/2016
        rawinfo.Nex            = datastruct.hdr.Config.NSetMeas; % Yun Jiang -- 1/4/2016
        rawinfo.TE             = datastruct.hdr.MeasYaps.alTE{1,1}; % Yun Jiang - -3/04/2018
        try
        rawinfo.InPlaneRot     = datastruct.hdr.Config.InPlaneRot; % Yun Jiang 03/30/2018
        catch
        end
        rawinfo.SliceThickness = datastruct.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness; % Yun JIang 03/30/2018
        rawinfo.npartitions    = datastruct.hdr.Config.NPar; % Yun Jiang 05252018

        try 
            rawinfo.idproj         = datastruct{1,1}.image.iceParam(20,:); % Yun Jiang - 2/13/2018
        catch
            rawinfo.idproj     =[];
        end;
        % PMU time
        rawinfo.pmutime = datastruct.image.pmutime;
        rawinfo.timestamp = datastruct.image.timestamp;
        try
            rawinfo.triggerdelay = datastruct.hdr.Phoenix.sPhysioImaging.sPhysioECG.lTriggerDelay; % in us
            rawinfo.scanWindow = datastruct.hdr.Phoenix.sPhysioImaging.sPhysioECG.lScanWindow;% in ms
        catch
            rawinfo.triggerdelay = [];
            rawinfo.scanWindow = [];
        end
        
        try
            noise = squeeze(datastruct.noise());
            noiseinfo = rawinfo;
            noise = reshape(noise,[noiseinfo.nADCpoints,noiseinfo.ncoils,noiseinfo.nshots,...
                numel(noise)/(noiseinfo.nADCpoints*noiseinfo.ncoils*noiseinfo.nshots)]);
            
        catch
            noise = [];
            noiseinfo = [];
        end;
        clear datastruct;
        
        if rawinfo.nslices > 1
        raw = reshape(raw,[rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.nshots,rawinfo.nslices,rawinfo.nframes]);
        else
        raw = reshape(raw,[rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.nshots,rawinfo.nframes]);
        end
            
    case 2
        % noise
        noise = squeeze(datastruct{1,1}.noise());
        noiseinfo.ncoils  = length(datastruct{1,1}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
        noiseinfo.nshots  = datastruct{1,1}.hdr.Config.NLinMeas;
        noiseinfo.nframes = datastruct{1,1}.hdr.Config.NSetMeas;
        noiseinfo.matrix  = datastruct{1,1}.hdr.Config.BaseResolution;
        noiseinfo.fieldofview = datastruct{1,1}.hdr.Config.RoFOV;
        noiseinfo.repetitiontime = datastruct{1,1}.hdr.Config.TR;
        noiseinfo.inversiontime  = datastruct{1,1}.hdr.Config.TI;
        noiseinfo.nslices = datastruct{1,1}.hdr.Config.NSlcMeas;
        noiseinfo.naverages = datastruct{1,1}.hdr.Config.NAveMeas;
        noiseinfo.nrepetitions = datastruct{1,1}.hdr.Config.NRepMeas;
        noiseinfo.nADCpoints = size(noise,1);
        
        
        % raw data
        raw = squeeze(datastruct{1,2}.image());
        % Acquision parameters
        rawinfo.ncoils         = length(datastruct{1,2}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
        rawinfo.nshots         = datastruct{1,2}.hdr.Config.NLinMeas;
        rawinfo.nframes        = datastruct{1,2}.hdr.Config.NSetMeas;
        rawinfo.matrix         = datastruct{1,2}.hdr.Config.BaseResolution;
        rawinfo.fieldofview    = datastruct{1,2}.hdr.Config.RoFOV;
        rawinfo.repetitiontime = datastruct{1,2}.hdr.Config.TR;
        rawinfo.inversiontime  = datastruct{1,2}.hdr.Config.TI;
        rawinfo.nslices        = datastruct{1,2}.hdr.Config.NSlcMeas;
        rawinfo.naverages      = datastruct{1,2}.hdr.Config.NAveMeas;
        rawinfo.nrepetitions   = datastruct{1,2}.hdr.Config.NRepMeas;
        rawinfo.nADCpoints     = size(raw,1);
        rawinfo.TI             = datastruct{1,2}.hdr.Config.TI;% Yun Jiang -- 1/4/2016
        rawinfo.TR             = datastruct{1,2}.hdr.Config.TR;% Yun Jiang -- 1/4/2016
        rawinfo.Nex            = datastruct{1,2}.hdr.Config.NSetMeas; % Yun Jiang -- 1/4/2016
        rawinfo.TE             = datastruct{1,2}.hdr.MeasYaps.alTE{1,1}; % Yun Jiang - -3/04/2018
        rawinfo.SliceThickness = datastruct{1,2}.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness; % Yun JIang 03/30/2018
        rawinfo.npartitions    = datastruct{1,2}.hdr.Config.NPar; % Yun Jiang 05252018


        try 
            rawinfo.idproj         = datastruct{1,2}.image.iceParam(20,:); % Yun Jiang - 2/13/2018
            rawinfo.InPlaneRot     = datastruct{1,2}.hdr.Config.InPlaneRot; % Yun Jiang 03/30/2018

        catch
            rawinfo.idproj     =[];
            rawinfo.InPlaneRot =[];
        end;

        
        
        rawinfo.pmutime = datastruct{1,2}.image.pmutime;
        rawinfo.timestamp = datastruct{1,2}.image.timestamp;
        try
            rawinfo.triggerdelay = datastruct{1,2}.hdr.Phoenix.sPhysioImaging.sPhysioECG.lTriggerDelay; % in us
            rawinfo.scanWindow = datastruct{1,2}.hdr.Phoenix.sPhysioImaging.sPhysioECG.lScanWindow;% in ms
        catch
            rawinfo.triggerdelay = [];
            rawinfo.scanWindow = [];
        end
        clear datastruct;
        
        try
            raw = reshape(raw,[rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.nshots,rawinfo.nframes]);
        catch
            raw = raw;
            %raw = reshape(raw,[rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.nshots,...
            %    numel(raw)/(rawinfo.nADCpoints*rawinfo.ncoils*rawinfo.nshots)]);
        end
        noise = reshape(noise,[noiseinfo.nADCpoints,noiseinfo.ncoils,noiseinfo.matrix*2,...
            numel(noise)/(noiseinfo.nADCpoints*noiseinfo.ncoils*noiseinfo.matrix*2)]);
        
    case 3
        %for WIP: 1st datastruct is B1 prescan, 2nd datastruct is noise prescan or coilsenitiviy map; 3rd is the raw data
        % 1st B1 prescan from body array b1scan is not an out put right now
        % in this code.
        b1scan = squeeze(datastruct{1,1}.image());
        b1scaninfo.ncoils  = length(datastruct{1,1}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
        b1scaninfo.nshots  = datastruct{1,1}.hdr.Config.NLinMeas;
        b1scaninfo.nframes = datastruct{1,1}.hdr.Config.NSetMeas;
        b1scaninfo.matrix  = datastruct{1,1}.hdr.Config.BaseResolution;
        b1scaninfo.fieldofview = datastruct{1,1}.hdr.Config.RoFOV;
        b1scaninfo.repetitiontime = datastruct{1,1}.hdr.Config.TR;
        b1scaninfo.inversiontime  = datastruct{1,1}.hdr.Config.TI;
        b1scaninfo.nslices = datastruct{1,1}.hdr.Config.NSlcMeas;
        b1scaninfo.naverages = datastruct{1,1}.hdr.Config.NAveMeas;
        b1scaninfo.nrepetitions = datastruct{1,1}.hdr.Config.NRepMeas;
        b1scaninfo.nADCpoints = size(b1scan,1);
        
        
        % noise
        try
            noise = squeeze(datastruct{1,2}.noise());
            noiseinfo.ncoils  = length(datastruct{1,2}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
            noiseinfo.nshots  = datastruct{1,2}.hdr.Config.NLinMeas;
            noiseinfo.nframes = datastruct{1,2}.hdr.Config.NSetMeas;
            noiseinfo.matrix  = datastruct{1,2}.hdr.Config.BaseResolution;
            noiseinfo.fieldofview = datastruct{1,2}.hdr.Config.RoFOV;
            noiseinfo.repetitiontime = datastruct{1,2}.hdr.Config.TR;
            noiseinfo.inversiontime  = datastruct{1,2}.hdr.Config.TI;
            noiseinfo.nslices = datastruct{1,2}.hdr.Config.NSlcMeas;
            noiseinfo.naverages = datastruct{1,2}.hdr.Config.NAveMeas;
            noiseinfo.nrepetitions = datastruct{1,2}.hdr.Config.NRepMeas;
            noiseinfo.nADCpoints = size(noise,1);
        catch
            noise = squeeze(datastruct{1,1}.noise());
            noiseinfo.ncoils  = length(datastruct{1,1}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
            noiseinfo.nshots  = datastruct{1,1}.hdr.Config.NLinMeas;
            noiseinfo.nframes = datastruct{1,1}.hdr.Config.NSetMeas;
            noiseinfo.matrix  = datastruct{1,1}.hdr.Config.BaseResolution;
            noiseinfo.fieldofview = datastruct{1,1}.hdr.Config.RoFOV;
            noiseinfo.repetitiontime = datastruct{1,1}.hdr.Config.TR;
            noiseinfo.inversiontime  = datastruct{1,1}.hdr.Config.TI;
            noiseinfo.nslices = datastruct{1,1}.hdr.Config.NSlcMeas;
            noiseinfo.naverages = datastruct{1,1}.hdr.Config.NAveMeas;
            noiseinfo.nrepetitions = datastruct{1,1}.hdr.Config.NRepMeas;
            noiseinfo.nADCpoints = size(noise,1);
        end
        
        
        % raw data
        raw = squeeze(datastruct{1,3}.image());
        % Acquision parameters
        rawinfo.ncoils         = length(datastruct{1,3}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
        rawinfo.nshots         = datastruct{1,3}.hdr.Config.NLinMeas;
        rawinfo.nframes        = datastruct{1,3}.hdr.Config.NSetMeas;
        rawinfo.matrix         = datastruct{1,3}.hdr.Config.BaseResolution;
        rawinfo.fieldofview    = datastruct{1,3}.hdr.Config.RoFOV;
        rawinfo.repetitiontime = datastruct{1,3}.hdr.Config.TR;
        rawinfo.inversiontime  = datastruct{1,3}.hdr.Config.TI;
        rawinfo.nslices        = datastruct{1,3}.hdr.Config.NSlcMeas;
        rawinfo.naverages      = datastruct{1,3}.hdr.Config.NAveMeas;
        rawinfo.nrepetitions   = datastruct{1,3}.hdr.Config.NRepMeas;
        rawinfo.nADCpoints     = size(raw,1);
        rawinfo.TI             = datastruct{1,3}.hdr.Config.TI;% Yun Jiang -- 1/4/2016
        rawinfo.TR             = datastruct{1,3}.hdr.Config.TR;% Yun Jiang -- 1/4/2016
        rawinfo.Nex            = datastruct{1,3}.hdr.Config.NSetMeas; % Yun Jiang -- 1/4/2016
        rawinfo.TE             = datastruct{1,3}.hdr.MeasYaps.alTE{1,1}; % Yun Jiang - -3/04/2018
        rawinfo.SliceThickness = datastruct{1,3}.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness; % Yun JIang 03/30/2018
        rawinfo.npartitions    = datastruct{1,3}.hdr.Config.NPar; % Yun Jiang 05252018



        try 
            rawinfo.idproj         = datastruct{1,3}.image.iceParam(20,:); % Yun Jiang - 2/13/2018
            rawinfo.InPlaneRot     = datastruct{1,3}.hdr.Config.InPlaneRot; % Yun Jiang 03/30/2018        

        catch
            rawinfo.idproj     =[];
            rawinfo.InPlaneRot =[];
        end;
        
        rawinfo.pmutime = datastruct{1,3}.image.pmutime;
        rawinfo.timestamp = datastruct{1,3}.image.timestamp;
        try
            rawinfo.triggerdelay = datastruct{1,3}.hdr.Phoenix.sPhysioImaging.sPhysioECG.lTriggerDelay; % in us
            rawinfo.scanWindow = datastruct{1,3}.hdr.Phoenix.sPhysioImaging.sPhysioECG.lScanWindow;% in ms
        catch
            rawinfo.triggerdelay = [];
            rawinfo.scanWindow = [];
        end
        clear datastruct;
        
        
        try
            raw = reshape(raw,[rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.Nex,rawinfo.nslices,rawinfo.nrepetitions]);
        catch
            raw = raw;
            %raw = reshape(raw,[rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.Nex,rawinfo.nslices,...
            %    numel(raw)/(rawinfo.nADCpoints*rawinfo.ncoils*rawinfo.Nex*rawinfo.nslices)]);
        end
        noise = reshape(noise,[noiseinfo.nADCpoints,noiseinfo.ncoils,noiseinfo.matrix*2,...
            numel(noise)/(noiseinfo.nADCpoints*noiseinfo.ncoils*noiseinfo.matrix*2)]);
        
        
    case 4
        %for WIP: 1st datastruct is B1 prescan, 2nd datastruct is noise prescan or coilsenitiviy map; 3rd is the raw data
        % 1st B1 prescan from body array b1scan is not an out put right now
        % in this code.
        try
            b1scan = squeeze(datastruct{1,1}.image());
            b1scaninfo.ncoils  = length(datastruct{1,1}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
            b1scaninfo.nshots  = datastruct{1,1}.hdr.Config.NLinMeas;
            b1scaninfo.nframes = datastruct{1,1}.hdr.Config.NSetMeas;
            b1scaninfo.matrix  = datastruct{1,1}.hdr.Config.BaseResolution;
            b1scaninfo.fieldofview = datastruct{1,1}.hdr.Config.RoFOV;
            b1scaninfo.repetitiontime = datastruct{1,1}.hdr.Config.TR;
            b1scaninfo.inversiontime  = datastruct{1,1}.hdr.Config.TI;
            b1scaninfo.nslices = datastruct{1,1}.hdr.Config.NSlcMeas;
            b1scaninfo.naverages = datastruct{1,1}.hdr.Config.NAveMeas;
            b1scaninfo.nrepetitions = datastruct{1,1}.hdr.Config.NRepMeas;
            b1scaninfo.nADCpoints = size(b1scan,1);
        catch
            b1scan = squeeze(datastruct{1,2}.image());
            b1scaninfo.ncoils  = length(datastruct{1,2}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
            b1scaninfo.nshots  = datastruct{1,2}.hdr.Config.NLinMeas;
            b1scaninfo.nframes = datastruct{1,2}.hdr.Config.NSetMeas;
            b1scaninfo.matrix  = datastruct{1,2}.hdr.Config.BaseResolution;
            b1scaninfo.fieldofview = datastruct{1,2}.hdr.Config.RoFOV;
            b1scaninfo.repetitiontime = datastruct{1,2}.hdr.Config.TR;
            b1scaninfo.inversiontime  = datastruct{1,2}.hdr.Config.TI;
            b1scaninfo.nslices = datastruct{1,2}.hdr.Config.NSlcMeas;
            b1scaninfo.naverages = datastruct{1,2}.hdr.Config.NAveMeas;
            b1scaninfo.nrepetitions = datastruct{1,2}.hdr.Config.NRepMeas;
            b1scaninfo.nADCpoints = size(b1scan,1);
        end;
        
        
        % noise
        try
            noise = squeeze(datastruct{1,2}.noise());
            noiseinfo.ncoils  = length(datastruct{1,2}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
            noiseinfo.nshots  = datastruct{1,2}.hdr.Config.NLinMeas;
            noiseinfo.nframes = datastruct{1,2}.hdr.Config.NSetMeas;
            noiseinfo.matrix  = datastruct{1,2}.hdr.Config.BaseResolution;
            noiseinfo.fieldofview = datastruct{1,2}.hdr.Config.RoFOV;
            noiseinfo.repetitiontime = datastruct{1,2}.hdr.Config.TR;
            noiseinfo.inversiontime  = datastruct{1,2}.hdr.Config.TI;
            noiseinfo.nslices = datastruct{1,2}.hdr.Config.NSlcMeas;
            noiseinfo.naverages = datastruct{1,2}.hdr.Config.NAveMeas;
            noiseinfo.nrepetitions = datastruct{1,2}.hdr.Config.NRepMeas;
            noiseinfo.nADCpoints = size(noise,1);
        catch
            noise = squeeze(datastruct{1,1}.noise());
            noiseinfo.ncoils  = length(datastruct{1,1}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
            noiseinfo.nshots  = datastruct{1,1}.hdr.Config.NLinMeas;
            noiseinfo.nframes = datastruct{1,1}.hdr.Config.NSetMeas;
            noiseinfo.matrix  = datastruct{1,1}.hdr.Config.BaseResolution;
            noiseinfo.fieldofview = datastruct{1,1}.hdr.Config.RoFOV;
            noiseinfo.repetitiontime = datastruct{1,1}.hdr.Config.TR;
            noiseinfo.inversiontime  = datastruct{1,1}.hdr.Config.TI;
            noiseinfo.nslices = datastruct{1,1}.hdr.Config.NSlcMeas;
            noiseinfo.naverages = datastruct{1,1}.hdr.Config.NAveMeas;
            noiseinfo.nrepetitions = datastruct{1,1}.hdr.Config.NRepMeas;
            noiseinfo.nADCpoints = size(noise,1);
        end
        
        
        % raw data
        raw = squeeze(datastruct{1,4}.image());
        % Acquision parameters
        rawinfo.ncoils         = length(datastruct{1,4}.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList);
        rawinfo.nshots         = datastruct{1,4}.hdr.Config.NLinMeas;
        rawinfo.nframes        = datastruct{1,4}.hdr.Config.NSetMeas;
        rawinfo.matrix         = datastruct{1,4}.hdr.Config.BaseResolution;
        rawinfo.fieldofview    = datastruct{1,4}.hdr.Config.RoFOV;
        rawinfo.repetitiontime = datastruct{1,4}.hdr.Config.TR;
        rawinfo.inversiontime  = datastruct{1,4}.hdr.Config.TI;
        rawinfo.nslices        = datastruct{1,4}.hdr.Config.NSlcMeas;
        rawinfo.naverages      = datastruct{1,4}.hdr.Config.NAveMeas;
        rawinfo.nrepetitions   = datastruct{1,4}.hdr.Config.NRepMeas;
        rawinfo.nADCpoints     = size(raw,1);
        rawinfo.TI             = datastruct{1,4}.hdr.Config.TI;% Yun Jiang -- 1/4/2016
        rawinfo.TR             = datastruct{1,4}.hdr.Config.TR;% Yun Jiang -- 1/4/2016
        rawinfo.Nex            = datastruct{1,4}.hdr.Config.NSetMeas; % Yun Jiang -- 1/4/2016        
        rawinfo.TE             = datastruct{1,4}.hdr.MeasYaps.alTE{1,1}; % Yun Jiang - -3/04/2018
        rawinfo.SliceThickness = datastruct{1,4}.hdr.MeasYaps.sSliceArray.asSlice{1,1}.dThickness; % Yun JIang 03/30/2018
        rawinfo.npartitions    = datastruct{1,4}.hdr.Config.NPar; % Yun Jiang 05252018


        
        
        rawinfo.pmutime = datastruct{1,4}.image.pmutime;
        rawinfo.timestamp = datastruct{1,4}.image.timestamp;
        try
            rawinfo.triggerdelay = datastruct{1,4}.hdr.Phoenix.sPhysioImaging.sPhysioECG.lTriggerDelay; % in us
            rawinfo.scanWindow = datastruct{1,4}.hdr.Phoenix.sPhysioImaging.sPhysioECG.lScanWindow;% in ms
            rawinfo.InPlaneRot     = datastruct{1,4}.hdr.Config.InPlaneRot; % Yun Jiang 03/30/2018

        catch
            rawinfo.triggerdelay = [];
            rawinfo.scanWindow   = [];
            rawinfo.InPlaneRot   = [];
        end
        clear datastruct;
        
        
        try
            raw = reshape(raw,[rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.Nex,rawinfo.nslices,rawinfo.nrepetitions]);
        catch
            raw = raw;
            
            %raw = reshape(raw,[rawinfo.nADCpoints,rawinfo.ncoils,rawinfo.Nex,rawinfo.nslices,...
            %    numel(raw)/(rawinfo.nADCpoints*rawinfo.ncoils*rawinfo.Nex*rawinfo.nslices)]);
        end
        noise = reshape(noise,[noiseinfo.nADCpoints,noiseinfo.ncoils,noiseinfo.matrix*2,...
            numel(noise)/(noiseinfo.nADCpoints*noiseinfo.ncoils*noiseinfo.matrix*2)]);
        
    otherwise
        error('Seriously!...Wrong dimension of input data');
end





