function [Mx,My,Mz]=dictionary_TRUEFISP_withRelaxDecay(r,flip,tr,te,phase,ti,...
    spins,frames,Inv,waitingDuration)

%Rudy: phase twisting is commented out since bSSFP does not add spoilers.
%We commented them out to save sim time.

cnt = size(r,1);
nt = length(flip);

if (nt>frames),
    nt = frames;
end;

%df=0;

Mx = zeros(nt*length(te),cnt,'single');
My = zeros(nt*length(te),cnt,'single');
Mz = zeros(nt*length(te),cnt,'single');

% no need for phase twisting in trueFISP
% phirange = linspace(-phasetwist/2,phasetwist/2,spins);
% Z = zeros(3,3,spins);
% for s = 1:spins
%     Z(:,:,s) = zrot(phirange(s));
% end

onemat = ones(1,spins);


for n = 1:cnt %Rudy debug
    
    
    T1 = r(n,1)/1000;
    T2 = r(n,2)/1000;
    df = r(n,3);
    
    M0 = repmat([0 0 1]',1,spins,length(te)+1);
    M = repmat([0 0 0]',1,spins);
    
    for iEcho = 1:length(te) %Rudy debug
        %iEcho
        for k = 1:nt %Rudy debug
            % disp(k)  %Rudy debug
            if (k == 1)
                if Inv == 1
                    M = throt(pi,0)*M0(:,:,iEcho);
                end
                [Ati,Bti] = freeprecess(ti,T1,T2,df);
                Bti = Bti*onemat;
                M = Ati*M+Bti;
                
                %Rudy 240627: spoiling not used in TRUEFISP + commented out to save sim time
                % spoiling
                % for s=1:spins
                %     M(:,s)=Z(:,:,s)*M(:,s);
                % end
            end
            
            
            
            
            % alpha pulse
            M = throt(flip(k),phase(k))*M;
            [Ate,Bte] = freeprecess(te(iEcho),T1,T2,df);
            Bte = Bte*onemat;
            M = Ate*M+Bte;

            %Rudy: rephase signal (ADC phase == RF phase)
            % to speed up, this is taken care outside spin simulation
            % M = zrot(-phase(k))*M;

            Mx(nt*(iEcho-1)+k,n) = mean(M(1,:),2);
            My(nt*(iEcho-1)+k,n) = mean(M(2,:),2);
            Mz(nt*(iEcho-1)+k,n) = mean(M(3,:),2);
            
            %Rudy: dephase signal (RF oscillation)
            % to speed up, this is taken care outside spin simulation
            % M = zrot(phase(k))*M;

            [A,B]=freeprecess(tr(k)-te(iEcho),T1,T2,df);
            B = B*onemat;
            M=A*M+B;

            % Rudy 240627: spoiler for FISP taken away
            % for s=1:spins
            %     M(:,s)=Z(:,:,s)*M(:,s);
            % end
            
            
            
        end % timepoints
        
        % relaxation time during the waiting time
        [A,B]=freeprecess(waitingDuration,T1,T2,df);
        B = B*onemat;
        M0(:,:,iEcho+1)=A*M+B;
        % disp('done1')
    end
    %end
    % disp('done2')
end % T1 T2 entry
