% correlation matrix for dictionary


%% FISP
%NIST t2 array
t1= [1990, 1755, 1520, 1240, 1020, 790, 600, 460, 345, 230, 180, 125, 95, 65];
t2= [1065,  798,  580,  415,  308, 215, 150, 110,  80,  50,  40,  25, 20, 15];

for pair = 1:length(t1)
    t_diff = abs(r(:,1)-t1(pair)) + abs(r(:,2)-t2(pair));
    ii = find(t_diff == min(t_diff));
    idx(pair) = ii(1); %retain one solution only
end

% % xxNorm = sqrt(sum(xx.*conj(xx),2));
% % normAll = xxNorm*dictnorm; %clear xxNorm
% % innerProduct = conj(xx)*dd./normAll; 

% cmtxNorm = sqrt(sum(dict_subspace.*conj(dict_subspace),2));
dict_subspace = squeeze(dict(:,1,idx));

for line=1:length(t1)
    signal = imag(dict_subspace(:,line)) - j*real(dict_subspace(:,line));     % conversion needed to use old code based on invivo data;
    [tt1,tt2,bb0,iP] = patternmatchDEBUG(signal',1,r,0,dict,dictnorm,1);
    cmtx(:,line) = iP(idx);
end

tstring = strcat(string(t1), '/', string(t2));
figure, heatmap(tstring,tstring,abs(cmtx))








%% trueFISP (we need to include B0 considerations)
t1= [1990, 1755, 1520, 1240, 1020, 790, 600, 460, 345, 230, 180, 125, 95, 65];
t2= [1065,  798,  580,  415,  308, 215, 150, 110,  80,  50,  40,  25, 20, 15];
b0= [0,       0,    0,    0,    0,   0,   0,   0,   0,   0,   0,   0,  0,  0,];


[idxs,values, fp] = dict_entry_lookup(dict, r, t1, t2, b0);

%%
doPlot = 0;
% dictO = dict;

%%
idxV = idxs(:); %vector of all indexes we want to explore for correlation (grouped by B0 shim value)

% dict = dictO;
dict_subspace = squeeze(dict(:,1,idxV)); 

f = waitbar(0, 'Processing...'); tic;
for v=1:length(idxV) 
    waitbar(v/length(idxV), f, sprintf('entry: %d/%d', v, length(idxV) ));
    signal = imag(dict_subspace(:,v)) - j*real(dict_subspace(:,v));     % conversion needed to use old code based on invivo data;
    [tt1,tt2,bb0,iP] = patternmatchDEBUG(signal',1,r,0,dict,dictnorm,1,doPlot);
    cmtx(:,v) = iP(idxV);
end
close(f); 

%%
tstring = strcat(string(t1), '/', string(t2));
figure, heatmap(tstring,tstring,abs(squeeze(cmtx())))
%%
figure, heatmap(abs(cmtx))
%%
figure, plot(abs(cmtx(1,:)))


%%

figure, plot(imag(iP))
%%
doPlot = 1;
%select 1 fingerprint (one that fails)
set1 = 7; set2 = 4;
idx_sel = idx(set1,set2);
values(set1,set2)

x = dict(:,1,idx_sel);
x = imag(x) - j*real(x);     % conversion needed to use old code based on invivo data;
figure, plot(abs(x)), title(['fingerprint: ' values(set1,set2)])

[tt1,tt2,bb0,iP] = patternmatchDEBUG(x',1,r,0,dict,dictnorm,1,doPlot);
figure, plot(abs(iP)), title(['correlation: ' values(set1,set2)])

%%
%63 600/150/0
%146 354/80/0

possible_sol = find(abs(iP)>0.996);
figure, plot(abs(iP(possible_sol))), ylabel('Inner Prod'), xlabel('possible solutions')

figure, 
subplot(311), plot(r(possible_sol,1)), ylabel('T1 [ms]')
subplot(312), plot(r(possible_sol,2)), ylabel('T2 [ms]')
subplot(313), plot(r(possible_sol,3)), ylabel('B0 [Hz]'), xlabel('possible solutions')