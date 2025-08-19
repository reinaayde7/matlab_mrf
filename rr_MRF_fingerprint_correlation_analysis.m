clear all, close all, clc
d.fisp = load('E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\rr_dict_NIST_FISP_t1t2.mat');
d.truefisp = load('E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\rr_dict_NIST_trueFISP_t1t2.mat');

%% correlation between two entries and FISP vs. TRUEFISP
entry = [200 200 0; 200 100 0]; %t1,t2
% entry = [1200 200 0]; %t1,t2

[is_present, col_index] = ismember(entry, d.truefisp.r, 'rows');
d.truefisp.f = squeeze(d.truefisp.dict(:,1,col_index));
[is_present, col_index] = ismember(entry, d.fisp.r, 'rows');
d.fisp.f = squeeze(d.fisp.dict(:,1,col_index));

%%
figure,
subplot(211),plot(abs(d.truefisp.f))
subplot(212), plot(abs(d.fisp.f))


xx = d.truefisp.f;
xxNorm = sqrt(sum(xx.*conj(xx),1));
d.truefisp.f = xx./xxNorm; %normalized fingerprint

xx = d.fisp.f;
xxNorm = sqrt(sum(xx.*conj(xx),1));
d.fisp.f = xx./xxNorm; %normalized fingerprint


%%
crossC_truefisp = xcorr(abs(d.truefisp.f(:,1)),abs(d.truefisp.f(:,2)));
crossC_fisp = xcorr(abs(d.fisp.f(:,1)),abs(d.fisp.f(:,2)));
figure, plot((crossC_truefisp)), hold on, plot((crossC_fisp)), legend('trueFISP', 'FISP')

%%

d.fisp = load('E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\rr_dict_NIST_FISP_t1t2.mat');
d.truefisp = load('E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\rr_dict_NIST_trueFISP_t1t2B0v3.mat'); %takes long to import - big ass dataset
%% correlation for TRUEFISP, same T1/T2 entry but different B0 offsets
entry = [1200 200 0; 1200 200 5; 1200 200 10; 1200 200 15; 1200 200 20; 1200 200 25]; %t1,t2
[is_present, col_index] = ismember(entry, d.truefisp.r, 'rows');
d.truefisp.f = squeeze(d.truefisp.dict(:,1,col_index));

styles = {'k'; 'b'; 'r'; 'g'; 'y';'m'};
figure,
for i =1:6
    plot(abs(d.truefisp.f(:,i)), 'Color', styles{i}), hold on
end    
legend('0','5','10','15','20','25')
%%
entry = [1200 200 0; 1200 200 1; 1200 200 2; 1200 200 3; 1200 200 4; 1200 200 5]; %t1,t2
[is_present, col_index] = ismember(entry, d.truefisp.r, 'rows');
d.truefisp.f = squeeze(d.truefisp.dict(:,1,col_index));

styles = {'k'; 'b'; 'r'; 'g'; 'y';'m'};
figure,
for i =1:6
    plot(abs(d.truefisp.f(:,i)), 'Color', styles{i}), hold on
end    
legend('0','1','2','3','4','5')

%%
entry = [1200 200 10; 1200 200 11; 1200 200 12; 1200 200 13; 1200 200 14; 1200 200 15]; %t1,t2
[is_present, col_index] = ismember(entry, d.truefisp.r, 'rows');
d.truefisp.f = squeeze(d.truefisp.dict(:,1,col_index));

styles = {'k'; 'b'; 'r'; 'g'; 'y';'m'};
figure,
for i =1:6
    plot(abs(d.truefisp.f(:,i)), 'Color', styles{i}), hold on
end    
legend('10','11','12','13','14','15')