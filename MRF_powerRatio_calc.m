
% dictionary_name = 'dict_trueFISP_TE1_8_TR13_FAbody15_Nex600_b0100';
dictionary_name = 'dict_FISP_TE1_8_TR14_FAbody10_Nex600';
dfile = load(['E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\' dictionary_name '.mat'],'dict','r','log');
dict_10 = dfile.dict; r_10 = dfile.r; dict_log_10 = dfile.log; 
clearvars dfile        

%%
% dictionary_name = 'dict_FISP_TE1_8_TR14_FAbody15_Nex600';
dictionary_name = 'dict_FISP_TE1_8_TR14_FAbody15_Nex600';
dfile = load(['E:\POSTDOC_UoM\08_Project_MRF\_rr_simulated_dictionaries\' dictionary_name '.mat'],'dict','r','log');
dict_15 = dfile.dict; r_15 = dfile.r; dict_log_15 = dfile.log;
clearvars dfile
%%

t1 = 490;
t2=50;

r_idx_15 = find(r_15(:,1) == t1 & r_15(:,2) == t2);
r_idx_10 = find(r_10(:,1) == t1 & r_10(:,2) == t2);

s_fisp_10 = dict_10(:,1,r_idx_10);
powerfisp_10 = mean(abs(s_fisp_10).^2);

s_fisp_15 = dict_15(:,1,r_idx_15);
powerfisp_15 = mean(abs(s_fisp_15).^2);

%%
x = 1:600;
figure,
set(gcf, 'Color', 'white'), hold on, 
plot(x,abs(s_fisp_10), 'b', 'LineWidth',2), hold on
plot(x,abs(s_fisp_15), 'r', 'LineWidth',2)
a = sprintf('FISP 10 avg.pow: %.2f', mean(powerfisp_10)*1e2);
b = sprintf('FISP 15 avg.pow: %.2f', mean(powerfisp_15)*1e2);
legend(a,b)
title(sprintf('T1: %.0f ms, T2: %.0f ms', t1, t2))
ylim([0, 0.45])
xlabel('Image Time Points')
ylabel('|S|')


%%
for idx=1:size(dict_10,3)
    p10(idx) = mean(abs( dict_10(:,1,idx)).^2);
    p15(idx) = mean(abs( dict_15(:,1,idx)).^2);
end

%%
rp = p15./p10;
%%
x = r_10(:,1); y = r_10(:,2);
z = rp;

idx = find(z<1.01 & z>0.99);
%%
figure, scatter(x,y,100, z, 'filled'), hold on
scatter(x(idx),y(idx),20,z(idx),'k')
colormap('jet'), colorbar
xlabel('T1 [ms]'), ylabel('T2 [ms]')
title('Average Power Ratio FA1.5/1.0')


%% FISP v bSSFP

