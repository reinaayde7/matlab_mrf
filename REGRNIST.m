function [means_t1, means_t2]=REGRNIST(NIST, eval, th_t1, th_t2)

set(0,'defaultfigurecolor','w');
set(0,'defaultaxesfontname','arial');
set(0,'defaultaxesfontsize',13);

th_t1_1 = th_t1(1); th_t1_2 = th_t1(2);
th_t2_1 = th_t2(1); th_t2_2 = th_t2(2);
% n = size(rois, 2);

% means_t1 = zeros(n, 1);
% means_t2 = zeros(n, 1);
% for ii = 1:n
%     means_t1(ii) = mean(t1(rois{ii}));
%     means_t2(ii) = mean(t2(rois{ii}));
% end

means_t1 = eval.means_t1;
means_t2 = eval.means_t2;
% ------------------------ T1s ----------------------------------

%regression
coeffs = polyfit(log(NIST(th_t1_1:th_t1_2, 1)), log(means_t1(th_t1_1:th_t1_2)), 1);
yfit = polyval(coeffs,log(NIST(th_t1_1:th_t1_2, 1)));
SStot= sum( ( log(means_t1(th_t1_1:th_t1_2))-mean(log(means_t1(th_t1_1:th_t1_2))) ).^2 );
SSres=sum( ( log(means_t1(th_t1_1:th_t1_2))-yfit ).^2 );
Rsq = 1-SSres/SStot;

f=figure;
f.Position = [0 0 800 1000];
t = tiledlayout(3, 2, "TileSpacing", "tight");

ax = nexttile;

plot(log(NIST(th_t1_1:th_t1_2, 1)), log(NIST(th_t1_1:th_t1_2, 1)), LineWidth=2), hold on
plot(log(NIST(th_t1_1:th_t1_2, 1)), yfit, LineWidth=1.5, Color='k'), hold on
plot(log(NIST(th_t1_1:th_t1_2, 1)), log(means_t1(th_t1_1:th_t1_2)), '.', MarkerSize=30, Color='r');

text_str = sprintf('y = %.3fx + %.3f \nR^2 = %.3f', coeffs(1), coeffs(2), Rsq);
annotation('textbox', [0.15, 0.835, 0.21, 0.07], 'String', text_str, 'FontSize',13);
% annotation('textbox', [0.02, 0.95, 0.2, 0.1], 'String', text_str, ...
%            'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k');
xaxis = flipud(round(NIST(th_t1_1:th_t1_2, 1)));
xticks(flipud(log(NIST(th_t1_1:th_t1_2, 1))))
xticklabels(xaxis)
yticks(flipud(log(NIST(th_t1_1:th_t1_2, 1))))
yticklabels(xaxis)
xlim(fliplr([log(NIST(th_t1_1, 1)), log(NIST(th_t1_2, 1))]))
ylim(fliplr([log(NIST(th_t1_1, 1)), log(NIST(th_t1_2, 1))]))
title('T1 [ms]')
set(ax,'xticklabel',[])


% ------------------------ T2s ----------------------------------
coeffs = polyfit(log(NIST(th_t2_1:th_t2_2, 2)), log(means_t2(th_t2_1:th_t2_2)), 1);
yfit = polyval(coeffs,log(NIST(th_t2_1:th_t2_2, 2)));
SStot= sum( ( log(means_t2(th_t2_1:th_t2_2))-mean(log(means_t2(th_t2_1:th_t2_2))) ).^2 );
SSres=sum( ( log(means_t2(th_t2_1:th_t2_2))-yfit ).^2 );
Rsq = 1-SSres/SStot;

ax = nexttile;
plot(log(NIST(th_t2_1:th_t2_2, 2)), log(NIST(th_t2_1:th_t2_2, 2)), LineWidth=2), hold on
plot(log(NIST(th_t2_1:th_t2_2, 2)), yfit, LineWidth=1.5, Color='k'), hold on
plot(log(NIST(th_t2_1:th_t2_2, 2)), log(means_t2(th_t2_1:th_t2_2)), '.', MarkerSize=30, Color='r');

text_str = sprintf('y = %.3fx + %.3f \nR^2 = %.3f', coeffs(1), coeffs(2), Rsq);
annotation('textbox', [0.575, 0.835, 0.21, 0.07], 'String', text_str, 'FontSize',13);
% annotation('textbox', [0.02, 0.95, 0.2, 0.1], 'String', text_str, ...
%            'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k');
xaxis = flipud(round(NIST(th_t2_1:th_t2_2, 2)));
xticks(flipud(log(NIST(th_t2_1:th_t2_2, 2))))
xticklabels(xaxis)
yticks(flipud(log(NIST(th_t2_1:th_t2_2, 2))))
yticklabels(xaxis)
xlim(fliplr([log(NIST(th_t2_1, 2)), log(NIST(th_t2_2, 2))]))
ylim(fliplr([log(NIST(th_t2_1, 2)), log(NIST(th_t2_2, 2))]))
title('T2 [ms]')
set(ax,'xticklabel',[])

% ------------------------ error T1s ----------------------------------
eval.err_t1(th_t1_1:th_t1_2) = (eval.means_t1(th_t1_1:th_t1_2) - NIST((th_t1_1:th_t1_2), 1))./NIST((th_t1_1:th_t1_2), 1) *1e2;
ax = nexttile;
plot(log(NIST(th_t1_1:th_t1_2, 1)),eval.err_t1(th_t1_1:th_t1_2), '.', MarkerSize=30);
title('error T1 [%]')
yline(10, '--')
yline(-10, '--')
yline(0)
xaxis = flipud(round(NIST(th_t1_1:th_t1_2, 1)));
xticks(flipud(log(NIST(th_t1_1:th_t1_2, 1))))
% xticklabels(xaxis)
xlim(fliplr([log(NIST(th_t1_1, 1)), log(NIST(th_t1_2, 1))]))
ylim([-50 50])
pbaspect([1 0.5 1])
set(ax,'xticklabel',[])

% ------------------------ error T2s ----------------------------------
eval.err_t2(th_t2_1:th_t2_2) = (eval.means_t2(th_t2_1:th_t2_2) - NIST((th_t2_1:th_t2_2), 2))./NIST((th_t2_1:th_t2_2), 2) *1e2;
ax = nexttile;
plot(log(NIST(th_t2_1:th_t2_2, 2)),eval.err_t2(th_t2_1:th_t2_2), '.', MarkerSize=30);
title('error T2 [%]')
yline(10, '--')
yline(-10, '--')
yline(0)
xaxis = flipud(round(NIST(th_t2_1:th_t2_2, 2)));
xticks(flipud(log(NIST(th_t2_1:th_t2_2, 2))))
xticklabels(xaxis)
xlim(fliplr([log(NIST(th_t2_1, 2)), log(NIST(th_t2_2, 2))]))
ylim([-50 50])
pbaspect([1 0.5 1])
set(ax,'xticklabel',[])
% ------------------------ COV T1s ----------------------------------
ax = nexttile;
plot(log(NIST(th_t1_1:th_t1_2, 1)),eval.cov_t1(th_t1_1:th_t1_2), '.', MarkerSize=30);
title('COV T1 [%]')
yline(10, '--')
yline(-10, '--')
yline(0)
xaxis = flipud(round(NIST(th_t1_1:th_t1_2, 1)));
xticks(flipud(log(NIST(th_t1_1:th_t1_2, 1))))
xticklabels(xaxis)
xlim(fliplr([log(NIST(th_t1_1, 1)), log(NIST(th_t1_2, 1))]))
ylim([-50 50])
pbaspect([1 0.5 1])

% ------------------------ COV T2s ----------------------------------
ax = nexttile;
plot(log(NIST(th_t2_1:th_t2_2, 1)),eval.cov_t2(th_t2_1:th_t2_2), '.', MarkerSize=30);
title('COV T2 [%]')
yline(10, '--')
yline(-10, '--')
yline(0)
xaxis = flipud(round(NIST(th_t2_1:th_t2_2, 1)));
xticks(flipud(log(NIST(th_t2_1:th_t2_2, 1))))
xticklabels(xaxis)
xlim(fliplr([log(NIST(th_t2_1, 1)), log(NIST(th_t2_2, 1))]))
ylim([-100 100])
pbaspect([1 0.5 1])

end