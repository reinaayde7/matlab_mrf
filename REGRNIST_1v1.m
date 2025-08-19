function [eval1, eval2]=REGRNIST_1v1(NIST, eval1, eval2, th, labels)

set(0,'defaultfigurecolor','w');
set(0,'defaultaxesfontname','arial');
set(0,'defaultaxesfontsize',13);

th1 = th(1); th2 = th(2);
% n = size(rois, 2);

% means_t1 = zeros(n, 1);
% means_t2 = zeros(n, 1);
% for ii = 1:n
%     means_t1(ii) = mean(t1(rois{ii}));
%     means_t2(ii) = mean(t2(rois{ii}));
% end

means_t1_1 = eval1.means_t1;
means_t2_1 = eval1.means_t2;

means_t1_2 = eval2.means_t1;
means_t2_2 = eval2.means_t2;
% ------------------------ T1s ----------------------------------

%regression
coeffs1 = polyfit(log(NIST(th1:th2, 1)), log(means_t1_1(th1:th2)), 1);
yfit1 = polyval(coeffs1,log(NIST(th1:th2, 1)));
SStot= sum( ( log(means_t1_1(th1:th2))-mean(log(means_t1_1(th1:th2))) ).^2 );
SSres=sum( ( log(means_t1_1(th1:th2))-yfit1 ).^2 );
Rsq1 = 1-SSres/SStot;

coeffs2 = polyfit(log(NIST(th1:th2, 1)), log(means_t1_2(th1:th2)), 1);
yfit2 = polyval(coeffs2,log(NIST(th1:th2, 1)));
SStot= sum( ( log(means_t1_2(th1:th2))-mean(log(means_t1_2(th1:th2))) ).^2 );
SSres=sum( ( log(means_t1_2(th1:th2))-yfit2 ).^2 );
Rsq2 = 1-SSres/SStot;

f=figure;
f.Position = [0 0 800 1000];
t = tiledlayout(3, 2, "TileSpacing", "tight");

ax = nexttile;

plot(log(NIST(th1:th2, 1)), log(NIST(th1:th2, 1)), LineWidth=2), hold on
plot(log(NIST(th1:th2, 1)), yfit1, LineWidth=1.5, Color='r'), hold on
plot(log(NIST(th1:th2, 1)), log(means_t1_1(th1:th2)), '.', MarkerSize=30, Color='r'), hold on

plot(log(NIST(th1:th2, 1)), yfit2, LineWidth=1.5, Color='b'), hold on
plot(log(NIST(th1:th2, 1)), log(means_t1_2(th1:th2)), '.', MarkerSize=30, Color='b');

text_str = sprintf('y = %.3fx + %.3f \nR^2 = %.3f', coeffs1(1), coeffs1(2), Rsq1);
annotation('textbox', [0.14, 0.86, 0.21, 0.06], 'String', text_str, 'FontSize',13, 'EdgeColor','r');
text_str = sprintf('y = %.3fx + %.3f \nR^2 = %.3f', coeffs2(1), coeffs2(2), Rsq2);
annotation('textbox', [0.26, 0.665, 0.21, 0.06], 'String', text_str, 'FontSize',13, 'EdgeColor','b');
% annotation('textbox', [0.02, 0.95, 0.2, 0.1], 'String', text_str, ...
%            'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k');
xaxis = flipud(round(NIST(th1:th2, 1)));
xticks(flipud(log(NIST(th1:th2, 1))))
xticklabels(xaxis)
yticks(flipud(log(NIST(th1:th2, 1))))
yticklabels(xaxis)
xlim(fliplr([log(NIST(th1, 1)), log(NIST(th2, 1))]))
ylim(fliplr([log(NIST(th1, 1)), log(NIST(th2, 1))]))
title('T1 [ms]')
set(ax,'xticklabel',[])


% ------------------------ T2s ----------------------------------
coeffs1 = polyfit(log(NIST(th1:th2, 2)), log(means_t2_1(th1:th2)), 1);
yfit1 = polyval(coeffs1,log(NIST(th1:th2, 2)));
SStot= sum( ( log(means_t2_1(th1:th2))-mean(log(means_t2_1(th1:th2))) ).^2 );
SSres=sum( ( log(means_t2_1(th1:th2))-yfit1 ).^2 );
Rsq1 = 1-SSres/SStot;

coeffs2 = polyfit(log(NIST(th1:th2, 2)), log(means_t2_2(th1:th2)), 1);
yfit2 = polyval(coeffs2,log(NIST(th1:th2, 2)));
SStot= sum( ( log(means_t2_2(th1:th2))-mean(log(means_t2_2(th1:th2))) ).^2 );
SSres=sum( ( log(means_t2_2(th1:th2))-yfit1 ).^2 );
Rsq2 = 1-SSres/SStot;

ax = nexttile;
plot(log(NIST(th1:th2, 2)), log(NIST(th1:th2, 2)), LineWidth=2), hold on
plot(log(NIST(th1:th2, 2)), yfit1, LineWidth=1.5, Color='r'), hold on
plot(log(NIST(th1:th2, 2)), log(means_t2_1(th1:th2)), '.', MarkerSize=30, Color='r'), hold on

plot(log(NIST(th1:th2, 2)), yfit2, LineWidth=1.5, Color='b'), hold on
plot(log(NIST(th1:th2, 2)), log(means_t2_2(th1:th2)), '.', MarkerSize=30, Color='b');

text_str = sprintf('y = %.3fx + %.3f \nR^2 = %.3f', coeffs1(1), coeffs1(2), Rsq1);
annotation('textbox', [0.565, 0.86, 0.21, 0.06], 'String', text_str, 'FontSize',13, 'EdgeColor','r');
text_str = sprintf('y = %.3fx + %.3f \nR^2 = %.3f', coeffs2(1), coeffs2(2), Rsq2);
annotation('textbox', [0.685, 0.665, 0.21, 0.06], 'String', text_str, 'FontSize',13, 'EdgeColor','b');


% annotation('textbox', [0.02, 0.95, 0.2, 0.1], 'String', text_str, ...
%            'FitBoxToText', 'on', 'BackgroundColor', 'w', 'EdgeColor', 'k');
xaxis = flipud(round(NIST(th1:th2, 2)));
xticks(flipud(log(NIST(th1:th2, 2))))
xticklabels(xaxis)
yticks(flipud(log(NIST(th1:th2, 2))))
yticklabels(xaxis)
xlim(fliplr([log(NIST(th1, 2)), log(NIST(th2, 2))]))
ylim(fliplr([log(NIST(th1, 2)), log(NIST(th2, 2))]))
title('T2 [ms]')
set(ax,'xticklabel',[])

% ------------------------ error T1s ----------------------------------
ax = nexttile;
eval1.err_t1(th1:th2) = (eval1.means_t1(th1:th2) - NIST((th1:th2), 1))./NIST((th1:th2), 1) *1e2;
eval2.err_t1(th1:th2) = (eval2.means_t1(th1:th2) - NIST((th1:th2), 1))./NIST((th1:th2), 1) *1e2;
plot(log(NIST(th1:th2, 1)),eval1.err_t1(th1:th2), '.', MarkerSize=30, Color='r'), hold on
plot(log(NIST(th1:th2, 1)),eval2.err_t1(th1:th2), '.', MarkerSize=30, Color='b');
title('error T1 [%]')
yline(10, '--')
yline(-10, '--')
yline(0)
xaxis = flipud(round(NIST(th1:th2, 1)));
xticks(flipud(log(NIST(th1:th2, 1))))
% xticklabels(xaxis)
xlim(fliplr([log(NIST(th1, 1)), log(NIST(th2, 1))]))
ylim([-30 30])
pbaspect([1 0.5 1])
set(ax,'xticklabel',[])

% ------------------------ error T2s ----------------------------------
eval1.err_t2(th1:th2) = (eval1.means_t2(th1:th2) - NIST((th1:th2), 2))./NIST((th1:th2), 2) *1e2;
eval2.err_t2(th1:th2) = (eval2.means_t2(th1:th2) - NIST((th1:th2), 2))./NIST((th1:th2), 2) *1e2;
ax = nexttile;
plot(log(NIST(th1:th2, 2)),eval1.err_t2(th1:th2), '.', MarkerSize=30, Color='r'), hold on
plot(log(NIST(th1:th2, 2)),eval2.err_t2(th1:th2), '.', MarkerSize=30, Color='b');
title('error T2 [%]')
yline(10, '--')
yline(-10, '--')
yline(0)
xaxis = flipud(round(NIST(th1:th2, 2)));
xticks(flipud(log(NIST(th1:th2, 2))))
xticklabels(xaxis)
xlim(fliplr([log(NIST(th1, 2)), log(NIST(th2, 2))]))
ylim([-30 30])
pbaspect([1 0.5 1])
set(ax,'xticklabel',[])
% ------------------------ COV T1s ----------------------------------
ax = nexttile;
plot(log(NIST(th1:th2, 1)),eval1.cov_t1(th1:th2), '.', MarkerSize=30, Color='r'), hold on 
plot(log(NIST(th1:th2, 1)),eval2.cov_t1(th1:th2), '.', MarkerSize=30, Color='b'), 
title('COV T1 [%]')
yline(10, '--')
yline(0)
xaxis = flipud(round(NIST(th1:th2, 1)));
xticks(flipud(log(NIST(th1:th2, 1))))
xticklabels(xaxis)
xtickangle(90)
xlim(fliplr([log(NIST(th1, 1)), log(NIST(th2, 1))]))
ylim([0 30])
pbaspect([1 0.5 1])

% ------------------------ COV T2s ----------------------------------
ax = nexttile;
plot(log(NIST(th1:th2, 2)),eval1.cov_t2(th1:th2), '.', MarkerSize=30, Color='r'), hold on
plot(log(NIST(th1:th2, 2)),eval2.cov_t2(th1:th2), '.', MarkerSize=30, Color='b');
title('COV T2 [%]')
yline(10, '--')
yline(0)
xaxis = flipud(round(NIST(th1:th2, 2)));
xticks(flipud(log(NIST(th1:th2, 2))))
xticklabels(xaxis)
xtickangle(90)
xlim(fliplr([log(NIST(th1, 2)), log(NIST(th2, 2))]))
ylim([0 30])
pbaspect([1 0.5 1])
legend(labels{1}, labels{2})

figure,
subplot(221)
boxplot(abs([eval1.err_t1(th1:th2),eval2.err_t1(th1:th2)]), 'Colors', 'rb'), title('|ERR| T1 boxplot')
subplot(222)
boxplot(abs([eval1.err_t2(th1:th2),eval2.err_t2(th1:th2)]), 'Colors', 'rb'), title('|ERR| T2 boxplot')
subplot(223)
boxplot([eval1.cov_t1(th1:th2),eval2.cov_t1(th1:th2)], 'Colors', 'rb'), title('COV T1 boxplot')
subplot(224)
boxplot([eval1.cov_t2(th1:th2),eval2.cov_t2(th1:th2)], 'Colors', 'rb'), title('COV T2 boxplot')

end