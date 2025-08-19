function [means_t1, means_t2]=plot_vs_NIST(NIST, t1, t2, rois)

n = size(rois, 2);

means_t1 = zeros(n, 1);
means_t2 = zeros(n, 1);
for ii = 1:n
    means_t1(ii) = mean(t1(rois{ii}));
    means_t2(ii) = mean(t2(rois{ii}));
end

% % % figure; 
% % % tiledlayout(1, 2);
% % % 
% % % nexttile;
% % % hold on;
% % % plot(NIST(1:n, 1), '.', MarkerSize=10);
% % % plot(means_t1, '.', MarkerSize=10);
% % % legend('NIST', 'MRF')
% % % title('T1 [ms]')
% % % 
% % % nexttile;
% % % hold on;
% % % plot(NIST(1:n, 2), '.', MarkerSize=10);
% % % plot(means_t2, '.', MarkerSize=10);
% % % legend('NIST', 'MRF')
% % % title('T2 [ms]')
% % % 
% % % figure;
% % % tiledlayout(1, 2);
% % % 
% % % nexttile;
% % % hold on;
% % % plot(NIST(1:n, 1), NIST(1:n, 1));
% % % plot(NIST(1:n, 1), means_t1, '.', MarkerSize=10);
% % % xlabel('NIST');
% % % ylabel('MRF');
% % % set(gca, 'XScale', 'log');
% % % set(gca, 'YScale', 'log');
% % % title('log(T1)')
% % % 
% % % nexttile;
% % % hold on;
% % % plot(NIST(1:n, 2), NIST(1:n, 2));
% % % plot(NIST(1:n, 2), means_t2, '.', MarkerSize=10);
% % % xlabel('NIST');
% % % ylabel('MRF');
% % % set(gca, 'XScale', 'log');
% % % set(gca, 'YScale', 'log');
% % % title('log(T2)')
% % % 
% % % 
% % % figure;
% % % tiledlayout(1, 2);
% % % 
% % % nexttile;
% % % plot(means_t1./NIST(1:n, 1), '.', MarkerSize=10);
% % % ylim([0 2]);
% % % title('T1'); 
% % % refline(0, 1);
% % % 
% % % nexttile; 
% % % plot(means_t2./NIST(1:n, 2), '.', MarkerSize=10);
% % % ylim([0 2]);
% % % title('T2');
% % % refline(0, 1);
% % % 
% % % 
coeffs = polyfit(log(NIST(1:n, 2)), log(means_t2), 1);
yfit = polyval(coeffs,log(NIST(1:n, 2)));
SStot= sum( ( log(means_t2)-mean(log(means_t2)) ).^2 );
SSres=sum( ( log(means_t2)-yfit ).^2 );
Rsq = 1-SSres/SStot;


figure; 
subplot(321)
plot(flipud(NIST(1:n, 1)), '.', MarkerSize=10), hold on
plot(flipud(means_t1), '.', MarkerSize=10);
legend('NIST', 'MRF', 'Location', 'northwest')
title('T1 [ms]'), grid on
subplot(322)
plot(flipud(NIST(1:n, 2)), '.', MarkerSize=10), hold on
plot(flipud(means_t2), '.', MarkerSize=10);
legend('NIST', 'MRF', 'Location', 'northwest')
title('T2 [ms]'), grid on

subplot(323)
plot(NIST(1:n, 1), NIST(1:n, 1)), hold on
plot(NIST(1:n, 1), means_t1, '.', MarkerSize=10)
xlabel('NIST');
ylabel('MRF');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
title('log(T1)'), grid on



subplot(324)
plot(log(NIST(1:n, 2)), log(NIST(1:n, 2))), hold on
plot(log(NIST(1:n, 2)), log(means_t2), '.', MarkerSize=10);

text_str = ['y=' coeffs(1) 'x+' coeffs(2)];
annotation('textbox', [0.02, 0.95, 0.2, 0.1], 'String', text_str);

xlabel('NIST');
ylabel('MRF');
xticks(flipud(log(NIST(1:n, 2))))
xticklabels(flipud(round(NIST(1:n, 2))))
yticklabels('')


% text(1, max(y), text_str, 'VerticalAlignment', 'top', ...
%      'HorizontalAlignment', 'left', 'BackgroundColor', 'w', ...
%      'EdgeColor', 'k');

% title('log(T2)'), grid on

subplot(325)
plot(flipud((means_t1-NIST(1:n, 1))./NIST(1:n, 1)), '.', MarkerSize=10);
ylim([-1;1])
title('error T1 [%]'), grid on
subplot(326)
plot(flipud((means_t2-NIST(1:n, 2))./NIST(1:n, 2)), '.', MarkerSize=10);
ylim([-1;1])
title('error T2 [%]'), grid on


end