clear;
close all;
load ex_eusipco_2010 % includes all the signals and parameters for the AMIPAPA

[m, w] = apa_all_AEC(x, d, miu, ord, p, dlt, a, h1, 1, 1);

% Plot the misalignment over time
t = linspace(0, N / 8000, N);
% figure;
% plot(t, m, 'k', 'LineWidth', 0.5);
% xlabel('\bf\fontsize{12}Time (seconds)');
% ylabel('\bf\fontsize{12}Misalignment (dB)');
% grid on;
%plot the filter weights and echo path
figure;
% plot(w);
% hold on;
plot(h1, 'b');
legend('Echo Path');
% legend('Filter Weights', 'Echo Path');
