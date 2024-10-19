function plotTypes11BifDiagFullVaryPhiI_v1_play
clc;
clear all;
close all;

root_txt = 'E:\paper2_Raoul\Sim_two_neurons_Raoul\Types11BifDiagFullVaryPhiI';

figure(1); hold on;
%% Plot ING
% for i = 0:1:11
%     file_txt = strcat(root_txt, '\v5\ING\ING_Types11BifDiagFullVaryPhiI', num2str(i), '.mat');
%     load(file_txt);
% 
%     figure(1);hold on;
%     plot(OneoverThetaPhiI_lin, f1, 'r-', 'LineWidth', 4);
%     plot(OneoverThetaPhiI_lin, f2, 'g*', 'LineWidth', 4);
%     plot(OneoverThetaPhiI_lin, f3, 'b-', 'LineWidth', 4);
%     plot(OneoverThetaPhiI_lin, f4, 'k-', 'LineWidth', 4);
%     plot(OneoverThetaPhiI_lin, f5, 'y-', 'LineWidth', 4);
% 
%     figure(2)
%     subplot(3,2,1);hold on
%     plot(OneoverThetaPhiI_lin, f1, 'r*', 'LineWidth', 20);
%     subplot(3,2,2);hold on
%     plot(OneoverThetaPhiI_lin, f2, 'g*', 'LineWidth', 16);
%     subplot(3,2,3);hold on
%     plot(OneoverThetaPhiI_lin, f3, 'b*', 'LineWidth', 12);
%     subplot(3,2,4);hold on
%     plot(OneoverThetaPhiI_lin, f4, 'k*', 'LineWidth', 8);
%     subplot(3,2,5);hold on
%     plot(OneoverThetaPhiI_lin, f5, 'y*', 'LineWidth', 4);
% 
%     
%     if (i == 0)
%         OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
%         f1_total = f1;
%         f2_total = f2;
%         f3_total = f3;
%         f4_total = f4;
%         f5_total = f5;
%     else
%         OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
%         f1_total = [f1_total f1];
%         f2_total = [f2_total f2];
%         f3_total = [f3_total f3];
%         f4_total = [f4_total f4];
%         f5_total = [f5_total f5];
%     end    
% end
% 
% plot(OneoverThetaPhiI_total, f1_total, 'r*');
% % plot(OneoverThetaPhiI_total, f2_total, 'g*');
% % plot(OneoverThetaPhiI_total, f3_total, 'b*');
% % plot(OneoverThetaPhiI_total, f4_total, 'k*');
% % % plot(OneoverThetaPhiI_total, f5_total, 'y*');
% 
% % Filter
% N_filter = 10;
% i_filter = 0;
% [rows, cols] = size(OneoverThetaPhiI_total);
% for i = 1:1:cols
%     if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
%         i_filter = i_filter + 1;
%         OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
%         f1_filter(1, i_filter) = f1_total(1, i);
%         f2_filter(1, i_filter) = f2_total(1, i);
%         f3_filter(1, i_filter) = f3_total(1, i);
%         f4_filter(1, i_filter) = f4_total(1, i);
%         f5_filter(1, i_filter) = f5_total(1, i);
%     end
% end
% 
% % Combine
% [rows, cols] = size(OneoverThetaPhiI_filter);
% 
% i_comb = 0;
% OneoverThetaPhiI_comb = 0;
% f_comb = 0;
% for i = 1:1:cols
%     if (isnan(f2_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f2_filter(1, i);
%     elseif (isnan(f3_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f3_filter(1, i);
%     elseif (isnan(f4_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f4_filter(1, i);
%     end
% end
% 
% % plot(OneoverThetaPhiI_comb, f_comb, 'b-', 'LineWidth', 4);

%% Plot Analytic ING
% thetaPhiI = 1./linspace(0.46, 0.59, 100);
% 
% m_tau = 0.4;
% epsilonII = -1.00;  % From I to I
% 
% T = m_tau + thetaPhiI + log(exp(-m_tau) - (1 - exp(-thetaPhiI)).*epsilonII);
% 
% plot(1./thetaPhiI, thetaPhiI./thetaPhiI.*1./T, 'b-', 'LineWidth', 4)
% 
% 
% % 
% % % % %% Plot PING
% % % % for i = 1:1:10
% % % %     file_txt = strcat(root_txt, '\PING\PING_Types11BifDiagFullVaryPhiI', num2str(i), '.mat');
% % % %     load(file_txt);
% % % %     
% % % %     plot(OneoverThetaPhiI_lin, f1, 'r*');
% % % %     plot(OneoverThetaPhiI_lin, f2, 'r*');
% % % %     plot(OneoverThetaPhiI_lin, f3, 'r*');
% % % %     plot(OneoverThetaPhiI_lin, f4, 'r*');
% % % %     plot(OneoverThetaPhiI_lin, f5, 'r*');
% % % % end
% % 
%% Plot PINGING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v5\PINGING\PINGING_Types11BifDiagFullVaryPhiI', num2str(i), '.mat');
    load(file_txt);

    figure(1);
    plot(OneoverThetaPhiI_lin, f1, 'r*');
    plot(OneoverThetaPhiI_lin, f2, 'g*');   % Top lines
    plot(OneoverThetaPhiI_lin, f3, 'b*');   % Top lines    
    plot(OneoverThetaPhiI_lin, f4, 'k*');   % Buttom lines
    plot(OneoverThetaPhiI_lin, f5, 'y*');
    
    figure(2);
    subplot(3,2,1);hold on
    plot(OneoverThetaPhiI_lin, f1, 'r*');
    subplot(3,2,2);hold on
    plot(OneoverThetaPhiI_lin, f2, 'g*');   % Top lines
    subplot(3,2,3);hold on
    plot(OneoverThetaPhiI_lin, f3, 'b*');   % Top lines    
    subplot(3,2,4);hold on
    plot(OneoverThetaPhiI_lin, f4, 'k*');   % Buttom lines
    subplot(3,2,5);hold on
    plot(OneoverThetaPhiI_lin, f5, 'y*');

    if (i == 0)
        OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
    else
        OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
    end    

end

% plot(OneoverThetaPhiI_total, f1_total, 'r*');
plot(OneoverThetaPhiI_total, f2_total, 'g-', 'LineWidth', 4);
plot(OneoverThetaPhiI_total, f3_total, 'k-', 'LineWidth', 4);
    plot(OneoverThetaPhiI_total, f4_total, 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);
% plot(OneoverThetaPhiI_total, f5_total, 'y*');


% Filter
N_filter = 1;
i_filter = 0;
[rows, cols] = size(OneoverThetaPhiI_total);
for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
    end
end

% Combine
[rows, cols] = size(OneoverThetaPhiI_filter);

i_comb = 0;
i_comb1 = 0;
OneoverThetaPhiI_comb = 0;
OneoverThetaPhiI_comb1 = 0;
f_comb = 0;
f_comb1 = 0;
for i = 1:1:cols
    if (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    end
    
    if (isnan(f2_filter(1, i)) == 0)
        i_comb = i_comb + 1;

        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f2_filter(1, i);
    end
end

% plot(OneoverThetaPhiI_comb, f_comb, 'g-', 'LineWidth', 4);
plot(OneoverThetaPhiI_comb, f_comb, 'k-', 'LineWidth', 4);

%% Plot PING
thetaPhiI = 1./linspace(0.46, 0.59, 100);
thetaPhiE = 1./0.495;

m_tau = 0.4;
epsilonEI = -0.5;       % From I to E

m_H = -log(exp(-2.*m_tau) - (1 - exp(-thetaPhiE)).*epsilonEI);
T = 2.*m_tau + thetaPhiE - m_H;

plot(1./thetaPhiI, thetaPhiI./thetaPhiI.*1./T, 'r-', 'LineWidth', 4)

% plot([0.5312 0.6], [0.3709 0.3709], 'r-', 'LineWidth', 4);
% % plot([OneoverThetaPhiI_comb1 0.4956], [f_comb1 f_comb1(1, end)], 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);
% % % plot(OneoverThetaPhiI_comb1, f_comb1, 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);
% % % 
% make_me_pretty(gcf, ...
%     gca, 40, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% 
% set(gca,'XTick',[0.46 0.50 0.54 0.58],'XTickLabel',{'';'';'';''});
% set(gca,'YTick',[0.33 0.35 0.37 0.39 0.41],'YTickLabel',{'';'';'';''});

% maximize_a_fig(gcf);
ylim([0.33 0.42]);
xlim([0.46 0.59]);
% m_savefig('Types11BifDiagFullVaryPhiI_v2', 'eps');
end

function val = cmp(x, y, tol_eq)
% CMP Two-value comparison
%   val = cmp(x, y, tol_eq)
% Input
%   x           the first number.
%   y           the second number.
%   tol_eq      if the first and second numbers are different less than
%               tol_eq, we say that the two numbers are equal.
% Output
%   val         0   : two numbers are the same.
%               -1  : the first number is less than the second number.
%               1   : the first number is greater than the second number.

if  (abs(x-y)<tol_eq)
    val=0;
elseif (x<y)
    val=-1;
else
    val=1;
end

end