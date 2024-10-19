function plotplotTypes12BifDiagFullVarygI2I_v1
clc;
clear all;
close all;

root_txt = 'E:\paper2_Raoul\Sim_two_neurons_Raoul\Types12BifDiagFullVarygI2I';

figure(1); hold on;
%% Plot PINGING
for i = 0:1:9    
    file_txt = strcat(root_txt, '\PINGING\v0\PINGING_Types12Tau04BifDiagFullVarygI2I', num2str(i), '.mat');
    load(file_txt);
    
%     plot(epsilonII_lin, f1, 'r*');
%     plot(epsilonII_lin, f2, 'g*');
%     plot(epsilonII_lin, f3, 'b*');
%     plot(epsilonII_lin, f4, 'k*');
%     plot(epsilonII_lin, f5, 'c*');
%     plot(epsilonII_lin, f51, 'y*');
%     plot(epsilonII_lin, f15, 'o');

    if (i == 0)
        epsilonII_total = epsilonII_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
        f51_total = f51;
        f15_total = f15;
    else
        epsilonII_total = [epsilonII_total epsilonII_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
        f51_total = [f51_total f51];
        f15_total = [f15_total f15];
    end      
end

% % plot(epsilonII_total, f1_total, 'r*');
% plot(epsilonII_total, f2_total, 'g*');
% plot(epsilonII_total, f3_total, 'b*');
% % plot(epsilonII_total, f4_total, 'k*');
% % plot(epsilonII_total, f5_total, 'y*');
% plot(epsilonII_total, f51_total, 'm*');
% plot(epsilonII_total, f15_total, 'co');

% Filter
N_filter = 1;
i_filter = 0;
[rows, cols] = size(epsilonII_total);
epsilonII_filter = 0;
f1_filter = 0;
f2_filter = 0;
f3_filter = 0;
f4_filter = 0;
f5_filter = 0;
f51_filter = 0;
f15_filter = 0;

for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        epsilonII_filter(1, i_filter) = epsilonII_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
        f51_filter(1, i_filter) = f51_total(1, i);
        f15_filter(1, i_filter) = f15_total(1, i);
    end
end

% Combine
[rows, cols] = size(epsilonII_filter);

i_comb = 0;
epsilonII_comb = 0;
f_comb = 0;

i_comb1 = 0;
epsilonII_comb1 = 0;
f_comb1 = 0;
for i = 1:1:cols
    if (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        epsilonII_comb(1, i_comb) = epsilonII_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    elseif (isnan(f4_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        epsilonII_comb(1, i_comb) = epsilonII_filter(1, i);
        f_comb(1, i_comb) = f4_filter(1, i);
    elseif (isnan(f15_filter(1, i)) == 0)
        i_comb1 = i_comb1 + 1;

        epsilonII_comb1(1, i_comb1) = epsilonII_filter(1, i);
        f_comb1(1, i_comb1) = f15_filter(1, i);
    end
end

plot(epsilonII_comb, f_comb, 'g-', 'LineWidth', 4);
plot(epsilonII_comb1, f_comb1, 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);

%% Plot ING
for i = 0:1:9
    file_txt = strcat(root_txt, '\ING\v0\ING_Types12Tau04BifDiagFullVarygI2I', num2str(i), '.mat');
    load(file_txt);
    
%     plot(epsilonII_lin, f1, 'r*');
%     plot(epsilonII_lin, f2, 'g*');
%     plot(epsilonII_lin, f3, 'b*');
%     plot(epsilonII_lin, f4, 'k*');
%     plot(epsilonII_lin, f5, 'c*');
%     plot(epsilonII_lin, f51, 'y*');
%     plot(epsilonII_lin, f15, 'o');

    if (i == 0)
        epsilonII_total = epsilonII_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
        f51_total = f51;
        f15_total = f15;
    else
        epsilonII_total = [epsilonII_total epsilonII_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
        f51_total = [f51_total f51];
        f15_total = [f15_total f15];
    end    
end

% % plot(epsilonII_total, f1_total, 'r*');
% plot(epsilonII_total, f2_total, 'g*');
% plot(epsilonII_total, f3_total, 'b*');
% % plot(epsilonII_total, f4_total, 'k*');
% % plot(epsilonII_total, f5_total, 'y*');
% % plot(epsilonII_total, f51_total, 'm*');
% plot(epsilonII_total, f15_total, 'co');
% 
% Filter
N_filter = 10;
i_filter = 0;
[rows, cols] = size(epsilonII_total);
epsilonII_filter = 0;
f1_filter = 0;
f2_filter = 0;
f3_filter = 0;
f4_filter = 0;
f5_filter = 0;
f51_filter = 0;
f15_filter = 0;

for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        epsilonII_filter(1, i_filter) = epsilonII_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
        f51_filter(1, i_filter) = f51_total(1, i);
        f15_filter(1, i_filter) = f15_total(1, i);
    end
end

% Combine
[rows, cols] = size(epsilonII_filter);

i_comb = 0;
epsilonII_comb = 0;
f_comb = 0;
for i = 1:1:cols
    if (isnan(f2_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        epsilonII_comb(1, i_comb) = epsilonII_filter(1, i);
        f_comb(1, i_comb) = f2_filter(1, i);
    elseif (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        epsilonII_comb(1, i_comb) = epsilonII_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    elseif (isnan(f15_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        epsilonII_comb(1, i_comb) = epsilonII_filter(1, i);
        f_comb(1, i_comb) = f15_filter(1, i);
    end
end

plot(epsilonII_comb, f_comb, 'b-', 'LineWidth', 4);

% %% Plot Analytic ING
% m_tau = 0.4;
% 
% thetaPhiI = 1./linspace(0.48, 0.52, 100);
% thetaPhiE = 1./0.7455;
% epsilonII = -0.41514;  % From I to I
% 
% T = m_tau + thetaPhiI - thetaPhiI./pi.*atan(tan(pi./thetaPhiI.*m_tau).*exp(-2.*pi.*epsilonII./thetaPhiI));
% 
% plot(1./thetaPhiI, 1./T, 'b-', 'LineWidth', 4)
% 
% %% Plot PING
% for i = 0:1:10
%     file_txt = strcat(root_txt, '\PING\v4\PING_Types12Tau04BifDiagFullVaryPhiI', num2str(i), '.mat');
%     load(file_txt);
%     
% % %     plot(OneoverThetaPhiI_lin, f1, 'r*');
% % %     plot(OneoverThetaPhiI_lin, f2, 'r*');
% %     plot(OneoverThetaPhiI_lin, f3, 'r*');
% % %     plot(OneoverThetaPhiI_lin, f4, 'r*');
% % %     plot(OneoverThetaPhiI_lin, f5, 'r*');
% %     plot(OneoverThetaPhiI_lin, f51, 'r*');
% % %     plot(OneoverThetaPhiI_lin, f15, 'r*');
%     
% %         if (size(f51, 2) ~= 100)
% %             file_txt
% %         end
% 
%     if (i == 0)
%         OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
%         f1_total = f1;
%         f2_total = f2;
%         f3_total = f3;
%         f4_total = f4;
%         f5_total = f5;
%         f51_total = f51;
%         f15_total = f15;                
%     else
%         OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
%                     
%         f1_total = [f1_total f1];
%         f2_total = [f2_total f2];
%         f3_total = [f3_total f3];
%         f4_total = [f4_total f4];
%         f5_total = [f5_total f5];
%         f51_total = [f51_total f51];
%         f15_total = [f15_total f15];
%     end      
% end
% 
% % % % plot(OneoverThetaPhiI_total, f1_total, 'r*');
% % % plot(OneoverThetaPhiI_total, f2_total, 'g*');
% % plot(OneoverThetaPhiI_total, f3_total, 'r', 'LineWidth', 4);
% % % % plot(OneoverThetaPhiI_total, f4_total, 'k*');
% % % % plot(OneoverThetaPhiI_total, f5_total, 'y*');
% % % % plot(OneoverThetaPhiI_total, f51_total, 'm*');
% % % % plot(OneoverThetaPhiI_total, f15_total, 'c*');
% % 
% % plot([0.44 0.501], [0.6153 0.6153], 'r-', 'LineWidth', 4);
% 
%% Plot Analytic PING
thetaPhiI = 1/0.5;
thetaPhiE = 1./0.7455;

m_tau = 0.4;
epsilonEI = -0.2;       % From I to E

m_H = -log(exp(-2.*m_tau) - (1 - exp(-thetaPhiE)).*epsilonEI);
T = 2.*m_tau + thetaPhiE - m_H;

plot(epsilonII_comb, epsilonII_comb./epsilonII_comb.*1./T, 'r-', 'LineWidth', 4)

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',[-0.6 -0.5 -0.4 -0.3],'XTickLabel',{'';'';'';''});
set(gca,'YTick',[0.57 0.59 0.61 0.63 0.65],'YTickLabel',{'';'';'';'';''});

maximize_a_fig(gcf);
ylim([0.558 0.66]);
xlim([-0.62 -0.22]);
m_savefig('Types12BifDiagFullVarygI2I', 'eps');

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