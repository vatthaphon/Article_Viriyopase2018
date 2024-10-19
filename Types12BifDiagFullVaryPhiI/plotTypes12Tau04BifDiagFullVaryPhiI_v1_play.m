function plotTypes12Tau04BifDiagFullVaryPhiI_v1_play
clc;
clear all;
close all;

root_txt = 'E:\paper2_Raoul\Sim_two_neurons_Raoul\Types12BifDiagFullVaryPhiI';

figure(1); hold on;
% %% Plot ING
% for i = 0:1:10
%     file_txt = strcat(root_txt, '\ING\v4\ING_Types12Tau04BifDiagFullVaryPhiI', num2str(i), '.mat');
%     load(file_txt);
%     
% % %     plot(OneoverThetaPhiI_lin, f1, 'b*');
% % %     plot(OneoverThetaPhiI_lin, f2, 'b*');
% %     plot(OneoverThetaPhiI_lin, f3, 'b*');
% % %     plot(OneoverThetaPhiI_lin, f4, 'b*');
% % %     plot(OneoverThetaPhiI_lin, f5, 'b*');
% %     plot(OneoverThetaPhiI_lin, f51, 'b*');
% % %     plot(OneoverThetaPhiI_lin, f15, 'b*');
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
% % % plot(OneoverThetaPhiI_total, f1_total, 'r*');
% % plot(OneoverThetaPhiI_total, f2_total, 'g*');
% % plot(OneoverThetaPhiI_total, f3_total, 'b*');
% % % plot(OneoverThetaPhiI_total, f4_total, 'k*');
% % % plot(OneoverThetaPhiI_total, f5_total, 'y*');
% % plot(OneoverThetaPhiI_total, f51_total, 'm*');
% % % plot(OneoverThetaPhiI_total, f15_total, 'c*');
% 
% % Filter
% N_filter = 10;
% i_filter = 0;
% [rows, cols] = size(OneoverThetaPhiI_total);
% OneoverThetaPhiI_filter = 0;
% f1_filter = 0;
% f2_filter = 0;
% f3_filter = 0;
% f4_filter = 0;
% f5_filter = 0;
% f51_filter = 0;
% f15_filter = 0;
% 
% for i = 1:1:cols
%     if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
%         i_filter = i_filter + 1;
%         OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
%         f1_filter(1, i_filter) = f1_total(1, i);
%         f2_filter(1, i_filter) = f2_total(1, i);
%         f3_filter(1, i_filter) = f3_total(1, i);
%         f4_filter(1, i_filter) = f4_total(1, i);
%         f5_filter(1, i_filter) = f5_total(1, i);
%         f51_filter(1, i_filter) = f51_total(1, i);
%         f15_filter(1, i_filter) = f15_total(1, i);
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
%     if (isnan(f3_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f3_filter(1, i);
%     elseif (isnan(f2_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f2_filter(1, i);        
%     elseif (isnan(f51_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
% 
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f51_filter(1, i);
%     end
% end
% 
% % plot(OneoverThetaPhiI_comb, f_comb, 'b-', 'LineWidth', 4);
% 
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
% %% Plot Analytic PING
% thetaPhiI = 1./linspace(0.48, 0.52, 100);
% thetaPhiE = 1./0.7455;
% 
% m_tau = 0.4;
% epsilonEI = -0.2;       % From I to E
% 
% m_H = -log(exp(-2.*m_tau) - (1 - exp(-thetaPhiE)).*epsilonEI);
% T = 2.*m_tau + thetaPhiE - m_H;
% 
% plot(1./thetaPhiI, thetaPhiI./thetaPhiI.*1./T, 'r-', 'LineWidth', 4)
% 
% 
% % % Filter
% % N_filter = 10;
% % i_filter = 0;
% % [rows, cols] = size(OneoverThetaPhiI_total);
% % OneoverThetaPhiI_filter = 0;
% % f1_filter = 0;
% % f2_filter = 0;
% % f3_filter = 0;
% % f4_filter = 0;
% % f5_filter = 0;
% % f51_filter = 0;
% % f15_filter = 0;
% % 
% % for i = 1:1:cols
% %     if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
% %         i_filter = i_filter + 1;
% %         OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
% %         f1_filter(1, i_filter) = f1_total(1, i);
% %         f2_filter(1, i_filter) = f2_total(1, i);
% %         f3_filter(1, i_filter) = f3_total(1, i);
% %         f4_filter(1, i_filter) = f4_total(1, i);
% %         f5_filter(1, i_filter) = f5_total(1, i);
% % %         f51_filter(1, i_filter) = f51_total(1, i);
% % %         f15_filter(1, i_filter) = f15_total(1, i);
% %     end
% % end
% % 
% % % Combine
% % [rows, cols] = size(OneoverThetaPhiI_filter);
% % 
% % i_comb = 0;
% % OneoverThetaPhiI_comb = 0;
% % f_comb = 0;
% % for i = 1:1:cols
% %     if (isnan(f2_filter(1, i)) == 0)
% %         i_comb = i_comb + 1;
% %         
% %         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
% %         f_comb(1, i_comb) = f2_filter(1, i);
% %     elseif (isnan(f3_filter(1, i)) == 0)
% %         i_comb = i_comb + 1;
% %         
% %         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
% %         f_comb(1, i_comb) = f3_filter(1, i);        
% %     end
% % end
% % 
% % % plot(OneoverThetaPhiI_comb, f_comb, 'r-', 'LineWidth', 4);
% % % % plot([0.501 0.5051], [0.6169 0.6169], 'r-', 'LineWidth', 4);
% % % % plot([0.5051 0.509], [0.6169 0.6175], 'r-', 'LineWidth', 4);
% % % % 
% % % % 
%% Plot PINGING
for i = 0:1:10
    file_txt = strcat(root_txt, '\PINGING\v4\PINGING_Types12Tau04BifDiagFullVaryPhiI', num2str(i), '.mat');
    load(file_txt);

    figure(1);hold on;    
    plot(OneoverThetaPhiI_lin, f1, 'r*');
    plot(OneoverThetaPhiI_lin, f2, 'g*');
    plot(OneoverThetaPhiI_lin, f3, 'b*');
    plot(OneoverThetaPhiI_lin, f4, 'k*');
    plot(OneoverThetaPhiI_lin, f5, 'y*');
    plot(OneoverThetaPhiI_lin, f51, 'r*');
    plot(OneoverThetaPhiI_lin, f15, 'r*');  

    figure(2);hold on;
    subplot(3,3,1);hold on
    plot(OneoverThetaPhiI_lin, f1, 'r*');
    subplot(3,3,2);hold on
    plot(OneoverThetaPhiI_lin, f2, 'g*');
    subplot(3,3,3);hold on
    plot(OneoverThetaPhiI_lin, f3, 'b*');
    subplot(3,3,4);hold on
    plot(OneoverThetaPhiI_lin, f4, 'k*');
    subplot(3,3,5);hold on
    plot(OneoverThetaPhiI_lin, f5, 'y*');
    subplot(3,3,6);hold on
    plot(OneoverThetaPhiI_lin, f51, 'r*');
    subplot(3,3,7);hold on
    plot(OneoverThetaPhiI_lin, f15, 'r*');  
    
    if (i == 0)
        OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
        f51_total = f51;
        f15_total = f15;
    else
        OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
        f51_total = [f51_total f51];
        f15_total = [f15_total f15];
    end      
end

% % plot(OneoverThetaPhiI_total, f1_total, 'r*');
% plot(OneoverThetaPhiI_total, f2_total, 'g*');
% plot(OneoverThetaPhiI_total, f3_total, 'g-', 'LineWidth', 4);
% % plot(OneoverThetaPhiI_total, f4_total, 'k*');
% % plot(OneoverThetaPhiI_total, f5_total, 'y*');
% plot(OneoverThetaPhiI_total, f51_total, 'g-', 'LineWidth', 4);
% % plot(OneoverThetaPhiI_total, f15_total, 'c*');

figure(1)
% plot([0.5004 0.501], [0.615 0.6155], 'g-', 'LineWidth', 4);

% % Filter
% N_filter = 1;
% i_filter = 0;
% [rows, cols] = size(OneoverThetaPhiI_total);
% OneoverThetaPhiI_filter = 0;
% f1_filter = 0;
% f2_filter = 0;
% f3_filter = 0;
% f4_filter = 0;
% f5_filter = 0;
% f51_filter = 0;
% f15_filter = 0;
% 
% for i = 1:1:cols
%     if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
%         i_filter = i_filter + 1;
%         OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
%         f1_filter(1, i_filter) = f1_total(1, i);
%         f2_filter(1, i_filter) = f2_total(1, i);
%         f3_filter(1, i_filter) = f3_total(1, i);
%         f4_filter(1, i_filter) = f4_total(1, i);
%         f5_filter(1, i_filter) = f5_total(1, i);
%         f51_filter(1, i_filter) = f51_total(1, i);
%         f15_filter(1, i_filter) = f15_total(1, i);
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
%     if (isnan(f3_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f3_filter(1, i);
%     elseif (isnan(f4_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f4_filter(1, i);
%     elseif (isnan(f5_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f5_filter(1, i);            
%     elseif (isnan(f15_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
% 
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f15_filter(1, i);
%     end
% end
% 
% plot(OneoverThetaPhiI_comb, f_comb, 'g-*', 'LineWidth', 4);

% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)

% set(gca,'XTick',[0.48 0.49 0.50 0.51 0.52],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0.59 0.61 0.63],'YTickLabel',{'';'';''});

% maximize_a_fig(gcf);
% ylim([0.585 0.645]);
xlim([0.48 0.52]);
% m_savefig('Types12Tai04BifDiagFullVaryPhiI_v2', 'eps');

%% PINGING sce3 analytic
thetaPhiI = 1./linspace(0.501, 0.52, 100);
thetaPhiE = 1./0.7455;

m_tau = 0.4;
epsilonEI = -0.2;       % From I to E
epsilonIE = 0.1;       % From E to I
epsilonII = -0.41514;  % From I to I

[rows, cols] = size(thetaPhiI);

for i = 1:1:cols
    
ThetaE = thetaPhiE;
ThetaI = thetaPhiI(1, i);
delta_Theta = ThetaE - ThetaI;

f = @(delta_phi)RHS_sce3(delta_phi, m_tau, delta_Theta, ThetaE, ThetaI, epsilonIE, epsilonEI, epsilonII);
delta_phi = fzero(f, 0.6);

T = m_tau + delta_phi + ThetaI - m_get_HLIF(m_tau + delta_phi - delta_Theta, ThetaE, epsilonEI);

f_sce3_anal(i) = 1/T;

end

plot(1./thetaPhiI, f_sce3_anal, 'g-')

%% PINGING sce5 then sce1 analytic
thetaPhiI = 1./linspace(0.48, 0.50, 100);
thetaPhiE = 1./0.7455;

m_tau = 0.4;
epsilonEI = -0.2;       % From I to E
epsilonIE = 0.1;       % From E to I
epsilonII = -0.41514;  % From I to I

[rows, cols] = size(thetaPhiI);

for i = 1:1:cols
    
ThetaE = thetaPhiE;
ThetaI = thetaPhiI(1, i);
delta_Theta = ThetaE - ThetaI;

f = @(delta_phi)RHS_sce51(delta_phi, m_tau, delta_Theta, ThetaE, ThetaI, epsilonIE, epsilonEI, epsilonII);
delta_phi = fzero(f, 0.6);

T = 2.*m_tau + ThetaI - m_get_Hsine_second_part(ThetaE - delta_phi + m_tau, ThetaI, epsilonIE) + ThetaE...
    - m_get_HLIF(2.*m_tau + ThetaI - m_get_Hsine_second_part(ThetaE - delta_phi + m_tau, ThetaI, epsilonIE), ThetaE, epsilonEI);

f_sce51_anal(i) = 1/T;

end

plot(1./thetaPhiI, f_sce51_anal, 'k-')

end

function val = RHS_sce51(delta_phi, tau, delta_Theta, ThetaE, ThetaI, epsilonE2I, epsilonI2E, epsilonI2I)

val = delta_phi + m_get_Hsine_second_part(tau, ThetaI, epsilonI2I) ...
    - m_get_HLIF(2.*tau + ThetaI - m_get_Hsine_second_part(ThetaE - delta_phi + tau, ThetaI, epsilonE2I), ThetaE, epsilonI2E);

end


function val = RHS_sce3(delta_phi, tau, delta_Theta, ThetaE, ThetaI, epsilonE2I, epsilonI2E, epsilonI2I)

tmp = m_get_Hsine(tau - delta_phi + delta_Theta, ThetaI, epsilonE2I);

val = delta_phi - m_get_HLIF(tau + delta_phi - delta_Theta, ThetaE, epsilonI2E) + m_get_Hsine(tmp + delta_phi - delta_Theta, ThetaI, epsilonI2I);

end

function val = m_get_Hsine(phi, Theta, epsilon)

val = (Theta./pi).*atan(tan(pi./Theta.*phi).*exp(-2.*pi.*epsilon./Theta));

end

function val = m_get_Hsine_second_part(phi, Theta, epsilon)

val = (Theta./pi).*atan(tan(pi./Theta.*phi).*exp(-2.*pi.*epsilon./Theta)) + Theta;

end

function val = m_get_HLIF(phi, Theta, epsilon)

val = -log(exp(-phi) - (1 - exp(-Theta)).*epsilon);

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