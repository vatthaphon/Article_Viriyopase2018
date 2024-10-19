function plotTypes12Tau04BifDiagFullVaryPhiE_play
clc;
clear all;
close all;

root_txt = 'E:\paper2_Raoul\Sim_two_neurons_Raoul\Types12BifDiagFullVaryPhiE';

figure(1); hold on;
%% Plot ING
% for i = 1:1:11
%     file_txt = strcat(root_txt, '\ING\v2\ING_Types12Tau04BifDiagFullVaryPhiE', num2str(i), '.mat');
%     load(file_txt);
%     
% %     plot(OneoverThetaPhiE_lin, f1, 'b*');
% %     plot(OneoverThetaPhiE_lin, f2, 'b*');
% %     plot(OneoverThetaPhiE_lin, f3, 'b*');
% %     plot(OneoverThetaPhiE_lin, f4, 'b*');
% %     plot(OneoverThetaPhiE_lin, f5, 'b*');
% %     plot(OneoverThetaPhiE_lin, f51, 'b*');
% %     plot(OneoverThetaPhiE_lin, f15, 'b*');
%     
%     if (i == 1)
%         OneoverThetaPhiE_total = OneoverThetaPhiE_lin;
%         f1_total = f1;
%         f2_total = f2;
%         f3_total = f3;
%         f4_total = f4;
%         f5_total = f5;
%         f51_total = f51;
%         f15_total = f15;
%     else
%         OneoverThetaPhiE_total = [OneoverThetaPhiE_total OneoverThetaPhiE_lin];
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
% % % plot(OneoverThetaPhiE_total, f1_total, 'r*');
% % plot(OneoverThetaPhiE_total, f2_total, 'g*');
% % plot(OneoverThetaPhiE_total, f3_total, 'b*');
% % % plot(OneoverThetaPhiE_total, f4_total, 'k*');
% % % plot(OneoverThetaPhiE_total, f5_total, 'y*');
% % plot(OneoverThetaPhiE_total, f51_total, 'm*');
% % % plot(OneoverThetaPhiE_total, f15_total, 'c*');
% 
% % Filter
% N_filter = 10;
% i_filter = 0;
% [rows, cols] = size(OneoverThetaPhiE_total);
% for i = 1:1:cols
%     if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
%         i_filter = i_filter + 1;
%         OneoverThetaPhiE_filter(1, i_filter) = OneoverThetaPhiE_total(1, i);
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
% [rows, cols] = size(OneoverThetaPhiE_filter);
% 
% i_comb = 0;
% OneoverThetaPhiE_comb = 0;
% f_comb = 0;
% for i = 1:1:cols
%     if (isnan(f3_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
%         f_comb(1, i_comb) = f3_filter(1, i);
%     elseif (isnan(f2_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
% 
%         OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
%         f_comb(1, i_comb) = f2_filter(1, i);
%     elseif (isnan(f51_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
% 
%         OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
%         f_comb(1, i_comb) = f51_filter(1, i);        
%     end
% end
% 
% % plot(OneoverThetaPhiE_comb, f_comb, 'b-', 'LineWidth', 4);
% 
%% Plot Analytic ING
m_tau = 0.4;

thetaPhiE = 1./linspace(0.71, 0.78, 100);
thetaPhiI = 1/0.5;   % Phase threshold
epsilonII = -0.41514;  % From I to I

T = m_tau + thetaPhiI - thetaPhiI./pi.*atan(tan(pi./thetaPhiI.*m_tau).*exp(-2.*pi.*epsilonII./thetaPhiI));

plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*1./T, 'b-', 'LineWidth', 4)

%% Plot PING
% for i = 1:1:11
%     file_txt = strcat(root_txt, '\PING\v2\PING_Types12Tau04BifDiagFullVaryPhiE', num2str(i), '.mat');
%     load(file_txt);
%     
% % %     plot(OneoverThetaPhiE_lin, f1, 'r*');
% % %     plot(OneoverThetaPhiE_lin, f2, 'r*');
% %     plot(OneoverThetaPhiE_lin, f3, 'r*');
% % %     plot(OneoverThetaPhiE_lin, f4, 'r*');
% % %     plot(OneoverThetaPhiE_lin, f5, 'r*');
% %     plot(OneoverThetaPhiE_lin, f51, 'r*');
% % %     plot(OneoverThetaPhiE_lin, f15, 'r*');   
% 
%     if (i == 1)
%         OneoverThetaPhiE_total = OneoverThetaPhiE_lin;
%         f1_total = f1;
%         f2_total = f2;
%         f3_total = f3;
%         f4_total = f4;
%         f5_total = f5;
%         f51_total = f51;
%         f15_total = f15;
%     else
%         OneoverThetaPhiE_total = [OneoverThetaPhiE_total OneoverThetaPhiE_lin];
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
% % % plot(OneoverThetaPhiE_total, f1_total, 'r*');
% % plot(OneoverThetaPhiE_total, f2_total, 'g*');
% % % plot(OneoverThetaPhiE_total, f3_total, 'b*');
% % % plot(OneoverThetaPhiE_total, f4_total, 'k*');
% % % plot(OneoverThetaPhiE_total, f5_total, 'y*');
% % % plot(OneoverThetaPhiE_total, f51_total, 'm*');
% % % plot(OneoverThetaPhiE_total, f15_total, 'c*');
% 
% % Filter
% N_filter = 10;
% i_filter = 0;
% [rows, cols] = size(OneoverThetaPhiE_total);
% OneoverThetaPhiE_filter = 0;
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
%         OneoverThetaPhiE_filter(1, i_filter) = OneoverThetaPhiE_total(1, i);
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
% [rows, cols] = size(OneoverThetaPhiE_filter);
% 
% i_comb = 0;
% OneoverThetaPhiE_comb = 0;
% f_comb = 0;
% for i = 1:1:cols
%     if (isnan(f2_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
%         f_comb(1, i_comb) = f2_filter(1, i);
%     end
% end
% 
% % plot(OneoverThetaPhiE_comb, f_comb, 'r-x', 'LineWidth', 4);
% 
%% Plot Analytic PING
m_tau = 0.4;

thetaPhiE = 1./linspace(0.71, 0.78, 100);
epsilonEI = -0.2;       % From I to E

m_H = -log(exp(-2.*m_tau) - (1 - exp(-thetaPhiE)).*epsilonEI);
T = 2.*m_tau + thetaPhiE - m_H;

plot(1./thetaPhiE, 1./T, 'r-', 'LineWidth', 4)

%% Plot PINGING
for i = 1:1:11
    file_txt = strcat(root_txt, '\PINGING\v2\PINGING_Types12Tau04BifDiagFullVaryPhiE', num2str(i), '.mat');
    load(file_txt);
    
% %     plot(OneoverThetaPhiE_lin, f1, 'g*');
% %     plot(OneoverThetaPhiE_lin, f2, 'g*');
%     plot(OneoverThetaPhiE_lin, f3, 'g*');
% %     plot(OneoverThetaPhiE_lin, f4, 'g*');
% %     plot(OneoverThetaPhiE_lin, f5, 'g*');
%     plot(OneoverThetaPhiE_lin, f51, 'g*');
% %     plot(OneoverThetaPhiE_lin, f15, 'g*');    


    figure(1);hold on;    
    plot(OneoverThetaPhiE_lin, f1, 'r*');
    plot(OneoverThetaPhiE_lin, f2, 'g*');
    plot(OneoverThetaPhiE_lin, f3, 'b*');
    plot(OneoverThetaPhiE_lin, f4, 'k*');
    plot(OneoverThetaPhiE_lin, f5, 'y*');
    plot(OneoverThetaPhiE_lin, f51, 'r*');
    plot(OneoverThetaPhiE_lin, f15, 'r*');  

    figure(2);hold on;
    subplot(3,3,1);hold on
    plot(OneoverThetaPhiE_lin, f1, 'r*');
    subplot(3,3,2);hold on
    plot(OneoverThetaPhiE_lin, f2, 'g*');
    subplot(3,3,3);hold on
    plot(OneoverThetaPhiE_lin, f3, 'b*');
    subplot(3,3,4);hold on
    plot(OneoverThetaPhiE_lin, f4, 'k*');
    subplot(3,3,5);hold on
    plot(OneoverThetaPhiE_lin, f5, 'y*');
    subplot(3,3,6);hold on
    plot(OneoverThetaPhiE_lin, f51, 'r*');
    subplot(3,3,7);hold on
    plot(OneoverThetaPhiE_lin, f15, 'r*');  

    if (i == 1)
        OneoverThetaPhiE_total = OneoverThetaPhiE_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
        f51_total = f51;
        f15_total = f15;
    else
        OneoverThetaPhiE_total = [OneoverThetaPhiE_total OneoverThetaPhiE_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
        f51_total = [f51_total f51];
        f15_total = [f15_total f15];
    end    
end

% % plot(OneoverThetaPhiE_total, f1_total, 'r*');
% plot(OneoverThetaPhiE_total, f2_total, 'g*');
% plot(OneoverThetaPhiE_total, f3_total, 'g-', 'LineWidth', 4);
% % plot(OneoverThetaPhiE_total, f4_total, 'k*');
% % plot(OneoverThetaPhiE_total, f5_total, 'y*');
% plot(OneoverThetaPhiE_total, f51_total, 'g-', 'LineWidth', 4);
% % plot(OneoverThetaPhiE_total, f15_total, 'c*');
% 
% plot([0.7436 0.7455], [0.6141 0.6148], 'g-', 'LineWidth', 4)

% % Filter
% N_filter = 10;
% i_filter = 0;
% [rows, cols] = size(OneoverThetaPhiE_total);
% OneoverThetaPhiE_filter = 0;
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
%         OneoverThetaPhiE_filter(1, i_filter) = OneoverThetaPhiE_total(1, i);
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
% [rows, cols] = size(OneoverThetaPhiE_filter);
% 
% i_comb = 0;
% OneoverThetaPhiE_comb = 0;
% f_comb = 0;
% for i = 1:1:cols
%     if (isnan(f3_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
%         f_comb(1, i_comb) = f3_filter(1, i);
%     elseif (isnan(f2_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
%         f_comb(1, i_comb) = f2_filter(1, i);        
%     elseif (isnan(f51_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
% 
%         OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
%         f_comb(1, i_comb) = f51_filter(1, i);
%     end
% end
% 
% plot(OneoverThetaPhiE_comb, f_comb, 'g-', 'LineWidth', 4);

% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% 
% set(gca,'XTick',[0.71 0.73 0.75 0.77],'XTickLabel',{'';'';'';''});
% set(gca,'YTick',[0.59 0.61 0.63],'YTickLabel',{'';'';''});
% 
% ylim([0.585 0.645]);
figure(1)
xlim([0.71 0.78]);
% 
% maximize_a_fig(gcf);
% m_savefig('Types12Tau04BifDiagFullVaryPhiE_v2', 'eps');

%% PINGING sce3 analytic
thetaPhiE = 1./linspace(0.71, 0.7436, 100);
thetaPhiI = 1./0.5;

m_tau = 0.4;
epsilonEI = -0.2;       % From I to E
epsilonIE = 0.1;       % From E to I
epsilonII = -0.41514;  % From I to I

[rows, cols] = size(thetaPhiE);

for i = 1:1:cols
    
ThetaI = thetaPhiI;
ThetaE = thetaPhiE(1, i);
delta_Theta = ThetaE - ThetaI;

f = @(delta_phi)RHS_sce3(delta_phi, m_tau, delta_Theta, ThetaE, ThetaI, epsilonIE, epsilonEI, epsilonII);
delta_phi = fzero(f, 0.61);

T = m_tau + delta_phi + ThetaI - m_get_HLIF(m_tau + delta_phi - delta_Theta, ThetaE, epsilonEI);

f_sce3_anal(i) = 1/T;

end

plot(1./thetaPhiE, f_sce3_anal, 'g-')

%% PINGING sce5 then sce1 analytic
thetaPhiE = 1./linspace(0.7445, 0.78, 100);
thetaPhiI = 1./0.5;

m_tau = 0.4;
epsilonEI = -0.2;       % From I to E
epsilonIE = 0.1;       % From E to I
epsilonII = -0.41514;  % From I to I

[rows, cols] = size(thetaPhiE);

for i = 1:1:cols
    
ThetaI = thetaPhiI;
ThetaE = thetaPhiE(1, i);
delta_Theta = ThetaE - ThetaI;

f = @(delta_phi)RHS_sce51(delta_phi, m_tau, delta_Theta, ThetaE, ThetaI, epsilonIE, epsilonEI, epsilonII);
delta_phi = fzero(f, 0.62);

T = 2.*m_tau + ThetaI - m_get_Hsine_second_part(ThetaE - delta_phi + m_tau, ThetaI, epsilonIE) + ThetaE...
    - m_get_HLIF(2.*m_tau + ThetaI - m_get_Hsine_second_part(ThetaE - delta_phi + m_tau, ThetaI, epsilonIE), ThetaE, epsilonEI);

f_sce51_anal(i) = 1/T;

end

plot(1./thetaPhiE, f_sce51_anal, 'k-')
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