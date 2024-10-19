function plotTypes12BifDiagFullVaryPhiI
clc;
clear all;
close all;

root_txt = 'C:\paper2_Raoul\Sim_two_neurons_Raoul\Types12BifDiagFullVaryPhiE';
LW = 20;

scenarios_color = {'y*', 'm*', 'c*', 'r*', 'g*'};
% scenarios_color = {'y', 'm', 'c', 'r', 'g'};

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
% %     plot(OneoverThetaPhiI_lin, f51, 'b*');f
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
N_OneoverThetaPhiE_lin = 100;
thetaPhiI = 1/0.5;
for i = 1:1:10
% for i = 5
    file_txt = strcat(root_txt, '\PINGING\v4\PINGING_Types12Tau04BifDiagFullVaryPhiE', num2str(i), '.mat');
    load(file_txt);
    
    %% If there is more than one fixed point for each scenarios.
    for j = 1:1:5
        for jj = 1:1:N_OneoverThetaPhiE_lin
            if (sum(isnan(firstIte_FP_DeltaPsi(j, jj, :))) <= 3)
                'More than one fixed point'
            end
        end
        
        for k = 1:1:5
            for jj = 1:1:N_OneoverThetaPhiE_lin
            if (sum(isnan(secondIte_FP_DeltaPsi(j, k, jj, :))) <= 3)
                'More than one fixed point'
            end                
            end
        end
    end
    
    figure(1); hold on;
    for j = 1:1:5
        % One iteration
        tmp_deltaPsi = [];
        for jj = 1:1:N_OneoverThetaPhiE_lin
            tmp_deltaPsi(jj) = firstIte_FP_DeltaPsi(j, jj, 1);
            
            if (isnan(tmp_deltaPsi(jj)) == 0)
                switch(j)
                    case 1
                    case 2                       
                        tmp_f(jj) = period_scenario2(gammaE, thetaVE, 1/OneoverThetaPhiE_lin(jj), tmp_deltaPsi(jj), tau, epsilonEI);
                    case 3
                        tmp_f(jj) = period_scenario3(gammaE, thetaVE, 1/OneoverThetaPhiE_lin(jj), tmp_deltaPsi(jj), tau, epsilonEI);
                    case 4
                    case 5
                    otherwise
                        fprintf('Invalid scenarios\n' );
                end
            else
                tmp_f(jj) = NaN;
            end
        end
        
        figure(1);hold on;
        plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, j)), 'LineWidth', 7)
        figure(2);hold on;
        plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, j)), 'LineWidth', 7)
        
        % Two iterations
        for k = 1:1:5
            tmp_deltaPsi = [];
            for jj = 1:1:N_OneoverThetaPhiE_lin
                tmp_deltaPsi(jj) = secondIte_FP_DeltaPsi(j, k, jj, 1);
                
                if (isnan(tmp_deltaPsi(jj)) == 0)
                    if ((j == 5) && (k == 1))
                        tmp_f(jj) = period_scenario5then1(gammaE, gammaI, thetaVE, thetaVI, 1/OneoverThetaPhiE_lin(jj), thetaPhiI, tmp_deltaPsi(jj), tau, epsilonEI, epsilonIE);
                    end
                else
                    tmp_f(jj) = NaN;
                end
            end
        
            figure(1);hold on;
            plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, j)), 'LineWidth', 14)
            plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, k)), 'LineWidth', 3)            
            
            figure(2);hold on;
            plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, j)), 'LineWidth', 14)
            plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, k)), 'LineWidth', 3)            
            
        end
    end
    
%     plot()
%     if (i == 1)
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
end

% % plot(OneoverThetaPhiI_total, f1_total, 'r*');
% plot(OneoverThetaPhiI_total, f2_total, 'g*');
% plot(OneoverThetaPhiI_total, f3_total, 'g-', 'LineWidth', 8);
% % plot(OneoverThetaPhiI_total, f4_total, 'k*');
% % plot(OneoverThetaPhiI_total, f5_total, 'y*');
% plot(OneoverThetaPhiI_total, f51_total, 'g-', 'LineWidth', 8, 'Color', [101/255 163/255 92/255]);
% % plot(OneoverThetaPhiI_total, f15_total, 'c*');
% 
% plot([0.5004 0.501], [0.615 0.6155], 'g-', 'LineWidth', 8);

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
% xlim([0.48 0.52]);
% m_savefig('Types12Tai04BifDiagFullVaryPhiI_v4', 'eps');

end

function f = period_scenario5then1(gammaE, gammaI, thetaVE, thetaVI, thetaPhiE, thetaPhiI, delta_psi, tau, eIE, eEI)

I_I = (gammaI*thetaVI)./(1 - exp(-gammaI*thetaPhiI));   
I_E = (gammaE*thetaVE)./(1 - exp(-gammaE*thetaPhiE));   

Hsine1 = Hepsilon(gammaI, thetaPhiI + tau - delta_psi, eEI, I_I, thetaVI, thetaPhiI, 2);

HLIF = Hepsilon(gammaE, 2*tau + thetaPhiI - Hsine1, eIE, I_E, thetaVE, thetaPhiE, 1);

Hsine2 = Hepsilon(gammaI, thetaPhiI + tau - delta_psi, eEI, I_I, thetaVI, thetaPhiI, 2);

T = 2*tau +  thetaPhiE + thetaPhiI - Hsine2 - HLIF;
f = 1/T;

end


function f = period_scenario2(gammaE, thetaVE, thetaPhiE, delta_psi, tau, eIE)

I_E = (gammaE*thetaVE)./(1 - exp(-gammaE*thetaPhiE));   % Current to drive the voltage of the LIF from 0 to thetaVE within thetaPhiE

HLIF = Hepsilon(gammaE, tau + delta_psi, eIE, I_E, 1, thetaPhiE, 1);

T = tau +  delta_psi + thetaPhiE - HLIF;
f = 1/T;

end

function f = period_scenario3(gammaE, thetaVE, thetaPhiE, delta_psi, tau, eIE)

I_E = (gammaE*thetaVE)./(1 - exp(-gammaE*thetaPhiE));   % Current to drive the voltage of the LIF from 0 to thetaVE within thetaPhiE

HLIF = Hepsilon(gammaE, tau + delta_psi, eIE, I_E, 1, thetaPhiE, 1);

T = tau +  delta_psi + thetaPhiE - HLIF;
f = 1/T;

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