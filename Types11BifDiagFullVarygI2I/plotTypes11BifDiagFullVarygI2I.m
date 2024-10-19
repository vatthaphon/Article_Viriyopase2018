function plotTypes11BifDiagFullVarygI2E
clc;
clear all;
close all;


root_txt = 'E:\paper2_Raoul\Sim_two_neurons_Raoul\Types11BifDiagFullVarygI2I';

figure(1); hold on;
%% Plot ING
for i = 0:1:9
    
    file_txt = strcat(root_txt, '\v0\ING\ING_Types11BifDiagFullVarygI2I', num2str(i), '.mat');
    load(file_txt);
    
%     plot(epsilonII_lin, f1, 'r*', 'LineWidth', 4);
%     plot(epsilonII_lin, f2, 'g*', 'LineWidth', 4);
%     plot(epsilonII_lin, f3, 'b*', 'LineWidth', 4);
%     plot(epsilonII_lin, f4, 'k*', 'LineWidth', 4);
%     plot(epsilonII_lin, f5, 'y*', 'LineWidth', 4);

    if (i == 0)
        epsilonII_total = epsilonII_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
    else
        epsilonII_total = [epsilonII_total epsilonII_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
    end    
end

% % plot(epsilonEI_total, f1_total, 'r*');
% plot(epsilonII_total, f2_total, 'g*');
% plot(epsilonII_total, f3_total, 'b*');
% plot(epsilonII_total, f4_total, 'k*');
% % plot(epsilonII_total, f5_total, 'y*');

% Filter
N_filter = 10;
i_filter = 0;
[rows, cols] = size(epsilonII_total);
for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        epsilonII_filter(1, i_filter) = epsilonII_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
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
    elseif (isnan(f4_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        epsilonII_comb(1, i_comb) = epsilonII_filter(1, i);
        f_comb(1, i_comb) = f4_filter(1, i);
    end
end

plot(epsilonII_comb, f_comb, 'b-', 'LineWidth', 4);

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
% % % %% Plot PING
% % % for i = 1:1:10
% % %     file_txt = strcat(root_txt, '\PING\PING_Types11BifDiagFullVaryPhiI', num2str(i), '.mat');
% % %     load(file_txt);
% % %     
% % %     plot(OneoverThetaPhiI_lin, f1, 'r*');
% % %     plot(OneoverThetaPhiI_lin, f2, 'r*');
% % %     plot(OneoverThetaPhiI_lin, f3, 'r*');
% % %     plot(OneoverThetaPhiI_lin, f4, 'r*');
% % %     plot(OneoverThetaPhiI_lin, f5, 'r*');
% % % end
% 
%% Plot PING
thetaPhiI = 1/0.532;
thetaPhiE = 1./0.495;

m_tau = 0.4;
epsilonEI = -0.5;       % From I to E

epsilonII = linspace(-2.0, 0, 100);       % From I to E

m_H = -log(exp(-2.*m_tau) - (1 - exp(-thetaPhiE)).*epsilonEI);
T = 2.*m_tau + thetaPhiE - m_H;

plot(epsilonII, epsilonII./epsilonII.*1./T, 'r-', 'LineWidth', 4)

%% Plot PINGING
for i = 0:1:9
    file_txt = strcat(root_txt, '\v0\PINGING\PINGING_Types11BifDiagFullVarygI2I', num2str(i), '.mat');
    load(file_txt);
    
    load(file_txt);
    
%     plot(epsilonII_lin, f1, 'ro', 'LineWidth', 4);
%     plot(epsilonII_lin, f2, 'go', 'LineWidth', 4);
%     plot(epsilonII_lin, f3, 'bo', 'LineWidth', 4);
%     plot(epsilonII_lin, f4, 'ko', 'LineWidth', 4);
%     plot(epsilonII_lin, f5, 'yo', 'LineWidth', 4);
    
    if (i == 0)
        epsilonII_total = epsilonII_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
    else
        epsilonII_total = [epsilonII_total epsilonII_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
    end
end

% % plot(epsilonII_total, f1_total, 'r*');
% plot(epsilonII_total, f2_total, 'g-', 'LineWidth', 4);
% plot(epsilonII_total, f3_total, 'k-', 'LineWidth', 4);
% plot(epsilonII_total, f4_total, 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);
% % plot(OneoverThetaPhiI_total, f5_total, 'y*');

% Filter
N_filter = 1;
i_filter = 0;
[rows, cols] = size(epsilonII_total);
for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        epsilonII_filter(1, i_filter) = epsilonII_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
    end
end

% Combine
[rows, cols] = size(epsilonII_filter);

i_comb = 0;
i_comb1 = 0;
epsilonII_comb = 0;
epsilonII_comb1 = 0;
f_comb = 0;
f_comb1 = 0;
for i = 1:1:cols
    if (isnan(f2_filter(1, i)) == 0)
        i_comb = i_comb + 1;

        epsilonII_comb(1, i_comb) = epsilonII_filter(1, i);
        f_comb(1, i_comb) = f2_filter(1, i);
    end
    
    if (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        epsilonII_comb(1, i_comb) = epsilonII_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    end
        
    if (isnan(f4_filter(1, i)) == 0)
        i_comb1 = i_comb1 + 1;

        epsilonII_comb1(1, i_comb1) = epsilonII_filter(1, i);
        f_comb1(1, i_comb1) = f4_filter(1, i);
    end    
end

plot(epsilonII_comb1, f_comb1, 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);
plot(epsilonII_comb, f_comb, 'g-', 'LineWidth', 4);

%% Plot Analytic PINGING
% m_tau = 0.4;
% 
% % sce. 2
% thetaPhiI = 1./linspace(0.5657, 0.598, 100);
% thetaPhiE = 1./0.495;
% 
% % New test, this agrees with the old test below.
% epsilonEI = -0.5;       % From I to E
% epsilonII = -1.0;  % From I to I
% epsilonIE = 0.1;       % From E to I
% 
% delta_THETA = thetaPhiE-thetaPhiI;
% 
% first_term = exp(-m_tau);
% 
% H1 = -log(exp(-m_tau) - (1 - exp(-thetaPhiI)).*epsilonII);
% second_term = exp(-H1 - delta_THETA);
% 
% H1 = -log(exp(-m_tau) - (1 - exp(-thetaPhiI)).*epsilonII);
% first_third_term = exp(-H1 - delta_THETA) - exp(-m_tau);
% 
% gamma_THETA_E_epsilonI2E = (1 - exp(-thetaPhiE)).*epsilonEI;
% gamma_THETA_I_epsilonE2I = (1 - exp(-thetaPhiI)).*epsilonIE;
% second_third_term = 4.*exp(-delta_THETA).*gamma_THETA_E_epsilonI2E.*gamma_THETA_I_epsilonE2I;
% third_term = sqrt(first_third_term.*first_third_term + second_third_term);
% 
% up1 = first_term - second_term + third_term;
% up2 = first_term - second_term - third_term;
% 
% gamma_THETA_E_epsilonI2E = (1 - exp(-thetaPhiE)).*epsilonEI;
% down = 2.*exp(-delta_THETA).*gamma_THETA_E_epsilonI2E;
% delta_phi1 = log(up1./down);
% delta_phi2 = log(up2./down);
% 
% H1 = -log(exp(-(m_tau + delta_phi1 - delta_THETA)) - (1 - exp(-thetaPhiE)).*epsilonEI);
% H2 = -log(exp(-(m_tau + delta_phi2 - delta_THETA)) - (1 - exp(-thetaPhiE)).*epsilonEI);
% T1 = m_tau + delta_phi1 + thetaPhiI - H1;
% T2 = m_tau + delta_phi2 + thetaPhiI - H2;
% 
% %% Old test
% % m_tau = 0.4;
% % epsilonEI = -0.5;       % From I to E
% % epsilonIE = 0.1;       % From E to I
% % epsilonII = -1.00;  % From I to I
% % 
% % thetaPhiE = 1./0.495;
% % thetaPhiI = 1./linspace(0.5657, 0.598, 100);
% % DeltaTheta = thetaPhiE - thetaPhiI;
% % 
% % H = -log(exp(-m_tau) - (1 - exp(-thetaPhiI)).*epsilonII);
% % E = H + DeltaTheta;
% % 
% % A = exp(-DeltaTheta).*(1 - exp(-thetaPhiE)).*epsilonEI;
% % B = exp(-E) - exp(-m_tau);
% % C = -(1 - exp(-thetaPhiI)).*epsilonIE;
% % 
% % DeltaPhi1 = log((-B + sqrt(B.*B - 4.*A.*C))./(2.*A));
% % DeltaPhi2 = log((-B - sqrt(B.*B - 4.*A.*C))./(2.*A));
% % 
% % PhiE1 = DeltaPhi1 + thetaPhiI;
% % PhiE2 = DeltaPhi2 + thetaPhiI;
% % 
% % T1 = m_tau + PhiE1 + log(exp(-(m_tau + DeltaPhi1 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI);
% % T2 = m_tau + PhiE2 + log(exp(-(m_tau + DeltaPhi2 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI);
% 
% % figure(3);hold on
% % plot(1./thetaPhiI, 1./T1, 'b*')
% % plot(1./thetaPhiI, 1./T2, 'r*')
% 
% % sce. 3
% thetaPhiI = 1./linspace(0.5226, 0.5646, 100);
% thetaPhiE = 1./0.495;
% 
% DeltaTheta = thetaPhiE - thetaPhiI;
% 
% C = (1 - exp(-thetaPhiE)).*epsilonEI;
% D = (1 - exp(-thetaPhiI)).*epsilonII;
% E = (1 - exp(-thetaPhiI)).*epsilonIE;
% 
% A = C;
% B = exp(-m_tau) - D - exp(-m_tau + DeltaTheta);
% C = -E.*exp(DeltaTheta);
% 
% DeltaPhi1 = log((-B + sqrt(B.*B - 4.*A.*C))./(2.*A));
% DeltaPhi2 = log((-B - sqrt(B.*B - 4.*A.*C))./(2.*A));
% 
% % m_tau + DeltaPhi1 - DeltaTheta
% 
% T1 = thetaPhiE + log(exp(-(m_tau + DeltaPhi1 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI) + m_tau + DeltaPhi1 - DeltaTheta;
% T2 = thetaPhiE + log(exp(-(m_tau + DeltaPhi2 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI) + m_tau + DeltaPhi2 - DeltaTheta;
% 
% % figure(3);hold on
% plot(1./thetaPhiI, 1./T1, 'b*')
% plot(1./thetaPhiI, 1./T2, 'r*')

make_me_pretty(gcf, ...
    gca, 40, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',[-1.3 -1.1 -0.9 -0.7],'XTickLabel',{'';'';'';''});
set(gca,'YTick',[0.35 0.37 0.39 0.41],'YTickLabel',{'';'';'';''});

maximize_a_fig(gcf);
ylim([0.344 0.42]);
xlim([-1.392 -0.6687]);
% m_savefig('Types11BifDiagFullVarygI2I', 'eps');
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