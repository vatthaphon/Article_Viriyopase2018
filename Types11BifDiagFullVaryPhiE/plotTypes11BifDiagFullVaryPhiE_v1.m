function plotTypes11BifDiagFullVaryPhiE
clc;
clear all;
close all;

m_tau = 0.4;
root_txt = 'C:\paper2_Raoul\Sim_two_neurons_Raoul\Types11BifDiagFullVaryPhiE';

LWSize = 40;

figure(1); hold on;
%% Plot ING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v5\ING\ING_Types11BifDiagFullVaryPhiE', num2str(i), '.mat');
    load(file_txt);
    
%     plot(OneoverThetaPhiE_lin, f1, 'b*');
%     plot(OneoverThetaPhiE_lin, f2, 'b*');
%     plot(OneoverThetaPhiE_lin, f3, 'b*');
%     plot(OneoverThetaPhiE_lin, f4, 'b*');
%     plot(OneoverThetaPhiE_lin, f5, 'b*')
    
    if (i == 0)
        OneoverThetaPhiE_total = OneoverThetaPhiE_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
    else
        OneoverThetaPhiE_total = [OneoverThetaPhiE_total OneoverThetaPhiE_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
    end
end

% % plot(OneoverThetaPhiE_total, f1_total, 'r*');
% plot(OneoverThetaPhiE_total, f2_total, 'g*');
% plot(OneoverThetaPhiE_total, f3_total, 'b*');
% plot(OneoverThetaPhiE_total, f4_total, 'k*');
% % plot(OneoverThetaPhiE_total, f5_total, 'y*');

% Filter
N_filter = 10;
i_filter = 0;
[rows, cols] = size(OneoverThetaPhiE_total);
for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        OneoverThetaPhiE_filter(1, i_filter) = OneoverThetaPhiE_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
    end
end

% Combine
[rows, cols] = size(OneoverThetaPhiE_filter);

i_comb = 0;
OneoverThetaPhiE_comb = 0;
f_comb = 0;
for i = 1:1:cols
    if (isnan(f2_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
        f_comb(1, i_comb) = f2_filter(1, i);    
    elseif (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    elseif (isnan(f4_filter(1, i)) == 0)
        i_comb = i_comb + 1;

        OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
        f_comb(1, i_comb) = f4_filter(1, i);
    end
end

% plot(OneoverThetaPhiE_comb, f_comb, 'b-', 'LineWidth', 4);

%% Plot Analytic ING
thetaPhiE = 1./linspace(0.42, 0.52, 100);
thetaPhiI = 1/0.495;
epsilonII = -1.00;  % From I to I

T = m_tau + thetaPhiI + log(exp(-m_tau) - (1 - exp(-thetaPhiI)).*epsilonII);
plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*1./T, 'b-', 'LineWidth', LWSize, 'Color', [102/255 102/255 255/255]);
% plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*1./T, 'b-', 'LineWidth', LWSize, 'Color', [153/255 153/255 255/255]);

%% Plot PING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v5\PING\PING_Types11BifDiagFullVaryPhiE', num2str(i), '.mat');
    load(file_txt);
    
%     plot(OneoverThetaPhiE_lin, f1, 'r*');
%     plot(OneoverThetaPhiE_lin, f2, 'r*');
%     plot(OneoverThetaPhiE_lin, f3, 'r*');
%     plot(OneoverThetaPhiE_lin, f4, 'r*');
%     plot(OneoverThetaPhiE_lin, f5, 'r*');

    if (i == 0)
        OneoverThetaPhiE_total = OneoverThetaPhiE_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
    else
        OneoverThetaPhiE_total = [OneoverThetaPhiE_total OneoverThetaPhiE_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
    end
end

% % plot(OneoverThetaPhiE_total, f1_total, 'r*');
% % plot(OneoverThetaPhiE_total, f2_total, 'g*');
% % plot(OneoverThetaPhiE_total, f3_total, 'b*');
% plot(OneoverThetaPhiE_total, f4_total, 'k*');
% % plot(OneoverThetaPhiE_total, f5_total, 'y*');

% Filter
N_filter = 10;
i_filter = 0;
[rows, cols] = size(OneoverThetaPhiE_total);
for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        OneoverThetaPhiE_filter(1, i_filter) = OneoverThetaPhiE_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
    end
end

% Combine
[rows, cols] = size(OneoverThetaPhiE_filter);

i_comb = 0;
OneoverThetaPhiE_comb = 0;
f_comb = 0;
for i = 1:1:cols
    if (isnan(f4_filter(1, i)) == 0)
        i_comb = i_comb + 1;

        OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
        f_comb(1, i_comb) = f4_filter(1, i);
    end
end

% f_comb(OneoverThetaPhiE_comb > 0.464) = NaN;
% plot(OneoverThetaPhiE_comb, f_comb, 'r-', 'LineWidth', 4);

% 
% plot(OneoverThetaPhiE_comb, f_comb, '-', 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);

%% Plot Analytic PING

thetaPhiE = 1./linspace(0.42, 0.52, 100);
epsilonEI = -0.5;       % From I to E

m_H = -log(exp(-2.*m_tau) - (1 - exp(-thetaPhiE)).*epsilonEI);
T = 2.*m_tau + thetaPhiE - m_H;

plot(1./thetaPhiE, 1./T, 'r-', 'LineWidth', LWSize)

%% Plot PINGING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v5\PINGING\PINGING_Types11BifDiagFullVaryPhiE', num2str(i), '.mat');
    load(file_txt);
    
%     plot(OneoverThetaPhiE_lin, f1, 'g*');
%     plot(OneoverThetaPhiE_lin, f2, 'r*');
%     plot(OneoverThetaPhiE_lin, f3, 'b*');
%     plot(OneoverThetaPhiE_lin, f4, 'r*');
%     plot(OneoverThetaPhiE_lin, f5, 'b*');
    
    if (i == 0)
        OneoverThetaPhiE_total = OneoverThetaPhiE_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
    else
        OneoverThetaPhiE_total = [OneoverThetaPhiE_total OneoverThetaPhiE_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
    end    
end

% % plot(OneoverThetaPhiE_total, f1_total, 'r*');
% plot(OneoverThetaPhiE_total, f2_total, 'g*');
% plot(OneoverThetaPhiE_total, f3_total, 'b*');
% plot(OneoverThetaPhiE_total, f4_total, '-', 'LineWidth', LWSize, 'Color', [101/255 163/255 92/255]);
plot(OneoverThetaPhiE_total, f4_total, '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255]);
% % plot(OneoverThetaPhiE_total, f5_total, 'y*');
% plot([0.4624 0.4638], [0.3511 0.352], '-', 'LineWidth', LWSize, 'Color', [101/255 163/255 92/255]); %% patch
plot([0.4624 0.4638], [0.3511 0.352], '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255]); %% patch

% Filter
N_filter = 1;
i_filter = 0;
[rows, cols] = size(OneoverThetaPhiE_total);
for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        OneoverThetaPhiE_filter(1, i_filter) = OneoverThetaPhiE_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
    end
end

% Combine
[rows, cols] = size(OneoverThetaPhiE_filter);

i_comb = 0;
OneoverThetaPhiE_comb = 0;
f_comb = 0;
for i = 1:1:cols
    if (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    end
    if (isnan(f2_filter(1, i)) == 0)
        i_comb = i_comb + 1;

        OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
        f_comb(1, i_comb) = f2_filter(1, i);
    end
end

% plot(OneoverThetaPhiE_filter, f2_filter, 'g-', 'LineWidth', 4);
% plot(OneoverThetaPhiE_filter, f3_filter, 'g-', 'LineWidth', 4);

plot(OneoverThetaPhiE_comb, f_comb, 'g-', 'LineWidth', LWSize);

%% Plot Analytic PINGING1
m_tau = 0.4;
epsilonEI = -0.5;       % From I to E
epsilonIE = 0.1;       % From E to I
epsilonII = -1.00;  % From I to I

thetaPhiE = 1./linspace(0.4194, 0.438, 10);
thetaPhiI = 1/0.495;
DeltaTheta = thetaPhiE - thetaPhiI;

H = -log(exp(-m_tau) - (1 - exp(-thetaPhiI)).*epsilonII);
E = H + DeltaTheta;

A = exp(-DeltaTheta).*(1 - exp(-thetaPhiE)).*epsilonEI;
B = exp(-E) - exp(-m_tau);
C = -(1 - exp(-thetaPhiI)).*epsilonIE;

DeltaPhi1 = log((-B + sqrt(B.*B - 4.*A.*C))./(2.*A));
DeltaPhi2 = log((-B - sqrt(B.*B - 4.*A.*C))./(2.*A));

PhiE1 = DeltaPhi1 + thetaPhiI;
PhiE2 = DeltaPhi2 + thetaPhiI;

T1 = m_tau + PhiE1 + log(exp(-(m_tau + DeltaPhi1 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI);
T2 = m_tau + PhiE2 + log(exp(-(m_tau + DeltaPhi2 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI);

figure(31);
subplot(2,1,1);hold on
plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*DeltaTheta, 'r');
% plot(1./thetaPhiE, 1./T1, 'r-*');
plot(1./thetaPhiE, DeltaPhi1, 'r-*');
plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*(DeltaTheta - m_tau), 'b');

subplot(2,1,2);hold on
plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*DeltaTheta, 'r');
% plot(1./thetaPhiE, 1./T2, 'r-*');
plot(1./thetaPhiE, DeltaPhi2, 'r-*');
plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*(DeltaTheta - m_tau), 'b');

%% Plot Analytic PINGING2
thetaPhiE = 1./linspace(0.4388, 0.4695, 10);
% thetaPhiE = 1./0.45;
thetaPhiI = 1/0.495;
DeltaTheta = thetaPhiE - thetaPhiI;

C = (1 - exp(-thetaPhiE)).*epsilonEI;
D = (1 - exp(-thetaPhiI)).*epsilonII;
E = (1 - exp(-thetaPhiI)).*epsilonIE;

A = C;
B = exp(-m_tau) - D - exp(-m_tau + DeltaTheta);
C = -E.*exp(DeltaTheta);

DeltaPhi1 = log((-B + sqrt(B.*B - 4.*A.*C))./(2.*A));
DeltaPhi2 = log((-B - sqrt(B.*B - 4.*A.*C))./(2.*A));

% m_tau + DeltaPhi1 - DeltaTheta

T1 = thetaPhiE + log(exp(-(m_tau + DeltaPhi1 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI) + m_tau + DeltaPhi1 - DeltaTheta;
T2 = thetaPhiE + log(exp(-(m_tau + DeltaPhi2 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI) + m_tau + DeltaPhi2 - DeltaTheta;

% figure(32);
% subplot(2,1,1);hold on
figure(1);
% plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*DeltaTheta, 'r');
% plot(1./thetaPhiE, 1./T1, 'r-o');
% plot(1./thetaPhiE, DeltaPhi1, 'r-*');
% plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*(DeltaTheta - m_tau), 'b');

% subplot(2,1,2);hold on
% plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*DeltaTheta, 'r');
figure(1);
% plot(1./thetaPhiE, 1./T2, 'r-*');
% plot(1./thetaPhiE, DeltaPhi2, 'r-*');
% plot(1./thetaPhiE, thetaPhiE./thetaPhiE.*(DeltaTheta - m_tau), 'b');

make_me_pretty(gcf, ...
    gca, 40, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',[0.42 0.44 0.46 0.48 0.50 0.52],'XTickLabel',{'';'';'';'';'';''});
set(gca,'YTick',[0.32 0.34 0.36 0.38],'YTickLabel',{'';'';'';''});

ylim([0.32 0.39]);
xlim([0.42 0.52]);

maximize_a_fig(gcf);
% % m_savefig('Types11BifDiagFullVaryPhiE_v2', 'eps');
% % m_savefig('Types11BifDiagFullVaryPhiE_v3', 'eps');
% m_savefig('Types11BifDiagFullVaryPhiE_v4', 'eps');

%% Calculate analytically DeltaPhi of PING
tau = 0.4;
ThetaI = 1/0.495;
ThetaE = 1/0.52;

epsiII = -1.0; % I->I
epsiEI = -0.5; % I->E

UP = exp(-tau) - (1 - exp(-ThetaI))*epsiII;
DOWN = exp(-2*tau) - (1 - exp(-ThetaE))*epsiEI;

DeltaPhi = log(UP/DOWN)
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