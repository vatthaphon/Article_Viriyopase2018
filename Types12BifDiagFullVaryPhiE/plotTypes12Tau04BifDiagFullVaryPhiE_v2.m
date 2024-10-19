function plotTypes12BifDiagFullVaryPhiI
clc;
clear all;
close all;

root_txt = '/Volumes/Data/paper2_Raoul/Sim_two_neurons_Raoul/Types12BifDiagFullVaryPhiE';
LW = 25;

scenarios_color = {'y*', 'm*', 'c*', 'r*', 'g*'};
% scenarios_color = {'y', 'm', 'c', 'r', 'g'};

OneoverThetaPhiI_max_show = 0.5182;
OneoverThetaPhiI_min_show = 0.5 - (OneoverThetaPhiI_max_show - 0.5);

OneoverThetaPhiE_min_show = 0.7058;
OneoverThetaPhiE_max_show = 0.74 + (0.74 - OneoverThetaPhiE_min_show);

%% Plot ING
Analytic_ING_freq = [];
Analytic_ING_OneoverThetaPhiE_lin = [];

ING_freq = [];
ING_OneoverThetaPhiE_lin = [];
for i = 1:1:10
    file_txt = strcat(root_txt, '/ING/v3/ING_Types12Tau04BifDiagFullVaryPhiE', num2str(i), '.mat');
    load(file_txt);
    
    %% Analytic ING
    T = tau + thetaPhiI - thetaPhiI./pi.*atan(tan(pi./thetaPhiI.*tau).*exp(-2.*pi.*epsilonII./thetaPhiI));
    
    for OneoverThetaPhiE = OneoverThetaPhiE_lin
        Analytic_ING_OneoverThetaPhiE_lin = [Analytic_ING_OneoverThetaPhiE_lin OneoverThetaPhiE];
        Analytic_ING_freq = [Analytic_ING_freq 1./T];        
    end
    
    %% Sim ING
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
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
        
%         figure(1);hold on;
%         plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, j)), 'LineWidth', 7)
%         figure(2);hold on;
%         plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, j)), 'LineWidth', 7)
        
        ING_OneoverThetaPhiE_lin = [ING_OneoverThetaPhiE_lin OneoverThetaPhiE_lin];
        ING_freq = [ING_freq tmp_f];
                        
        % Two iterations
        for k = 1:1:5
            tmp_deltaPsi = [];
            for jj = 1:1:N_OneoverThetaPhiE_lin
                tmp_deltaPsi(jj) = secondIte_FP_DeltaPsi(j, k, jj, 1);
                
                if (isnan(tmp_deltaPsi(jj)) == 0)  
                    G2_id = strcat(num2str(j), num2str(k)); 
                    switch(G2_id)
                        case '15'
                        case '51'
                            
                            tmp_f(jj) = period_scenario5then1(gammaE, gammaI, thetaVE, thetaVI, 1/OneoverThetaPhiE_lin(jj), thetaPhiI, tmp_deltaPsi(jj), tau, epsilonEI, epsilonIE);                            
                        case strcat(num2str(j), num2str(j))
                        otherwise
                            
                            fprintf('G2 Invalid scenarios\n' );
                    end
                else
                    tmp_f(jj) = NaN;
                end
            end
        
%             figure(1);hold on;
%             plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, j)), 'LineWidth', 14)
%             plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, k)), 'LineWidth', 3)            
%             
%             figure(2);hold on;
%             plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, j)), 'LineWidth', 14)
%             plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, k)), 'LineWidth', 3)            
            
            ING_OneoverThetaPhiE_lin = [ING_OneoverThetaPhiE_lin OneoverThetaPhiE_lin];
            ING_freq = [ING_freq tmp_f];            
        end
    end    
end

figure(3);hold on;
plot(ING_OneoverThetaPhiE_lin, ING_freq./ING_freq.*max(max(ING_freq)), 'b-', 'LineWidth', LW, 'Color', [102/255 102/255 255/255])
% plot(Analytic_ING_OneoverThetaPhiE_lin(isnan(ING_OneoverThetaPhiE_lin) == 0) , Analytic_ING_freq(isnan(ING_OneoverThetaPhiE_lin) == 0), 'ro', 'LineWidth', 0.1)


%% Plot Analytic PING
Analytic_PING_freq = [];
Analytic_PING_OneoverThetaPhiE_lin = [];
Analytic_PING_freq_ext = [];
for i = 1:1:10    
    file_txt = strcat(root_txt, '/PINGING/v4/PINGING_Types12Tau04BifDiagFullVaryPhiE', num2str(i), '.mat');
    load(file_txt);
       
    T_I = tau + thetaPhiI - Hepsilon(gammaI, tau, epsilonII, NaN, thetaVI, thetaPhiI, 2);
    
    thetaPhiE = 1./OneoverThetaPhiE_lin;
    N_thetaPhiE = size(thetaPhiE, 2);
    
    
    for j = 1:1:N_thetaPhiE
        tmp_thetaPhiE = thetaPhiE(1, j);            
        
        m_H = -log(exp(-2.*tau) - (1 - exp(-tmp_thetaPhiE)).*epsilonEI);
        T_E = 2.*tau + tmp_thetaPhiE - m_H;
        
        Analytic_PING_OneoverThetaPhiE_lin = [Analytic_PING_OneoverThetaPhiE_lin 1/tmp_thetaPhiE];
%% This piece of code is not neccessary when we assume that in pure PING, 1/Theta_I is so low that the I neuron always spike by the infinite kick        
%% from the E neuron.
%         if (T_E <= T_I)
%             Analytic_PING_freq = [Analytic_PING_freq 1/T_E];                        
%             Analytic_PING_freq_ext = [Analytic_PING_freq_ext NaN];                        
%         else
%             Analytic_PING_freq = [Analytic_PING_freq NaN];                        
%             Analytic_PING_freq_ext = [Analytic_PING_freq_ext 1/T_E];                        
%         end
%         
%         Analytic_PING_freq = [Analytic_PING_freq 1/T_E];
%         Analytic_PING_freq_ext = [Analytic_PING_freq_ext NaN];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Analytic_PING_freq = [Analytic_PING_freq 1/T_E];
    end
end

figure(3);hold on;
plot(Analytic_PING_OneoverThetaPhiE_lin, Analytic_PING_freq, 'r-', 'LineWidth', LW)
% plot(Analytic_PING_OneoverThetaPhiE_lin, Analytic_PING_freq_ext, 'r--', 'LineWidth', 8)


%% Plot PINGING
PINGING_freq = [];
PINGING_OneoverThetaPhiE_lin = [];
PINGING_PINGfreq = [];
PINGING_PINGOneoverThetaPhiE_lin = [];
for i = 1:1:10
    file_txt = strcat(root_txt, '/PINGING/v4/PINGING_Types12Tau04BifDiagFullVaryPhiE', num2str(i), '.mat');
    load(file_txt);
    
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
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
        
%         figure(1);hold on;
%         plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, j)), 'LineWidth', 7)
%         figure(2);hold on;
%         plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, j)), 'LineWidth', 7)
        
        PINGING_OneoverThetaPhiE_lin = [PINGING_OneoverThetaPhiE_lin OneoverThetaPhiE_lin];
        PINGING_freq = [PINGING_freq tmp_f];        
        
        % Two iterations
        for k = 1:1:5
            tmp_deltaPsi = [];
            for jj = 1:1:N_OneoverThetaPhiE_lin
                tmp_deltaPsi(jj) = secondIte_FP_DeltaPsi(j, k, jj, 1);
                
                if (isnan(tmp_deltaPsi(jj)) == 0)
                    if ((j == 5) && (k == 1))
                        tmp_f(jj) = period_scenario5then1(gammaE, gammaI, thetaVE, thetaVI, 1/OneoverThetaPhiE_lin(jj), thetaPhiI, tmp_deltaPsi(jj), tau, epsilonEI, epsilonIE);
                        
                        % Determine when do they spike after the arrival of
                        % the E input.
                        deltaPsi_after_sce5 = cal_delta_psi_after_scenario5(gammaE, gammaI, thetaVE, thetaVI, 1/OneoverThetaPhiE_lin(jj), thetaPhiI, tmp_deltaPsi(jj), tau, epsilonEI, epsilonIE);
                        duration_to_spike = thetaPhiI + (deltaPsi_after_sce5 + (1/OneoverThetaPhiE_lin(jj) - thetaPhiI) - tau);
                        
                        tmp_T = 1/tmp_f(jj);
                        
                        figure(4);hold on
                        plot(1, duration_to_spike/tmp_T, '*');                        
%                         plot(1, duration_to_spike, '*');                        
%                         plot(1, tmp_f(jj), '*');                        
                    end
                else
                    tmp_f(jj) = NaN;
                end
            end
        
%             figure(1);hold on;
%             plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, j)), 'LineWidth', 14)
%             plot(OneoverThetaPhiE_lin, tmp_deltaPsi, char(scenarios_color(1, k)), 'LineWidth', 3)            
%             
%             figure(2);hold on;
%             plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, j)), 'LineWidth', 14)
%             plot(OneoverThetaPhiE_lin, tmp_f, char(scenarios_color(1, k)), 'LineWidth', 3)  
            
            PINGING_OneoverThetaPhiE_lin = [PINGING_OneoverThetaPhiE_lin OneoverThetaPhiE_lin];
            PINGING_freq = [PINGING_freq tmp_f];                        

            PINGING_PINGOneoverThetaPhiE_lin = [PINGING_PINGOneoverThetaPhiE_lin OneoverThetaPhiE_lin];
            PINGING_PINGfreq = [PINGING_PINGfreq tmp_f];                                    
        end
    end
end

figure(3);hold on;
plot(PINGING_OneoverThetaPhiE_lin, PINGING_freq, 'g-', 'LineWidth', LW)
% plot(PINGING_PINGOneoverThetaPhiE_lin, PINGING_PINGfreq, 'g-', 'LineWidth', LW, 'Color', [101/255 163/255 92/255])
plot(PINGING_PINGOneoverThetaPhiE_lin, PINGING_PINGfreq, 'g-', 'LineWidth', LW, 'Color', [0/255 102/255 92/255])

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',[0.71 0.73 0.75 0.77],'XTickLabel',{'';'';'';''});
set(gca,'YTick',[0.6 0.64 0.68],'YTickLabel',{'';'';''});

% maximize_a_fig(gcf);
ylim([0.58 0.70]);
xlim([OneoverThetaPhiE_min_show OneoverThetaPhiE_max_show]);
% m_savefig('Types12Tai04BifDiagFullVaryPhiE_v5', 'eps');
% m_savefig('Types12Tai04BifDiagFullVaryPhiE_v5', 'png');

% m_savefig('Types12Tai04BifDiagFullVaryPhiE_v6', 'eps');

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

function next_delta_psi = cal_delta_psi_after_scenario5(gammaE, gammaI, thetaVE, thetaVI, thetaPhiE, thetaPhiI, delta_psi, tau, eIE, eEI)

I_I = (gammaI*thetaVI)./(1 - exp(-gammaI*thetaPhiI));   

Hsine2 = Hepsilon(gammaI, thetaPhiI + tau - delta_psi, eEI, I_I, thetaVI, thetaPhiI, 2);

next_delta_psi = tau - Hsine2 - (thetaPhiE - thetaPhiI);

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