function plotTypes11QIFTau04BifDiagFullVaryPhiEPhiI
clc;
clear all;
close all;

root_txt = '/Volumes/Data/paper2_Raoul/Sim_two_neurons_Raoul/TypesQIF11BifDiagFullVaryPhiEPhiI';
LW = 20;
% LWSize = 30;
LWSize = 4 + 10;

scenarios_color = {'y*', 'm*', 'c*', 'r*', 'g*', 'b*'};
% scenarios_color = {'y', 'm', 'c', 'r', 'g', 'b'};
scenarios_color_id = [1 2 3 4 5 6];

RGB_code(1, :) = [1 1 0];
RGB_code(2, :) = [1 0 1];
RGB_code(3, :) = [0 1 1];
RGB_code(4, :) = [1 0 0];
RGB_code(5, :) = [0 1 0];
RGB_code(6, :) = [0 0 1];

PINGING_FaceAlpha = 0.5 + 0.0;
PING_FaceAlpha = 0.7 - 0.0;
ING_FaceAlpha = 0.7 - 0.0;

% ThetaPhiI_min = 2.2357 - 0.06;
% ThetaPhiI_max = 2.2357 + 0.06;
% 
% ThetaPhiE_min = 2.027 - 0.06;
% ThetaPhiE_max = 2.027 + 0.06;

% N_OneoverThetaPhiE_lin = 2;
% N_OneoverThetaPhiI_lin = 10;
% 
% OneoverThetaPhiI_max_show = 1./ThetaPhiI_min;
% OneoverThetaPhiI_min_show = 1./ThetaPhiI_max;
% 
% OneoverThetaPhiE_min_show = 1./ThetaPhiE_max;
% OneoverThetaPhiE_max_show = 1./ThetaPhiE_min;

version = 'v5';

N_OneoverThetaPhiI_lin_mid = 50;
% N_OneoverThetaPhiI_lin_mid = 5;

N_OneoverThetaPhiE_lin_mid = 10;
% N_OneoverThetaPhiE_lin_mid = 2;

file_txt = strcat(root_txt, '/ING/', version, '/ING_TypesQIF11Tau04BifDiagFullVaryPhiEVaryPhiI1.mat');
load(file_txt);
OneoverThetaPhiE_min = min(OneoverThetaPhiE_lin);

file_txt = strcat(root_txt, '/ING/', version, '/ING_TypesQIF11Tau04BifDiagFullVaryPhiEVaryPhiI10.mat');
load(file_txt);
OneoverThetaPhiE_max = max(OneoverThetaPhiE_lin);
OneoverThetaPhiI_min = min(OneoverThetaPhiI_lin);
OneoverThetaPhiI_max = max(OneoverThetaPhiI_lin);

N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
N_OneoverThetaPhiI_lin = size(firstIte_FP_DeltaPsi, 3);

OneoverThetaPhiE_mean = (OneoverThetaPhiE_min + OneoverThetaPhiE_max)/2
OneoverThetaPhiI_mean = (OneoverThetaPhiI_min + OneoverThetaPhiI_max)/2
N_FP = size(firstIte_FP_DeltaPsi, 4);


ING_FP_varyThetaPhiI = [];
ING_FP_varyThetaPhiI_i = 1;

% ING_FP_varyThetaPhiE
% 
% PING_FP_varyThetaPhiI
% PING_FP_varyThetaPhiE
% 
% PINGING_FP_varyThetaPhiI
% PINGING_FP_varyThetaPhiE

%% Plot ING
for i = 1:1:10
    file_txt = strcat(root_txt, '/ING/', version, '/ING_TypesQIF11Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
    N_OneoverThetaPhiI_lin = size(firstIte_FP_DeltaPsi, 3);
    N_FP = size(firstIte_FP_DeltaPsi, 4);
    
    %% First Iteration
    for ii = 1:1:5
        for j = 1:1:N_OneoverThetaPhiE_lin
            OneoverThetaPhiE = OneoverThetaPhiE_lin(1, j);
            for k = 1:1:N_OneoverThetaPhiI_lin
                OneoverThetaPhiI = OneoverThetaPhiI_lin(1, k);
                
                FP = squeeze(firstIte_FP_DeltaPsi(ii, j, k, :));
                [N_fixepoint, nCols] = size(FP);
                for l = 1:1:N_fixepoint
                    if (isnan(FP(l, 1)) == 1)
                    else                                                
                        switch(ii)
                            case 1
                                display('ING:1');
                            case 2
                                f = period_scenario2(gammaE, thetaVE, 1/OneoverThetaPhiE, FP(l, 1), tau, epsilonEI);
%                                 display('ING:2');
                            case 3
                                f = period_scenario3(gammaE, thetaVE, 1/OneoverThetaPhiE, FP(l, 1), tau, epsilonEI);
%                                 display('ING:3');                
                            case 4
                                display('ING:4');
                            case 5
                                display('ING:5');
                            otherwise
                                fprintf('Invalid scenarios\n' );
                        end
                        
                        FP(l, 1) = f;
                    end
                end
                
                sorted_FP = sort(FP);
                firstIte_FP_DeltaPsi(ii, j, k, :) = sorted_FP;                
            end
        end
    end

    %% First Iteration
    for i_FP = 1:1:N_FP        
        for ii = 1:1:5
% % %             if (i_FP == 1)
% %                 data = squeeze(firstIte_FP_DeltaPsi(ii, :, :, i_FP));
% %                 
% %                 figure(1); hold on;
% %                 hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
% %                 set(hSurface,'FaceColor', [102/255 102/255 255/255], 'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
% % %             end
        end        
    end

    %% Second Iteration
    for ii = 1:1:5
        for iii = 1:1:5
            for j = 1:1:N_OneoverThetaPhiE_lin
                OneoverThetaPhiE = OneoverThetaPhiE_lin(1, j);
                for k = 1:1:N_OneoverThetaPhiI_lin
                    OneoverThetaPhiI = OneoverThetaPhiI_lin(1, k);
                    
                    FP = squeeze(secondIte_FP_DeltaPsi (ii, iii, j, k, :));
                    [N_fixepoint, nCols] = size(FP);
                    for l = 1:1:N_fixepoint
                        if (isnan(FP(l, 1)) == 1)
                        else
                            if ((ii == 5) && (iii == 1))
                                f = period_scenario5then1(gammaE, gammaI, 1, 1, 1/OneoverThetaPhiE, 1/OneoverThetaPhiI, FP(l, 1), tau, epsilonEI, epsilonIE);
                                FP(l, 1) = f;                               
                            else
                                if  ((ii == iii)  || ((ii == 1) && (iii == 5)))
                                else
                                    msg = strcat('ING:', num2str(ii), ',', num2str(iii));
                                    display(msg);
                                end
                                
                                FP(l, 1) = nan;                                
                            end
                        end
                    end
                    
                    sorted_FP = sort(FP);
                    secondIte_FP_DeltaPsi (ii, iii, j, k, :) = sorted_FP;     
                    
                    if (sum(isnan(sorted_FP) == 0) > 1)
%                         'd'
                    end
                end
            end
        end
    end

    %% Second Iteration
%     for i_FP = 1:1:N_FP
    for i_FP = 1:1:1
        for ii = 5
            for iii = 1
                data = squeeze(secondIte_FP_DeltaPsi (ii, iii, :, :, i_FP));                                
                
%                 figure(1); hold on;
%                 hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
%                 set(hSurface,'FaceColor', [102/255 102/255 255/255], 'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
                
                
                if (i == 5)
                    figure(2); hold on;
                    plot(OneoverThetaPhiI_lin, data(N_OneoverThetaPhiE_lin_mid, :), 'b-', 'LineWidth', LWSize, 'Color', [102/255 102/255 255/255])                    
                end
                
                figure(3); hold on;
                plot(OneoverThetaPhiE_lin, data(:, N_OneoverThetaPhiI_lin_mid), 'b-', 'LineWidth', LWSize, 'Color', [102/255 102/255 255/255])                    
            end
        end
    end
end

%% Plot PING
for i = 1:1:10
    file_txt = strcat(root_txt, '/PINGING/', version, '/PINGING_TypesQIF11Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
    N_OneoverThetaPhiI_lin = size(firstIte_FP_DeltaPsi, 3);
    
    f_PING = NaN(N_OneoverThetaPhiE_lin, N_OneoverThetaPhiI_lin);
    
    for j = 1:1:N_OneoverThetaPhiE_lin
        OneoverThetaPhiE = OneoverThetaPhiE_lin(1, j);
        thetaPhiE = 1/OneoverThetaPhiE;
        I_E = (pi/thetaPhiE)*(pi/thetaPhiE);
        
        tmp = Hepsilon(1, 2.*tau, epsilonEI, I_E, 1, thetaPhiE, 3);
        
        T = 2.*tau + thetaPhiE - tmp;
        
        f = 1./T;
        
        for k = 1:1:N_OneoverThetaPhiI_lin
            OneoverThetaPhiI = OneoverThetaPhiI_lin(1, k);
            
            f_PING(j, k) = f;
            
        end
        
    end
    
%     figure(1); hold on;
%     hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, f_PING);
%     c_code = [RGB_code(scenarios_color_id(1, 4), 1) RGB_code(scenarios_color_id(1, 4), 2) RGB_code(scenarios_color_id(1, 4), 3)];
%     set(hSurface,'FaceColor', c_code, 'FaceAlpha', PING_FaceAlpha, 'EdgeColor', 'none');
    
    
    if (i == 5)
        figure(2); hold on;
        plot(OneoverThetaPhiI_lin, f_PING(N_OneoverThetaPhiE_lin_mid, :), 'r-', 'LineWidth', LWSize)                    
    end
                
    figure(3); hold on;
    plot(OneoverThetaPhiE_lin, f_PING(:, N_OneoverThetaPhiI_lin_mid), 'r-', 'LineWidth', LWSize)                    
end

deltaPsi_after_sce5_max_f = -1000;
deltaPsi_after_sce5_thetaPhiI = NaN;
deltaPsi_after_sce5_thetaPhiE = NaN;

%% Plot PINGING
for i = 1:1:10
    file_txt = strcat(root_txt, '/PINGING/', version, '/PINGING_TypesQIF11Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
    N_OneoverThetaPhiI_lin = size(firstIte_FP_DeltaPsi, 3);
    N_FP = size(firstIte_FP_DeltaPsi, 4);
    
    %% First Iteration
    for ii = 1:1:5
        for j = 1:1:N_OneoverThetaPhiE_lin
            OneoverThetaPhiE = OneoverThetaPhiE_lin(1, j);
            for k = 1:1:N_OneoverThetaPhiI_lin
                OneoverThetaPhiI = OneoverThetaPhiI_lin(1, k);
                
                FP = squeeze(firstIte_FP_DeltaPsi(ii, j, k, :));
                [N_fixepoint, nCols] = size(FP);
                for l = 1:1:N_fixepoint
                    if (isnan(FP(l, 1)) == 1)
                        FP(l, 1) = nan;
                    else                                                
                        switch(ii)
                            case 1
                                display('ING:1');
                            case 2
                                f = period_scenario2(gammaE, thetaVE, 1/OneoverThetaPhiE, FP(l, 1), tau, epsilonEI);
                            case 3
                                f = period_scenario3(gammaE, thetaVE, 1/OneoverThetaPhiE, FP(l, 1), tau, epsilonEI);
                            case 4
                                display('ING:4');
                            case 5
                                display('ING:5');
                            otherwise
                                fprintf('Invalid scenarios\n' );
                        end
                        
                        FP(l, 1) = f;
                    end
                end
                
                sorted_FP = sort(FP);
                firstIte_FP_DeltaPsi(ii, j, k, :) = sorted_FP;     
                
                if (sum(sum((isnan(FP) == 0)) > 1))
                    'sum(sum((isnan(FP) == 0)) > 1'
                end
            end
        end
    end

    %% First Iteration
    for i_FP = 1:1:N_FP
        for ii = 1:1:5
            data = squeeze(firstIte_FP_DeltaPsi(ii, :, :, i_FP));
            
            data_freq_low = data;
            data_freq_high = data;
            
            data_freq_low(data_freq_low > 0.42505) = NaN;
            data_freq_high(data_freq_high < 0.42505) = NaN;
            
            data = data_freq_low;
%             figure(1); hold on;
%             hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
%                         
%             c_code = [RGB_code(scenarios_color_id(1, 5), 1) RGB_code(scenarios_color_id(1, 5), 2) RGB_code(scenarios_color_id(1, 5), 3)];
%             set(hSurface,'FaceColor', c_code, 'FaceAlpha', 0.7, 'EdgeColor', 'none');

            
            if (i == 5)
                figure(2); hold on;                
                plot(OneoverThetaPhiI_lin, data(N_OneoverThetaPhiE_lin_mid, :), 'g-', 'LineWidth', LWSize)
                
                if (sum(isnan(data(N_OneoverThetaPhiE_lin_mid, :)) == 0) > 0)
                    [OneoverThetaPhiI_lin_ING_f, I] = min(data(N_OneoverThetaPhiE_lin_mid, :));
                    OneoverThetaPhiI_lin_ING_I = OneoverThetaPhiI_lin(1, I);                
                end
            end
            
            figure(3); hold on;
            plot(OneoverThetaPhiE_lin, data(:, N_OneoverThetaPhiI_lin_mid), 'g-', 'LineWidth', LWSize)
            
            if (i_FP == 2)
                if (i == 5)
                    if (sum(isnan(data(:, N_OneoverThetaPhiI_lin_mid)) == 0) > 0)
                        [OneoverThetaPhiE_lin_ING_f, I] = max(data(:, N_OneoverThetaPhiI_lin_mid));
                        OneoverThetaPhiE_lin_ING_I = OneoverThetaPhiE_lin(1, I);                                            
                    end 
                end
            end
        
            
            if (sum(sum((isnan(data_freq_high) == 0))) > 0)
                % (E, I)
                [nRows, nCols] = size(data_freq_high);
                for nRows_i = 1:1:nRows
                    for nCols_i = 1:1:nCols
                        if ((isnan(data_freq_high(nRows_i, nCols_i)) == 1) && (nCols_i > 1))
                            good_data = data_freq_high(nRows_i, nCols_i - 1);
                            for tmp_i = 0:1:2
                                if ((nCols_i + tmp_i) <= nCols)
                                    data_freq_high(nRows_i, nCols_i + tmp_i) = good_data;
                                end
                            end
                            
                            break;
                        end
                    end
                end
            end
            
            data = data_freq_high;
%             figure(1); hold on;
%             hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
%                         
%             c_code = [RGB_code(scenarios_color_id(1, 5), 1) RGB_code(scenarios_color_id(1, 5), 2) RGB_code(scenarios_color_id(1, 5), 3)];
%             set(hSurface,'FaceColor', c_code, 'FaceAlpha', 0.7, 'EdgeColor', 'none');

            
            if (i == 5)
                figure(2); hold on;                
                plot(OneoverThetaPhiI_lin, data(N_OneoverThetaPhiE_lin_mid, :), 'g-', 'LineWidth', LWSize)
                
                if (sum(isnan(data(N_OneoverThetaPhiE_lin_mid, :)) == 0) > 0)
                    [tmp, I] = min(data(N_OneoverThetaPhiE_lin_mid, :));
                    if (tmp < 0.42505)
                        [OneoverThetaPhiI_lin_ING_f, I] = min(data(N_OneoverThetaPhiE_lin_mid, :));
                        OneoverThetaPhiI_lin_ING_I = OneoverThetaPhiI_lin(1, I);                
                    end
                end
            end
            
            figure(3); hold on;
            plot(OneoverThetaPhiE_lin, data(:, N_OneoverThetaPhiI_lin_mid), 'g-', 'LineWidth', LWSize)
            
            if (i_FP == 1)
                if (i == 5)
                    if (sum(isnan(data(:, N_OneoverThetaPhiI_lin_mid)) == 0) > 0)
                        [OneoverThetaPhiE_lin_ING_f, I] = max(data(:, N_OneoverThetaPhiI_lin_mid));
                        OneoverThetaPhiE_lin_ING_I = OneoverThetaPhiE_lin(1, I);                                            
                    end 
                end
            end
            
            if (sum(isnan(data_freq_low(:, N_OneoverThetaPhiI_lin_mid)) == 0) > 0)
                [tmp, I] = min(data_freq_low(:, N_OneoverThetaPhiI_lin_mid));
                if (tmp < 0.42505)
                    [OneoverThetaPhiE_lin_PING_f_left, I] = min(data_freq_low(:, N_OneoverThetaPhiI_lin_mid));
                    OneoverThetaPhiE_lin_PING_I_left = OneoverThetaPhiE_lin(1, I);
                end
            end
        end
    end

    %% Second Iteration
    for ii = 1:1:5
        for iii = 1:1:5
            for j = 1:1:N_OneoverThetaPhiE_lin
                OneoverThetaPhiE = OneoverThetaPhiE_lin(1, j);
                for k = 1:1:N_OneoverThetaPhiI_lin
                    OneoverThetaPhiI = OneoverThetaPhiI_lin(1, k);
                    
                    FP = squeeze(secondIte_FP_DeltaPsi (ii, iii, j, k, :));
                    [N_fixepoint, nCols] = size(FP);
                    for l = 1:1:N_fixepoint
                        if (isnan(FP(l, 1)) == 1)
                        else
                            if ((ii == 5) && (iii == 1))
                                tmp_delta_psi = FP(l, 1);
                                f = period_scenario5then1(gammaE, gammaI, 1, 1, 1/OneoverThetaPhiE, 1/OneoverThetaPhiI, FP(l, 1), tau, epsilonEI, epsilonIE);
                                FP(l, 1) = f;                            
                                
                                % Determine when do they spike after the arrival of
                                % the E input.
                                deltaPsi_after_sce5 = cal_delta_psi_after_scenario5(gammaE, gammaI, 1, 1, 1/OneoverThetaPhiE, 1/OneoverThetaPhiI, tmp_delta_psi, tau, epsilonEI, epsilonIE);
                                duration_to_spike =  1/OneoverThetaPhiI + (deltaPsi_after_sce5 + (1/OneoverThetaPhiE -  1/OneoverThetaPhiI) - tau);
                                
                                tmp_T = 1/f;
                                
                                figure(6);hold on
                                plot(1, duration_to_spike/tmp_T, '*');                                
%                                 plot(1, duration_to_spike, '*');     
                                
                                if (duration_to_spike > deltaPsi_after_sce5_max_f)
                                    deltaPsi_after_sce5_max_f = duration_to_spike;
                                    deltaPsi_after_sce5_thetaPhiI = OneoverThetaPhiI;
                                    deltaPsi_after_sce5_thetaPhiE = OneoverThetaPhiE;
                                    
                                    deltaPsi_after_sce5_before_sce5 = tmp_delta_psi;
                                    deltaPsi_after_sce5_after_sce5 = deltaPsi_after_sce5;
                                    deltaPsi_after_sce5_after_sce1 = cal_delta_psi_after_scenario1(gammaE, gammaI, 1, 1, 1/OneoverThetaPhiE, 1/OneoverThetaPhiI, deltaPsi_after_sce5, tau, epsilonEI, epsilonIE, epsilonII);
                                    
                                    deltaPsi_after_sce5_duration_to_spike = duration_to_spike;
                                end
                            else
                                if  ((ii == iii)  || ((ii == 1) && (iii == 5)))
                                else
                                    msg = strcat('ING:', num2str(ii), ',', num2str(iii));
                                    display(msg);
                                end
                                
                                FP(l, 1) = nan;                                
                            end
                        end
                    end
                    
                    sorted_FP = sort(FP);
                    secondIte_FP_DeltaPsi (ii, iii, j, k, :) = sorted_FP;  
                    
                    if (sum(sum((isnan(FP) == 0)) > 1))
                        'sum(sum((isnan(FP) == 0)) > 1'
                    end
                end
            end
        end
    end

    %% Second Iteration
    for i_FP = 1:1:N_FP
        for ii = 5
            for iii = 1
                data = squeeze(secondIte_FP_DeltaPsi (ii, iii, :, :, i_FP));                                
                
%                 figure(1); hold on;
%                 hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
%                                 
% %                 if (i_FP == 2)
% %                     c_code = [RGB_code(scenarios_color_id(1, 5), 1) RGB_code(scenarios_color_id(1, 5), 2) RGB_code(scenarios_color_id(1, 5), 3)];
% %                     set(hSurface,'FaceColor', c_code, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
% %                 else
%                     set(hSurface,'FaceColor',[0/255 102/255 92/255],'FaceAlpha', 0.9, 'EdgeColor', 'none');
% %                 end
                
                if (i == 5)
                    figure(2); hold on;
                    plot(OneoverThetaPhiI_lin, data(N_OneoverThetaPhiE_lin_mid, :), '-*', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255])
                    
                    if (sum(isnan(data(N_OneoverThetaPhiE_lin_mid, :)) == 0) > 0)
                        [OneoverThetaPhiI_lin_PING_f, I] = max(data(N_OneoverThetaPhiE_lin_mid, :));                    
                        OneoverThetaPhiI_lin_PING_I = OneoverThetaPhiI_lin(1, I);                        
                    end
                end
                
                figure(3); hold on;
                plot(OneoverThetaPhiE_lin, data(:, N_OneoverThetaPhiI_lin_mid), '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255])
                
                if (i == 6)
                    if (i_FP == 1)
                        if (sum(isnan(data(:, N_OneoverThetaPhiI_lin_mid)) == 0) > 0)
                            [OneoverThetaPhiE_lin_PING_f_right, I] = min(data(:, N_OneoverThetaPhiI_lin_mid));
                            OneoverThetaPhiE_lin_PING_I_right = OneoverThetaPhiE_lin(1, I);                        
                        end
                    end
                end                
            end
        end
    end
end

deltaPsi_after_sce5_max_f
deltaPsi_after_sce5_thetaPhiI
deltaPsi_after_sce5_thetaPhiE
deltaPsi_after_sce5_before_sce5
deltaPsi_after_sce5_after_sce5
deltaPsi_after_sce5_after_sce1
tau
deltaPsi_after_sce5_duration_to_spike

%% Patches the curves
figure(2); hold on;
plot([OneoverThetaPhiI_lin_PING_I (OneoverThetaPhiI_lin_PING_I + OneoverThetaPhiI_lin_ING_I)/2], [OneoverThetaPhiI_lin_PING_f (OneoverThetaPhiI_lin_PING_f + OneoverThetaPhiI_lin_ING_f)/2], '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255]);
plot([(OneoverThetaPhiI_lin_PING_I + OneoverThetaPhiI_lin_ING_I)/2 OneoverThetaPhiI_lin_ING_I], [(OneoverThetaPhiI_lin_PING_f + OneoverThetaPhiI_lin_ING_f)/2 OneoverThetaPhiI_lin_ING_f], 'g-', 'LineWidth', LWSize)

figure(3); hold on;
% plot([OneoverThetaPhiE_lin_PING_I (OneoverThetaPhiE_lin_PING_I + OneoverThetaPhiE_lin_ING_I)/2], [OneoverThetaPhiE_lin_PING_f (OneoverThetaPhiE_lin_PING_f + OneoverThetaPhiE_lin_ING_f)/2], '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255]);
% plot([(OneoverThetaPhiE_lin_PING_I + OneoverThetaPhiE_lin_ING_I)/2 OneoverThetaPhiE_lin_ING_I], [(OneoverThetaPhiE_lin_PING_f + OneoverThetaPhiE_lin_ING_f)/2 OneoverThetaPhiE_lin_ING_f], 'g-', 'LineWidth', LWSize)

plot([OneoverThetaPhiE_lin_PING_I_left (OneoverThetaPhiE_lin_PING_I_left + OneoverThetaPhiE_lin_PING_I_right)/2], [OneoverThetaPhiE_lin_PING_f_left (OneoverThetaPhiE_lin_PING_f_left + OneoverThetaPhiE_lin_PING_f_right)/2], 'g-', 'LineWidth', LWSize);
plot([(OneoverThetaPhiE_lin_PING_I_left + OneoverThetaPhiE_lin_PING_I_right)/2 OneoverThetaPhiE_lin_PING_I_right], [(OneoverThetaPhiE_lin_PING_f_left + OneoverThetaPhiE_lin_PING_f_right)/2 OneoverThetaPhiE_lin_PING_f_right], '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255])


% 
% figure(1)
% axis tight
% view(-133, 20)
% set(gca,'YDir','reverse')
% zlim([0.4184 0.4229]);
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% set(gca,'XTick',[0.439 0.44 0.441],'XTickLabel',{'';'';''});
% set(gca,'YTick',[0.473 0.474 0.475 0.476],'YTickLabel',{'';'';'';''});
% set(gca,'ZTick',[0.419 0.420 0.421 0.422],'ZTickLabel',{'';'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiEVaryPhiI_v2_down', 'eps');

% zlim([0.4272 0.4355]);
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% set(gca,'XTick',[0.439 0.440 0.441],'XTickLabel',{'';'';''});
% set(gca,'YTick',[0.473 0.474 0.475 0.476],'YTickLabel',{'';'';'';''});
% set(gca,'ZTick',[0.428 0.430 0.432 0.434],'ZTickLabel',{'';'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiEVaryPhiI_v2_up', 'eps');

% figure(2); hold on;
% xlim([(0.4401 - 0.0018) (0.4401 + 0.0018)]);
% ylim([0.4183 0.42207]);
% set(gca,'XTick',[0.439 0.44 0.441],'XTickLabel',{'';'';''});
% set(gca,'YTick',[0.419 0.420 0.421 0.422],'YTickLabel',{'';'';'';''});
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% m_savefig('TypesQIF11BifDiagFullVaryPhiI_v2_down', 'eps');

% figure(3); hold on;
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% xlim([0.4744 - 0.0018 0.4744 + 0.0018]);
% ylim([0.4184 0.4214]);
% set(gca,'XTick',[0.473 0.474 0.475 0.476],'XTickLabel',{'';'';'';''});
% set(gca,'YTick',[0.419 0.420 0.421],'YTickLabel',{'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiE_v2_down', 'eps');
% 
% xlim([0.10029 0.10065]);
% ylim([0.1001 0.10061]);
% set(gca,'XTick',[0.1003 0.1004 0.1005 0.1006],'XTickLabel',{'';'';''});
% set(gca,'YTick',[0.1002 0.1004 0.1006],'YTickLabel',{'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiE_up', 'eps');
    
end

function f = period_scenario5then1(gammaE, gammaI, thetaVE, thetaVI, thetaPhiE, thetaPhiI, delta_psi, tau, eIE, eEI)

I_E = (pi*pi)/(thetaPhiE*thetaPhiE);
I_I = (pi*pi)/(thetaPhiI*thetaPhiI);

Hsine1 = Hepsilon(gammaI, thetaPhiI + tau - delta_psi, eEI, I_I, thetaVI, thetaPhiI, 3);

HLIF = Hepsilon(gammaE, 2*tau + thetaPhiI - Hsine1, eIE, I_E, thetaVE, thetaPhiE, 3);

Hsine2 = Hepsilon(gammaI, thetaPhiI + tau - delta_psi, eEI, I_I, thetaVI, thetaPhiI, 3);

T = 2*tau +  thetaPhiE + thetaPhiI - Hsine2 - HLIF;
f = 1/T;

end

function next_delta_psi = cal_delta_psi_after_scenario5(gammaE, gammaI, thetaVE, thetaVI, thetaPhiE, thetaPhiI, delta_psi, tau, eIE, eEI)

I_I = (pi*pi)/(thetaPhiI*thetaPhiI);

Hsine2 = Hepsilon(gammaI, thetaPhiI + tau - delta_psi, eEI, I_I, thetaVI, thetaPhiI, 3);

next_delta_psi = tau - Hsine2 - (thetaPhiE - thetaPhiI);

end

function next_delta_psi = cal_delta_psi_after_scenario1(gammaE, gammaI, thetaVE, thetaVI, thetaPhiE, thetaPhiI, delta_psi, tau, eIE, eEI, eII)

I_E = (pi*pi)/(thetaPhiE*thetaPhiE);
I_I = (pi*pi)/(thetaPhiI*thetaPhiI);

Hsine1 = Hepsilon(gammaE, thetaPhiE + delta_psi + tau, eIE, I_E, thetaVE, thetaPhiE, 3);

Hsine2 = Hepsilon(gammaI, tau, eII, I_I, thetaVI, thetaPhiI, 3);

next_delta_psi = Hsine1 - Hsine2 - (thetaPhiE - thetaPhiI);

end

function f = period_scenario2(gammaE, thetaVE, thetaPhiE, delta_psi, tau, eIE)

I_E = (pi*pi)/(thetaPhiE*thetaPhiE);

HLIF = Hepsilon(gammaE, tau + delta_psi, eIE, I_E, 1, thetaPhiE, 3);

T = tau +  delta_psi + thetaPhiE - HLIF;
f = 1/T;

end

function f = period_scenario3(gammaE, thetaVE, thetaPhiE, delta_psi, tau, eIE)

I_E = (pi*pi)/(thetaPhiE*thetaPhiE);

HLIF = Hepsilon(gammaE, tau + delta_psi, eIE, I_E, 1, thetaPhiE, 3);

T = tau +  delta_psi + thetaPhiE - HLIF;
f = 1/T;

end
