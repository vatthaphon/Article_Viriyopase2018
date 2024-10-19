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

version = 'v4';

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
            data = squeeze(firstIte_FP_DeltaPsi(ii, :, :, i_FP));
            
%             figure(1); hold on;
%             hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
%             
%             c_code = [0 0 1];
%             set(hSurface,'FaceColor', c_code, 'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');            
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
                
                figure(1); hold on;
                hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
                set(hSurface,'FaceColor', [102/255 102/255 255/255], 'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
                
                
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
    
    figure(1); hold on;
    hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, f_PING);
    c_code = [RGB_code(scenarios_color_id(1, 4), 1) RGB_code(scenarios_color_id(1, 4), 2) RGB_code(scenarios_color_id(1, 4), 3)];
    set(hSurface,'FaceColor', c_code, 'FaceAlpha', PING_FaceAlpha, 'EdgeColor', 'none');
    
    
    if (i == 5)
        figure(2); hold on;
        plot(OneoverThetaPhiI_lin, f_PING(N_OneoverThetaPhiE_lin_mid, :), 'r-', 'LineWidth', LWSize)                    
    end
                
    figure(3); hold on;
    plot(OneoverThetaPhiE_lin, f_PING(:, N_OneoverThetaPhiI_lin_mid), 'r-', 'LineWidth', LWSize)                    
end

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
            end
        end
    end

    %% First Iteration
    for i_FP = 1:1:N_FP
        for ii = 1:1:5
            data = squeeze(firstIte_FP_DeltaPsi(ii, :, :, i_FP));
            
            figure(1); hold on;
            hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
                        
            c_code = [RGB_code(scenarios_color_id(1, 5), 1) RGB_code(scenarios_color_id(1, 5), 2) RGB_code(scenarios_color_id(1, 5), 3)];
            set(hSurface,'FaceColor', c_code, 'FaceAlpha', 0.7, 'EdgeColor', 'none');

            
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
            
            if (i_FP == 1)
                if (i == 5)
                    if (sum(isnan(data(:, N_OneoverThetaPhiI_lin_mid)) == 0) > 0)
                        [OneoverThetaPhiE_lin_ING_f, I] = max(data(:, N_OneoverThetaPhiI_lin_mid));
                        OneoverThetaPhiE_lin_ING_I = OneoverThetaPhiE_lin(1, I);                                            
                    end 
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
                end
            end
        end
    end

    %% Second Iteration
    for i_FP = 1:1:N_FP
        for ii = 5
            for iii = 1
                data = squeeze(secondIte_FP_DeltaPsi (ii, iii, :, :, i_FP));                                
                
                figure(1); hold on;
                hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, data);
                
                
                if (i_FP == 2)
                    c_code = [RGB_code(scenarios_color_id(1, 5), 1) RGB_code(scenarios_color_id(1, 5), 2) RGB_code(scenarios_color_id(1, 5), 3)];
                    set(hSurface,'FaceColor', c_code, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
                else
                    set(hSurface,'FaceColor',[0/255 102/255 92/255],'FaceAlpha', 0.9, 'EdgeColor', 'none');
                end
                
                if (i == 5)
                    figure(2); hold on;
                    plot(OneoverThetaPhiI_lin, data(N_OneoverThetaPhiE_lin_mid, :), '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255])
                    
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
                            [OneoverThetaPhiE_lin_PING_f, I] = min(data(:, N_OneoverThetaPhiI_lin_mid));
                            OneoverThetaPhiE_lin_PING_I = OneoverThetaPhiE_lin(1, I);                        
                        end
                    end
                end
            end
        end
    end
end

% %% Patches the curves
% figure(2); hold on;
% plot([OneoverThetaPhiI_lin_PING_I (OneoverThetaPhiI_lin_PING_I + OneoverThetaPhiI_lin_ING_I)/2], [OneoverThetaPhiI_lin_PING_f (OneoverThetaPhiI_lin_PING_f + OneoverThetaPhiI_lin_ING_f)/2], '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255]);
% plot([(OneoverThetaPhiI_lin_PING_I + OneoverThetaPhiI_lin_ING_I)/2 OneoverThetaPhiI_lin_ING_I], [(OneoverThetaPhiI_lin_PING_f + OneoverThetaPhiI_lin_ING_f)/2 OneoverThetaPhiI_lin_ING_f], 'g-', 'LineWidth', LWSize)
% 
% figure(3); hold on;
% plot([OneoverThetaPhiE_lin_PING_I (OneoverThetaPhiE_lin_PING_I + OneoverThetaPhiE_lin_ING_I)/2], [OneoverThetaPhiE_lin_PING_f (OneoverThetaPhiE_lin_PING_f + OneoverThetaPhiE_lin_ING_f)/2], '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255]);
% plot([(OneoverThetaPhiE_lin_PING_I + OneoverThetaPhiE_lin_ING_I)/2 OneoverThetaPhiE_lin_ING_I], [(OneoverThetaPhiE_lin_PING_f + OneoverThetaPhiE_lin_ING_f)/2 OneoverThetaPhiE_lin_ING_f], 'g-', 'LineWidth', LWSize)
% 
% % plot([], [], 'g-', 'LineWidth', LWSize)
% % plot([], [], '-', 'LineWidth', LWSize, 'Color', [0/255 102/255 92/255])
% 
% figure(1)
% axis tight
% view(-133, 20)
% set(gca,'YDir','reverse')
% zlim([0.09912 0.09958]);
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% % set(gca,'XTick',[0.1006 0.1008 0.1010],'XTickLabel',{'';'';''});
% % set(gca,'YTick',[0.1003 0.1005 0.1007],'YTickLabel',{'';'';''});
% 
% set(gca,'XTick',[0.1007 0.1009],'XTickLabel',{'';'';''});
% set(gca,'YTick',[0.1004 0.1006],'YTickLabel',{'';'';''});
% 
% set(gca,'ZTick',[0.0992 0.0993 0.0994 0.0995],'ZTickLabel',{'';'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiEVaryPhiI_down', 'eps');
% 
% zlim([0.1000 0.1007])
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% % set(gca,'XTick',[0.1006 0.1008 0.1010],'XTickLabel',{'';'';''});
% % set(gca,'YTick',[0.1003 0.1005 0.1007],'YTickLabel',{'';'';''});
% 
% set(gca,'XTick',[0.1007 0.1009],'XTickLabel',{'';'';''});
% set(gca,'YTick',[0.1004 0.1006],'YTickLabel',{'';'';''});
% 
% set(gca,'ZTick',[0.1001 0.1003 0.1005 0.1007],'ZTickLabel',{'';'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiEVaryPhiI_up', 'eps');
% 
% figure(2); hold on;
% set(gca, 'XTickLabel', num2str(get(gca, 'XTick'), '%.7f'))
% set(gca, 'YTickLabel', num2str(get(gca, 'YTick'), '%.7f'))
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% xlim([0.10065 0.101001]);
% ylim([0.09916 0.09953]);
% 
% set(gca,'XTick',[0.1007 0.1009],'XTickLabel',{'';''});
% set(gca,'YTick',[0.0992 0.0993 0.0994 0.0995],'YTickLabel',{'';'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiI_down', 'eps');
% 
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% xlim([0.10065 0.101001]);
% ylim([0.10033 0.100415]);
% 
% set(gca,'XTick',[0.1006 0.1008 0.1010],'XTickLabel',{'';'';''});
% set(gca,'YTick',[0.10035 0.100375 0.10040],'YTickLabel',{'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiI_up', 'eps');
% 
% figure(3); hold on;
% make_me_pretty(gcf, ...
%     gca, 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12, ...
%     [], 12)
% xlim([0.10029 0.10065]);
% ylim([0.09914 0.0995]);
% set(gca,'XTick',[0.1004 0.1006],'XTickLabel',{'';''});
% set(gca,'YTick',[0.0992 0.0993 0.0994],'YTickLabel',{'';'';''});
% m_savefig('TypesQIF11BifDiagFullVaryPhiE_down', 'eps');
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
