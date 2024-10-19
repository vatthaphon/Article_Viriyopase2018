function plotTypes12BifDiagFullVaryPhiI
clc;
clear all;
close all;

root_txt = 'C:\paper2_Raoul\Sim_two_neurons_Raoul\Types12BifDiagFullVaryPhiEPhiI';
LW = 20;

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

ThetaPhiI_min = 1/0.54;
ThetaPhiI_max = 1/0.46;

ThetaPhiE_min = 1/0.85;
ThetaPhiE_max = 1/0.63;

N_OneoverThetaPhiE_lin = 10;
N_OneoverThetaPhiI_lin = 100;

OneoverThetaPhiI_max_show = 0.5182;
OneoverThetaPhiI_min_show = 0.5 - (OneoverThetaPhiI_max_show - 0.5);

OneoverThetaPhiE_min_show = 0.7058;
OneoverThetaPhiE_max_show = 0.74 + (0.74 - OneoverThetaPhiE_min_show);

%% Plot ING
N_OneoverThetaPhiE_lin_total = 0;
for i = 1:1:10
    file_txt = strcat(root_txt, '\ING\v1\ING_Types12Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
    N_OneoverThetaPhiI_lin = size(firstIte_FP_DeltaPsi, 3);
    
    N_OneoverThetaPhiE_lin_total = N_OneoverThetaPhiE_lin + N_OneoverThetaPhiE_lin_total;
end

% ING_freq = NaN(N_OneoverThetaPhiI_lin, N_OneoverThetaPhiE_lin_total); % [m, n] = size(ING_freq) where m = size(N_OneoverThetaPhiI_lin) and n = size(N_OneoverThetaPhiE_lin).
ING_freq  = [];
OneoverThetaPhiE_lin_total = [];
N_OneoverThetaPhiE_lin_id = 0;
for i = 1:1:10
    file_txt = strcat(root_txt, '\ING\v1\ING_Types12Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
    N_OneoverThetaPhiI_lin = size(firstIte_FP_DeltaPsi, 3);
    %% If there is more than one fixed point for each scenarios.
    for j = 1:1:5
        for jj = 1:1:N_OneoverThetaPhiE_lin
            for jjj = 1:1:N_OneoverThetaPhiI_lin
                if (sum(isnan(firstIte_FP_DeltaPsi(j, jj, jjj, :))) <= 3)
                    'More than one fixed point'
                end
            end
        end
        
        for k = 1:1:5
            for jj = 1:1:N_OneoverThetaPhiE_lin
                for jjj = 1:1:N_OneoverThetaPhiI_lin
                    if (sum(isnan(secondIte_FP_DeltaPsi(j, k, jj, jjj, :))) <= 3)
                        'More than one fixed point'
                    end
                end
            end
        end
    end
    
    for j = 1:1:5
        % One iteration
        tmp_deltaPsi = [];
        for jj = 1:1:N_OneoverThetaPhiE_lin
            for jjj = 1:1:N_OneoverThetaPhiI_lin
                tmp_deltaPsi(jjj, jj) = firstIte_FP_DeltaPsi(j, jj, jjj, 1);
                
                if (isnan(tmp_deltaPsi(jjj, jj)) == 0)
                    switch(j)
                        case 1
                        case 2
                            tmp_f(jjj, jj) = period_scenario2(gammaE, thetaVE, 1/OneoverThetaPhiE_lin(jj), tmp_deltaPsi(jjj, jj), tau, epsilonEI);
                        case 3
                            tmp_f(jjj, jj) = period_scenario3(gammaE, thetaVE, 1/OneoverThetaPhiE_lin(jj), tmp_deltaPsi(jjj, jj), tau, epsilonEI);
                        case 4
                        case 5
                        otherwise
                            fprintf('Invalid scenarios\n' );
                    end
                else
                    tmp_f(jjj, jj) = NaN;
                end
            end
        end
        
        %         figure(1);hold on;
        %         view(-44, 66)
        %         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin,tmp_deltaPsi);
        %         c_code = [RGB_code(scenarios_color_id(1, j), 1) RGB_code(scenarios_color_id(1, j), 2) RGB_code(scenarios_color_id(1, j), 3)];
        %         set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
        %         xlabel('E');ylabel('I')
        %         axis tight
        
        %         figure(2);hold on;
        %         view(-44, 66)
        %         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_f);
        %         c_code = [RGB_code(scenarios_color_id(1, j), 1) RGB_code(scenarios_color_id(1, j), 2) RGB_code(scenarios_color_id(1, j), 3)];
        %         set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
        %         xlabel('E');ylabel('I')
        %         axis tight
        
        %         figure(3);hold on;
        %         view(-44, 66)
        %         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_f);
        %         c_code = [RGB_code(scenarios_color_id(1, 6), 1) RGB_code(scenarios_color_id(1, 6), 2) RGB_code(scenarios_color_id(1, 6), 3)];
        %         set(hSurface,'FaceColor', c_code, 'FaceAlpha', ING_FaceAlpha, 'EdgeColor', c_code);
        %         xlabel('E');ylabel('I')
        %         axis tight
        
        ING_freq = [ING_freq tmp_f];
        OneoverThetaPhiE_lin_total = [OneoverThetaPhiE_lin_total OneoverThetaPhiE_lin];
        
        
        %         PINGING_OneoverThetaPhiI_lin = [PINGING_OneoverThetaPhiI_lin OneoverThetaPhiI_lin];
        %         PINGING_freq = [PINGING_freq tmp_f];
        %
        % Two iterations
        for k = 1:1:5
            tmp_deltaPsi = [];
            for jj = 1:1:N_OneoverThetaPhiE_lin
                for jjj = 1:1:N_OneoverThetaPhiI_lin
                    tmp_deltaPsi(jjj, jj) = secondIte_FP_DeltaPsi(j, k, jj, jjj, 1);
                    
                    if (isnan(tmp_deltaPsi(jjj, jj)) == 0)
                        if ((j == 5) && (k == 1))
                            tmp_f(jjj, jj) = period_scenario5then1(gammaE, gammaI, thetaVE, thetaVI, 1/OneoverThetaPhiE_lin(jj), 1/OneoverThetaPhiI_lin(jjj), tmp_deltaPsi(jjj, jj), tau, epsilonEI, epsilonIE);
                        end
                    else
                        tmp_f(jjj, jj) = NaN;
                    end
                end
            end
            
            %         figure(1);hold on;
            %         view(-44, 66)
            %         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_deltaPsi);
            %         c_code = [RGB_code(scenarios_color_id(1, j), 1) RGB_code(scenarios_color_id(1, j), 2) RGB_code(scenarios_color_id(1, j), 3)];
            %         set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
            %         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_deltaPsi);
            %
            %         c_code = [RGB_code(scenarios_color_id(1, k), 1) RGB_code(scenarios_color_id(1, k), 2) RGB_code(scenarios_color_id(1, k), 3)];
            %         set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
            %         xlabel('E');ylabel('I')
            %         axis tight
            
            %         figure(2);hold on;
            %         view(-44, 66)
            %         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_f);
            %         c_code = [RGB_code(scenarios_color_id(1, j), 1) RGB_code(scenarios_color_id(1, j), 2) RGB_code(scenarios_color_id(1, j), 3)];
            %         set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
            %         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_f);
            %
            %         c_code = [RGB_code(scenarios_color_id(1, k), 1) RGB_code(scenarios_color_id(1, k), 2) RGB_code(scenarios_color_id(1, k), 3)];
            %         set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
            %         xlabel('E');ylabel('I')
            %         axis tight
            
            %             figure(3);hold on;
            %             view(-44, 66)
            %             hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_f);
            %             c_code = [RGB_code(scenarios_color_id(1, 6), 1) RGB_code(scenarios_color_id(1, 6), 2) RGB_code(scenarios_color_id(1, 6), 3)];
            %             set(hSurface,'FaceColor', c_code, 'FaceAlpha', ING_FaceAlpha, 'EdgeColor', c_code);
            %             xlabel('E');ylabel('I')
            %             axis tight
            
            ING_freq = [ING_freq tmp_f];
            OneoverThetaPhiE_lin_total = [OneoverThetaPhiE_lin_total OneoverThetaPhiE_lin];
            
            %
            % %             figure(1);hold on;
            % %             plot(OneoverThetaPhiI_lin, tmp_deltaPsi, char(scenarios_color(1, j)), 'LineWidth', 14)
            % %             plot(OneoverThetaPhiI_lin, tmp_deltaPsi, char(scenarios_color(1, k)), 'LineWidth', 3)
            % %
            % %             figure(2);hold on;
            % %             plot(OneoverThetaPhiI_lin, tmp_f, char(scenarios_color(1, j)), 'LineWidth', 14)
            % %             plot(OneoverThetaPhiI_lin, tmp_f, char(scenarios_color(1, k)), 'LineWidth', 3)
            %
            %             PINGING_OneoverThetaPhiI_lin = [PINGING_OneoverThetaPhiI_lin OneoverThetaPhiI_lin];
            %             PINGING_freq = [PINGING_freq tmp_f];
            %
            %             PINGING_PINGOneoverThetaPhiI_lin = [PINGING_PINGOneoverThetaPhiI_lin OneoverThetaPhiI_lin];
            %             PINGING_PINGfreq = [PINGING_PINGfreq tmp_f];
        end
    end
end

[SortedOneoverThetaPhiE_lin_total, I_SortedOneoverThetaPhiE_lin_total] = sort(OneoverThetaPhiE_lin_total);
SortedING_freq = ING_freq(:, I_SortedOneoverThetaPhiE_lin_total);

N_tmp_I = size(ING_freq, 1);

ING_freq_show = NaN(N_tmp_I, N_OneoverThetaPhiE_lin_total);
for i = 1:1:N_tmp_I
    N_tmp_E_total_id = 0;
    for j = 1:1:N_OneoverThetaPhiE_lin_total
        for k = 1:1:30
            N_tmp_E_total_id = N_tmp_E_total_id + 1;
            if (isnan(SortedING_freq(i, N_tmp_E_total_id)) == 0)
                ING_freq_show(i, j) = SortedING_freq(i, N_tmp_E_total_id);
            end
        end
    end
end

OneoverThetaPhiE_lin = linspace(1/ThetaPhiE_max, 1/ThetaPhiE_min, N_OneoverThetaPhiE_lin_total);

figure(4);hold on;
view(-44, 66)

goodE = (OneoverThetaPhiE_min_show <= OneoverThetaPhiE_lin) & (OneoverThetaPhiE_lin <= OneoverThetaPhiE_max_show);
goodI = (OneoverThetaPhiI_min_show <= OneoverThetaPhiI_lin) & (OneoverThetaPhiI_lin <= OneoverThetaPhiI_max_show);

hSurface = surf(OneoverThetaPhiE_lin(goodE), OneoverThetaPhiI_lin(goodI), ING_freq_show(goodI, goodE));
c_code = [RGB_code(scenarios_color_id(1, 6), 1) RGB_code(scenarios_color_id(1, 6), 2) RGB_code(scenarios_color_id(1, 6), 3)];
% set(hSurface,'FaceColor', c_code, 'FaceAlpha', ING_FaceAlpha, 'EdgeColor', c_code);
set(hSurface,'FaceColor', [102/255 102/255 255/255], 'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
xlabel('E');ylabel('I')
axis tight

%% Plot Analytic PING
N_OneoverThetaPhiE_lin_total = 0;
for i = 1:1:10
    file_txt = strcat(root_txt, '\ING\v1\ING_Types12Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
    N_OneoverThetaPhiI_lin = size(firstIte_FP_DeltaPsi, 3);
    
    N_OneoverThetaPhiE_lin_total = N_OneoverThetaPhiE_lin + N_OneoverThetaPhiE_lin_total;
end

OneoverThetaPhiE_lin = linspace(1/ThetaPhiE_max, 1/ThetaPhiE_min, N_OneoverThetaPhiE_lin_total);


PING_freq  = NaN(N_OneoverThetaPhiI_lin, N_OneoverThetaPhiE_lin_total);
for i = 1:1:N_OneoverThetaPhiE_lin_total
    thetaPhiE = 1/OneoverThetaPhiE_lin(i);

    m_H = -log(exp(-2.*tau) - (1 - exp(-thetaPhiE)).*epsilonEI);
    T_E = 2.*tau + thetaPhiE - m_H;    
    
    for j = 1:1:N_OneoverThetaPhiI_lin
        PING_freq(j, i) = 1/T_E;
    end
end

figure(4);hold on;
view(-44, 66)

goodE = (OneoverThetaPhiE_min_show <= OneoverThetaPhiE_lin) & (OneoverThetaPhiE_lin <= OneoverThetaPhiE_max_show);
goodI = (OneoverThetaPhiI_min_show <= OneoverThetaPhiI_lin) & (OneoverThetaPhiI_lin <= OneoverThetaPhiI_max_show);

hSurface = surf(OneoverThetaPhiE_lin(goodE), OneoverThetaPhiI_lin(goodI), PING_freq(goodI, goodE));
c_code = [RGB_code(scenarios_color_id(1, 4), 1) RGB_code(scenarios_color_id(1, 4), 2) RGB_code(scenarios_color_id(1, 4), 3)];
% set(hSurface,'FaceColor', c_code, 'FaceAlpha', PING_FaceAlpha, 'EdgeColor', c_code);
set(hSurface,'FaceColor', c_code, 'FaceAlpha', PING_FaceAlpha, 'EdgeColor', 'none');
xlabel('E');ylabel('I')
axis tight

%% Plot PINGING
PINGING_freq_up_surface  = [];
OneoverThetaPhiE_lin_total_up_surface = [];

PINGING_freq_normal_surface_ING  = [];
OneoverThetaPhiE_lin_total_normal_surface_ING = [];

PINGING_freq_normal_surface_PING  = [];
OneoverThetaPhiE_lin_total_normal_surface_PING = [];
for i = 1:1:10
    file_txt = strcat(root_txt, '\PINGING\v1\PINGING_Types12Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
    N_OneoverThetaPhiE_lin = size(firstIte_FP_DeltaPsi, 2);
    N_OneoverThetaPhiI_lin = size(firstIte_FP_DeltaPsi, 3);
    %% If there is more than one fixed point for each scenarios.
    for j = 1:1:5
        for jj = 1:1:N_OneoverThetaPhiE_lin
            for jjj = 1:1:N_OneoverThetaPhiI_lin
                if (sum(isnan(firstIte_FP_DeltaPsi(j, jj, jjj, :))) <= 3)
                    'More than one fixed point'
                end
            end
        end
        
        for k = 1:1:5
            for jj = 1:1:N_OneoverThetaPhiE_lin
                for jjj = 1:1:N_OneoverThetaPhiI_lin
                    if (sum(isnan(secondIte_FP_DeltaPsi(j, k, jj, jjj, :))) <= 3)
                        'More than one fixed point'
                    end
                end
            end
        end
    end
    
    for j = 1:1:5
        % One iteration
        tmp_deltaPsi = [];
        for jj = 1:1:N_OneoverThetaPhiE_lin
            for jjj = 1:1:N_OneoverThetaPhiI_lin
                tmp_deltaPsi(jjj, jj) = firstIte_FP_DeltaPsi(j, jj, jjj, 1);
                
                if (isnan(tmp_deltaPsi(jjj, jj)) == 0)
                    switch(j)
                        case 1
                        case 2
                            tmp_f(jjj, jj) = period_scenario2(gammaE, thetaVE, 1/OneoverThetaPhiE_lin(jj), tmp_deltaPsi(jjj, jj), tau, epsilonEI);
                        case 3
                            tmp_f(jjj, jj) = period_scenario3(gammaE, thetaVE, 1/OneoverThetaPhiE_lin(jj), tmp_deltaPsi(jjj, jj), tau, epsilonEI);
                        case 4
                        case 5
                        otherwise
                            fprintf('Invalid scenarios\n' );
                    end
                else
                    tmp_f(jjj, jj) = NaN;
                end
            end
        end
        
%         figure(1);hold on;
%         view(-44, 66)
%         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin,tmp_deltaPsi);
%         c_code = [RGB_code(scenarios_color_id(1, j), 1) RGB_code(scenarios_color_id(1, j), 2) RGB_code(scenarios_color_id(1, j), 3)];
%         set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
%         xlabel('E');ylabel('I')
%         axis tight
%         
%         figure(2);hold on;
%         view(-44, 66)
%         hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_f);
%         c_code = [RGB_code(scenarios_color_id(1, j), 1) RGB_code(scenarios_color_id(1, j), 2) RGB_code(scenarios_color_id(1, j), 3)];
%         set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
%         xlabel('E');ylabel('I')
%         axis tight
        
        if (j == 2)
            PINGING_freq_up_surface = [PINGING_freq_up_surface tmp_f];
            OneoverThetaPhiE_lin_total_up_surface = [OneoverThetaPhiE_lin_total_up_surface OneoverThetaPhiE_lin];
        end
        
        if (j == 3)
            PINGING_freq_normal_surface_ING = [PINGING_freq_normal_surface_ING tmp_f];
            OneoverThetaPhiE_lin_total_normal_surface_ING = [OneoverThetaPhiE_lin_total_normal_surface_ING OneoverThetaPhiE_lin];
        end
        
        % Two iterations
        for k = 1:1:5
            tmp_deltaPsi = [];
            for jj = 1:1:N_OneoverThetaPhiE_lin
                for jjj = 1:1:N_OneoverThetaPhiI_lin
                    tmp_deltaPsi(jjj, jj) = secondIte_FP_DeltaPsi(j, k, jj, jjj, 1);
                    
                    if (isnan(tmp_deltaPsi(jjj, jj)) == 0)
                        if ((j == 5) && (k == 1))
                            tmp_f(jjj, jj) = period_scenario5then1(gammaE, gammaI, thetaVE, thetaVI, 1/OneoverThetaPhiE_lin(jj), 1/OneoverThetaPhiI_lin(jjj), tmp_deltaPsi(jjj, jj), tau, epsilonEI, epsilonIE);
                        end
                    else
                        tmp_f(jjj, jj) = NaN;
                    end
                end
            end
            
%             figure(1);hold on;
%             view(-44, 66)
%             hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_deltaPsi);
%             c_code = [RGB_code(scenarios_color_id(1, j), 1) RGB_code(scenarios_color_id(1, j), 2) RGB_code(scenarios_color_id(1, j), 3)];
%             set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
%             hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_deltaPsi);
%             
%             c_code = [RGB_code(scenarios_color_id(1, k), 1) RGB_code(scenarios_color_id(1, k), 2) RGB_code(scenarios_color_id(1, k), 3)];
%             set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
%             xlabel('E');ylabel('I')
%             axis tight
%             
%             figure(2);hold on;
%             view(-44, 66)
%             hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_f);
%             c_code = [RGB_code(scenarios_color_id(1, j), 1) RGB_code(scenarios_color_id(1, j), 2) RGB_code(scenarios_color_id(1, j), 3)];
%             set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
%             hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, tmp_f);
%             
%             c_code = [RGB_code(scenarios_color_id(1, k), 1) RGB_code(scenarios_color_id(1, k), 2) RGB_code(scenarios_color_id(1, k), 3)];
%             set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
%             xlabel('E');ylabel('I')
%             axis tight
            
            if ((j == 5) && (k == 1))
                PINGING_freq_normal_surface_PING = [PINGING_freq_normal_surface_PING tmp_f];
                OneoverThetaPhiE_lin_total_normal_surface_PING = [OneoverThetaPhiE_lin_total_normal_surface_PING OneoverThetaPhiE_lin];
            end
        end
    end
end

[UniqueOneoverThetaPhiE_lin_total_up_surface, ia, ic] = unique(OneoverThetaPhiE_lin_total_up_surface);
UniquePINGING_freq_up_surface = PINGING_freq_up_surface(:, ia);

figure(4);hold on;
view(-44, 66)

goodE = (OneoverThetaPhiE_min_show <= UniqueOneoverThetaPhiE_lin_total_up_surface) & (UniqueOneoverThetaPhiE_lin_total_up_surface <= OneoverThetaPhiE_max_show);
goodI = (OneoverThetaPhiI_min_show <= OneoverThetaPhiI_lin) & (OneoverThetaPhiI_lin <= OneoverThetaPhiI_max_show);

% Light green
hSurface = surf(UniqueOneoverThetaPhiE_lin_total_up_surface(goodE), OneoverThetaPhiI_lin(goodI), UniquePINGING_freq_up_surface(goodI, goodE));
c_code = [RGB_code(scenarios_color_id(1, 5), 1) RGB_code(scenarios_color_id(1, 5), 2) RGB_code(scenarios_color_id(1, 5), 3)];
% set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
set(hSurface,'FaceColor', [167/255 212/255 140/255], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
xlabel('E');ylabel('I')
axis tight

[UniqueOneoverThetaPhiE_lin_total_normal_surface_ING, ia, ic] = unique(OneoverThetaPhiE_lin_total_normal_surface_ING);
UniquePINGING_freq_normal_surface_ING = PINGING_freq_normal_surface_ING(:, ia);

figure(4);hold on;
view(-44, 66)

goodE = (OneoverThetaPhiE_min_show <= UniqueOneoverThetaPhiE_lin_total_normal_surface_ING) & (UniqueOneoverThetaPhiE_lin_total_normal_surface_ING <= OneoverThetaPhiE_max_show);
goodI = (OneoverThetaPhiI_min_show <= OneoverThetaPhiI_lin) & (OneoverThetaPhiI_lin <= OneoverThetaPhiI_max_show);

% Light green
hSurface = surf(UniqueOneoverThetaPhiE_lin_total_normal_surface_ING(goodE), OneoverThetaPhiI_lin(goodI), UniquePINGING_freq_normal_surface_ING(goodI, goodE) - 0.0004);
c_code = [RGB_code(scenarios_color_id(1, 5), 1) RGB_code(scenarios_color_id(1, 5), 2) RGB_code(scenarios_color_id(1, 5), 3)];
% set(hSurface,'FaceColor', c_code, 'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', c_code);
set(hSurface,'FaceColor', [167/255 212/255 140/255], 'FaceAlpha', 0.7, 'EdgeColor', 'none');
xlabel('E');ylabel('I')
axis tight

[UniqueOneoverThetaPhiE_lin_total_normal_surface_PING, ia, ic] = unique(OneoverThetaPhiE_lin_total_normal_surface_PING);
UniquePINGING_freq_normal_surface_PING = PINGING_freq_normal_surface_PING(:, ia);

figure(4);hold on;
view(-44, 66)

goodE = (OneoverThetaPhiE_min_show <= UniqueOneoverThetaPhiE_lin_total_normal_surface_PING) & (UniqueOneoverThetaPhiE_lin_total_normal_surface_PING <= OneoverThetaPhiE_max_show);
goodI = (OneoverThetaPhiI_min_show <= OneoverThetaPhiI_lin) & (OneoverThetaPhiI_lin <= OneoverThetaPhiI_max_show);

% Dark green
hSurface = surf(UniqueOneoverThetaPhiE_lin_total_normal_surface_PING(goodE), OneoverThetaPhiI_lin(goodI), UniquePINGING_freq_normal_surface_PING(goodI, goodE) - 0.00085);
% - 0.00085 to get rid of the glitch between the red and dark green surfaces.
% set(hSurface,'FaceColor',[101/255 163/255 92/255],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [101/255 163/255 92/255]);
set(hSurface,'FaceColor',[20/255 156/255 73/255],'FaceAlpha', 0.9, 'EdgeColor', 'none');

xlabel('E');ylabel('I')
axis tight

figure(4);hold on;
view(-35, 8)
xlabel('');ylabel('')
make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',[0.71 0.73 0.75 0.77],'XTickLabel',{'';'';'';''});
set(gca,'YTick',[0.49 0.50 0.51],'YTickLabel',{'';'';''});
set(gca,'ZTick',[0.6 0.62 0.64 0.66 0.68],'ZTickLabel',{'';'';'';'';''});

maximize_a_fig(gcf);
m_savefig('Types12Tai04BifDiagFullVaryPhiEVaryPhiI_v4_printing', 'eps');

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