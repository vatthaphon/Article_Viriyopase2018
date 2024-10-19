function plotTypes11BifDiagFullVaryPhiI
clc;
clear all;
close all;

root_txt = 'C:\paper2_Raoul\Sim_two_neurons_Raoul\Types11BifDiagFullVaryPhiI';

figure(1); hold on;
%% Plot ING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v1\ING\ING_Types11BifDiagFullVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
% %     plot(OneoverThetaPhiI_lin, f1, 'b-', 'LineWidth', 4);
% %     plot(OneoverThetaPhiI_lin, f2, 'r*', 'LineWidth', 4);
%     plot(OneoverThetaPhiI_lin, f3, 'b-', 'LineWidth', 4);
%     plot(OneoverThetaPhiI_lin, f4, 'r-', 'LineWidth', 4);
% %     plot(OneoverThetaPhiI_lin, f5, 'b-', 'LineWidth', 4);
    
    if (i == 0)
        OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
    else
        OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
    end    
end

% plot(OneoverThetaPhiI_total, f1_total, 'r*');
% plot(OneoverThetaPhiI_total, f2_total, 'g*');
% plot(OneoverThetaPhiI_total, f3_total, 'b*');
% plot(OneoverThetaPhiI_total, f4_total, 'k*');
% plot(OneoverThetaPhiI_total, f5_total, 'y*');

% Filter
N_filter = 10;
i_filter = 0;
[rows, cols] = size(OneoverThetaPhiI_total);
for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
    end
end

% Combine
[rows, cols] = size(OneoverThetaPhiI_filter);

i_comb = 0;
OneoverThetaPhiI_comb = 0;
f_comb = 0;
for i = 1:1:cols
    if (isnan(f2_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f2_filter(1, i);
    elseif (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    elseif (isnan(f4_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f4_filter(1, i);
    end
end

plot(OneoverThetaPhiI_comb, f_comb, 'b-', 'LineWidth', 4);

% % %% Plot PING
% % for i = 1:1:10
% %     file_txt = strcat(root_txt, '\PING\PING_Types11BifDiagFullVaryPhiI', num2str(i), '.mat');
% %     load(file_txt);
% %     
% %     plot(OneoverThetaPhiI_lin, f1, 'r*');
% %     plot(OneoverThetaPhiI_lin, f2, 'r*');
% %     plot(OneoverThetaPhiI_lin, f3, 'r*');
% %     plot(OneoverThetaPhiI_lin, f4, 'r*');
% %     plot(OneoverThetaPhiI_lin, f5, 'r*');
% % end
% 
%% Plot PINGING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v1\PINGING\PINGING_Types11BifDiagFullVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
% %     plot(OneoverThetaPhiI_lin, f1, 'g*');
%     plot(OneoverThetaPhiI_lin, f2, 'g*');   % Top lines
%     plot(OneoverThetaPhiI_lin, f3, 'g*');   % Top lines
%     
%     plot(OneoverThetaPhiI_lin, f4, 'g*');   % Buttom lines
% %     plot(OneoverThetaPhiI_lin, f5, 'g*');

    if (i == 1)
        OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
    else
        OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
    end    

end

% plot(OneoverThetaPhiI_total, f1_total, 'r*');
% plot(OneoverThetaPhiI_total, f2_total, 'g*');
% plot(OneoverThetaPhiI_total, f3_total, 'b*');
% plot(OneoverThetaPhiI_total, f4_total, 'k*');
% plot(OneoverThetaPhiI_total, f5_total, 'y*');


% Filter
N_filter = 1;
i_filter = 0;
[rows, cols] = size(OneoverThetaPhiI_total);
for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
    end
end

% Combine
[rows, cols] = size(OneoverThetaPhiI_filter);

i_comb = 0;
i_comb1 = 0;
OneoverThetaPhiI_comb = 0;
OneoverThetaPhiI_comb1 = 0;
f_comb = 0;
f_comb1 = 0;
for i = 1:1:cols
    if (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    end
    
    if (isnan(f2_filter(1, i)) == 0)
        i_comb = i_comb + 1;

        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f2_filter(1, i);
    end
    
    if (isnan(f4_filter(1, i)) == 0)
        i_comb1 = i_comb1 + 1;

        OneoverThetaPhiI_comb1(1, i_comb1) = OneoverThetaPhiI_filter(1, i);
        f_comb1(1, i_comb1) = f4_filter(1, i);
    end
    
end

plot(OneoverThetaPhiI_comb, f_comb, 'g-', 'LineWidth', 4);
plot([OneoverThetaPhiI_comb1 0.5493], [f_comb1 f_comb1(1, end)], 'r-', 'LineWidth', 4);
plot([OneoverThetaPhiI_comb1 0.4956], [f_comb1 f_comb1(1, end)], 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);
% plot(OneoverThetaPhiI_comb1, f_comb1, 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);
% 
make_me_pretty(gcf, ...
    gca, 40, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

% set(gca,'XTick',[0.46 0.48 0.50 0.52 0.54],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0.35 0.37 0.39 0.41],'YTickLabel',{'';'';'';''});
% % 
maximize_a_fig(gcf);
ylim([0.34 0.4105]);
% xlim([0.45 0.5493]);
xlim([0.46 0.54]);
% savefig('Types11BifDiagFullVaryPhiI', 'eps');
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