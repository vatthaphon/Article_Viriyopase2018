function plotTypes11BifDiagFullVaryPhiE
clc;
clear all;
close all;

root_txt = 'E:\paper2_Raoul\Sim_two_neurons_Raoul\Types11BifDiagFullVaryPhiE';

figure(1); hold on;
%% Plot ING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v1\ING\ING_Types11BifDiagFullVaryPhiE', num2str(i), '.mat');
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

plot(OneoverThetaPhiE_comb, f_comb, 'b-', 'LineWidth', 4);

%% Plot PING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v1\PING\PING_Types11BifDiagFullVaryPhiE', num2str(i), '.mat');
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

plot(OneoverThetaPhiE_comb, f_comb, 'r-', 'LineWidth', 4);

f_comb(OneoverThetaPhiE_comb < 0.4944) = NaN;

plot(OneoverThetaPhiE_comb, f_comb, '-', 'LineWidth', 4, 'Color', [101/255 163/255 92/255]);

%% Plot PINGING
for i = 0:1:11
    file_txt = strcat(root_txt, '\v1\PINGING\PINGING_Types11BifDiagFullVaryPhiE', num2str(i), '.mat');
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
    if (isnan(f4_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
%         
%         OneoverThetaPhiE_comb(1, i_comb) = OneoverThetaPhiE_filter(1, i);
%         f_comb(1, i_comb) = f4_filter(1, i);
%   
    end
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

plot(OneoverThetaPhiE_comb, f_comb, 'g-', 'LineWidth', 4);


% plot(OneoverThetaPhiE_total, f1_total, 'r*');
% plot(OneoverThetaPhiE_total, f2_total, 'g*');
% plot(OneoverThetaPhiE_total, f3_total, 'b*');
% plot(OneoverThetaPhiE_total, f4_total, 'k*');
% plot(OneoverThetaPhiE_total, f5_total, 'y*');

make_me_pretty(gcf, ...
    gca, 40, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

% set(gca,'XTick',[0.47 0.49 0.51 0.53 0.55],'XTickLabel',{'';'';'';'';''});
% set(gca,'YTick',[0.35 0.37 0.39 0.41],'YTickLabel',{'';'';'';''});

% ylim([0.34 0.42]);
xlim([0.4493 0.6]);

maximize_a_fig(gcf);
savefig('Types11BifDiagFullVaryPhiE', 'eps');
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