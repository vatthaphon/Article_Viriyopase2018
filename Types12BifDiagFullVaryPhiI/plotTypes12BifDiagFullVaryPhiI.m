function plotTypes12BifDiagFullVaryPhiI
clc;
clear all;
close all;

root_txt = 'E:\paper2_Raoul\Sim_two_neurons_Raoul\Types12BifDiagFullVaryPhiI';

figure(1); hold on;
%% Plot ING
for i = 1:1:9
    file_txt = strcat(root_txt, '\ING\v0\ING_Types12BifDiagFullVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
% %     plot(OneoverThetaPhiI_lin, f1, 'b*');
% %     plot(OneoverThetaPhiI_lin, f2, 'b*');
%     plot(OneoverThetaPhiI_lin, f3, 'b*');
% %     plot(OneoverThetaPhiI_lin, f4, 'b*');
% %     plot(OneoverThetaPhiI_lin, f5, 'b*');
%     plot(OneoverThetaPhiI_lin, f51, 'b*');
% %     plot(OneoverThetaPhiI_lin, f15, 'b*');

    if (i == 1)
        OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
        f51_total = f51;
        f15_total = f15;
    else
        OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
        f51_total = [f51_total f51];
        f15_total = [f15_total f15];
    end    
end

% Filter
N_filter = 10;
i_filter = 0;
[rows, cols] = size(OneoverThetaPhiI_total);
OneoverThetaPhiI_filter = 0;
f1_filter = 0;
f2_filter = 0;
f3_filter = 0;
f4_filter = 0;
f5_filter = 0;
f51_filter = 0;
f15_filter = 0;

for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
        f51_filter(1, i_filter) = f51_total(1, i);
        f15_filter(1, i_filter) = f15_total(1, i);
    end
end

% Combine
[rows, cols] = size(OneoverThetaPhiI_filter);

i_comb = 0;
OneoverThetaPhiI_comb = 0;
f_comb = 0;
for i = 1:1:cols
    if (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    elseif (isnan(f51_filter(1, i)) == 0)
        i_comb = i_comb + 1;

        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f51_filter(1, i);
    end
end

plot(OneoverThetaPhiI_comb, f_comb, 'b-', 'LineWidth', 4);

%% Plot PING
for i = 1:1:9
    file_txt = strcat(root_txt, '\PING\v0\PING_Types12BifDiagFullVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
% %     plot(OneoverThetaPhiI_lin, f1, 'r*');
% %     plot(OneoverThetaPhiI_lin, f2, 'r*');
%     plot(OneoverThetaPhiI_lin, f3, 'r*');
% %     plot(OneoverThetaPhiI_lin, f4, 'r*');
% %     plot(OneoverThetaPhiI_lin, f5, 'r*');
%     plot(OneoverThetaPhiI_lin, f51, 'r*');
% %     plot(OneoverThetaPhiI_lin, f15, 'r*');
    
    if (i == 1)
        OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
        f51_total = f51;
        f15_total = f15;
    else
        OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
        f51_total = [f51_total f51];
        f15_total = [f15_total f15];
    end      
end

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
%     elseif (isnan(f51_filter(1, i)) == 0)
%         i_comb = i_comb + 1;
% 
%         OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
%         f_comb(1, i_comb) = f51_filter(1, i);
%     end
% end
% 
% plot(OneoverThetaPhiI_comb, f_comb, 'r-', 'LineWidth', 4);
plot([0.501 0.5051], [0.6169 0.6169], 'r-', 'LineWidth', 4);
plot([0.5051 0.509], [0.6169 0.6175], 'r-', 'LineWidth', 4);


%% Plot PINGING
for i = 1:1:9
    file_txt = strcat(root_txt, '\PINGING\v0\PINGING_Types12BifDiagFullVaryPhiI', num2str(i), '.mat');
    load(file_txt);
    
% %     plot(OneoverThetaPhiI_lin, f1, 'g*');
% %     plot(OneoverThetaPhiI_lin, f2, 'g*');
%     plot(OneoverThetaPhiI_lin, f3, 'g*');
% %     plot(OneoverThetaPhiI_lin, f4, 'g*');
% %     plot(OneoverThetaPhiI_lin, f5, 'g*');
%     plot(OneoverThetaPhiI_lin, f51, 'g*');
% %     plot(OneoverThetaPhiI_lin, f15, 'g*');  

    if (i == 1)
        OneoverThetaPhiI_total = OneoverThetaPhiI_lin;
        f1_total = f1;
        f2_total = f2;
        f3_total = f3;
        f4_total = f4;
        f5_total = f5;
        f51_total = f51;
        f15_total = f15;
    else
        OneoverThetaPhiI_total = [OneoverThetaPhiI_total OneoverThetaPhiI_lin];
        f1_total = [f1_total f1];
        f2_total = [f2_total f2];
        f3_total = [f3_total f3];
        f4_total = [f4_total f4];
        f5_total = [f5_total f5];
        f51_total = [f51_total f51];
        f15_total = [f15_total f15];
    end      
end

% Filter
N_filter = 10;
i_filter = 0;
[rows, cols] = size(OneoverThetaPhiI_total);
OneoverThetaPhiI_filter = 0;
f1_filter = 0;
f2_filter = 0;
f3_filter = 0;
f4_filter = 0;
f5_filter = 0;
f51_filter = 0;
f15_filter = 0;

for i = 1:1:cols
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        i_filter = i_filter + 1;
        OneoverThetaPhiI_filter(1, i_filter) = OneoverThetaPhiI_total(1, i);
        f1_filter(1, i_filter) = f1_total(1, i);
        f2_filter(1, i_filter) = f2_total(1, i);
        f3_filter(1, i_filter) = f3_total(1, i);
        f4_filter(1, i_filter) = f4_total(1, i);
        f5_filter(1, i_filter) = f5_total(1, i);
        f51_filter(1, i_filter) = f51_total(1, i);
        f15_filter(1, i_filter) = f15_total(1, i);
    end
end

% Combine
[rows, cols] = size(OneoverThetaPhiI_filter);

i_comb = 0;
OneoverThetaPhiI_comb = 0;
f_comb = 0;
for i = 1:1:cols
    if (isnan(f3_filter(1, i)) == 0)
        i_comb = i_comb + 1;
        
        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f3_filter(1, i);
    elseif (isnan(f51_filter(1, i)) == 0)
        i_comb = i_comb + 1;

        OneoverThetaPhiI_comb(1, i_comb) = OneoverThetaPhiI_filter(1, i);
        f_comb(1, i_comb) = f51_filter(1, i);
    end
end

plot(OneoverThetaPhiI_comb, f_comb, 'g-', 'LineWidth', 4);

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',[0.501 0.503 0.505 0.507 0.509],'XTickLabel',{'';'';'';'';''});
set(gca,'YTick',[0.610 0.614 0.618 0.622],'YTickLabel',{'';'';'';''});

maximize_a_fig(gcf);
ylim([0.61 0.624]);
xlim([0.501 0.509]);
% savefig('Types12BifDiagFullVaryPhiI', 'eps');


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