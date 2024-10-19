function plotTypes12Tau04BifDiagFullVaryPhiEPhiI_v2
clc;
clear all;
close all;

root_txt = 'C:\paper2_Raoul\Sim_two_neurons_Raoul\Types12BifDiagFullVaryPhiEPhiI';

OneoverThetaPhiE_min = 0.71;
OneoverThetaPhiE_max = 0.776;
OneoverThetaPhiI_min = 0.48;
OneoverThetaPhiI_max = 0.52;

PINGING_FaceAlpha = 0.5 + 0.0;
PING_FaceAlpha = 0.7 - 0.0;
ING_FaceAlpha = 0.7 - 0.0;

figure(1); hold on;
%% Plot ING
for i = 1:1:20
    file_txt = strcat(root_txt, '\ING\v0\ING_Types12Tau04BifDiagFullVaryPhiEPhiI', num2str(i), '.mat');
    
    load(file_txt);
    
    OneoverThetaPhiE_total = OneoverThetaPhiE_lin;
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

figure(1);hold on;xlabel('E');ylabel('I');
% % surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f1_total);
% surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f2_total);
% surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f3_total);
% % surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f4_total);
% % surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f5_total);
% surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f15_total);
% surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f51_total);

% Filter
N_filter = 1;
i_filterE = 0;
[rowsE, colsE] = size(OneoverThetaPhiE_total);
[rowsI, colsI] = size(OneoverThetaPhiI_total);
for i = 1:1:colsE
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        if ((OneoverThetaPhiE_min <= OneoverThetaPhiE_total(1, i)) && (OneoverThetaPhiE_total(1, i) <= OneoverThetaPhiE_max))
            i_filterE = i_filterE + 1;
            OneoverThetaPhiE_filter(1, i_filterE) = OneoverThetaPhiE_total(1, i);
            
            i_filterI = 0;
            for j = 1:1:colsI
                if (cmp(rem(j, N_filter), 0, 1e-6) == 0)
                    if ((OneoverThetaPhiI_min <= OneoverThetaPhiI_total(1, j)) && (OneoverThetaPhiI_total(1, j) <= OneoverThetaPhiI_max))
                        i_filterI = i_filterI + 1;
                        OneoverThetaPhiI_filter(1, i_filterI) = OneoverThetaPhiI_total(1, j);
                        
                        f1_filter(i_filterE, i_filterI) = f1_total(i, j);
                        f2_filter(i_filterE, i_filterI) = f2_total(i, j);
                        f3_filter(i_filterE, i_filterI) = f3_total(i, j);
                        f4_filter(i_filterE, i_filterI) = f4_total(i, j);
                        f5_filter(i_filterE, i_filterI) = f5_total(i, j);
                        f51_filter(i_filterE, i_filterI) = f51_total(i, j);
                        f15_filter(i_filterE, i_filterI) = f15_total(i, j);
                        
                    end
                end
            end
            
        end
    end
end

% hSurface = surf(OneoverThetaPhiE_filter, OneoverThetaPhiI_filter, f2_filter');
% set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);
% hSurface = surf(OneoverThetaPhiE_filter, OneoverThetaPhiI_filter, f3_filter');
% set(hSurface,'FaceColor',[1 1 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);
% hSurface = surf(OneoverThetaPhiE_filter, OneoverThetaPhiI_filter, f15_filter');
% set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);

% Combine
[rowsE, colsE] = size(OneoverThetaPhiE_filter);
[rowsI, colsI] = size(OneoverThetaPhiI_filter);
f_comb = NaN(colsE, colsI);

for i = 1:1:colsE
    for j = 1:1:colsI
        if (isnan(f2_filter(i, j)) == 0)
            f_comb(i, j) = f2_filter(i, j);
        elseif (isnan(f3_filter(i, j)) == 0)
            f_comb(i, j) = f3_filter(i, j);
        elseif (isnan(f51_filter(i, j)) == 0)
            f_comb(i, j) = f15_filter(i, j);
        else
            f_comb(i, j) = f15_filter(i, j);
        end
    end
end

tmp_ING = f_comb';
hSurface = surf(OneoverThetaPhiE_filter, OneoverThetaPhiI_filter, f_comb');
set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', [0 0 1]);

%% Plot Analytic PING
m_tau = 0.4;

[rowE, colE] = size(OneoverThetaPhiE_filter);
[rowI, colI] = size(OneoverThetaPhiI_filter);

f_comb = NaN(colE, colI);
for i = 1:1:colE
    for j = 1:1:colI
        thetaPhiE = 1./OneoverThetaPhiE_filter(1, i);
        thetaPhiI = 1./OneoverThetaPhiI_filter(1, j);
        
        epsilonEI = -0.2;       % From I to E
        epsilonII = -0.41514;       % From I to I
        
        m_H = -log(exp(-2.*m_tau) - (1 - exp(-thetaPhiE)).*epsilonEI);
        T = 2.*m_tau + thetaPhiE - m_H;
        
        f_comb(i, j) = 1/T;
        
        T_I = m_tau + thetaPhiI - (thetaPhiI/pi)*atan(tan(pi/thetaPhiI*m_tau)*exp(-2*pi*(-epsilonII)/thetaPhiI));
        T_E = 2*m_tau + thetaPhiE - (-log(exp(-2*m_tau) - (1 - exp(-thetaPhiE))*epsilonEI));
        if (T_E > T_I)
            f_comb(i, j) = NaN;            
        end
    end
end

tmp_PING = f_comb';
hSurface = surf(OneoverThetaPhiE_filter, OneoverThetaPhiI_filter, f_comb');
set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', [1 0 0]);

%% Plot PINGING
for i = 1:1:20
    file_txt = strcat(root_txt, '\PINGING\v0\PINGING_Types12Tau04BifDiagFullVaryPhiEPhiI', num2str(i), '.mat');
    
    load(file_txt);
    
    OneoverThetaPhiE_total = OneoverThetaPhiE_lin;
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

figure(1);hold on;xlabel('E');ylabel('I');
% % surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f1_total);
% surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f2_total);
% surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f3_total);
% % surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f4_total);
% % surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f5_total);
% surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f15_total);
% surf(OneoverThetaPhiI_total, OneoverThetaPhiE_total, f51_total);

% Filter
N_filter = 1;
i_filterE = 0;
[rowsE, colsE] = size(OneoverThetaPhiE_total);
[rowsI, colsI] = size(OneoverThetaPhiI_total);
for i = 1:1:colsE
    if (cmp(rem(i, N_filter), 0, 1e-6) == 0)
        if ((OneoverThetaPhiE_min <= OneoverThetaPhiE_total(1, i)) && (OneoverThetaPhiE_total(1, i) <= OneoverThetaPhiE_max))
            i_filterE = i_filterE + 1;
            OneoverThetaPhiE_filter(1, i_filterE) = OneoverThetaPhiE_total(1, i);
            
            i_filterI = 0;
            for j = 1:1:colsI
                if (cmp(rem(j, N_filter), 0, 1e-6) == 0)
                    if ((OneoverThetaPhiI_min <= OneoverThetaPhiI_total(1, j)) && (OneoverThetaPhiI_total(1, j) <= OneoverThetaPhiI_max))
                        i_filterI = i_filterI + 1;
                        OneoverThetaPhiI_filter(1, i_filterI) = OneoverThetaPhiI_total(1, j);
                        
                        f1_filter(i_filterE, i_filterI) = f1_total(i, j);
                        f2_filter(i_filterE, i_filterI) = f2_total(i, j);
                        f3_filter(i_filterE, i_filterI) = f3_total(i, j);
                        f4_filter(i_filterE, i_filterI) = f4_total(i, j);
                        f5_filter(i_filterE, i_filterI) = f5_total(i, j);
                        f51_filter(i_filterE, i_filterI) = f51_total(i, j);
                        f15_filter(i_filterE, i_filterI) = f15_total(i, j);
                        
                    end
                end
            end
            
        end
    end
end

% Combine
[rowsE, colsE] = size(OneoverThetaPhiE_filter);
[rowsI, colsI] = size(OneoverThetaPhiI_filter);
f_comb = NaN(colsE, colsI);

for i = 1:1:colsE
    for j = 1:1:colsI
        if (isnan(f2_filter(i, j)) == 0)
            f_comb(i, j) = f3_filter(i, j);
        elseif (isnan(f3_filter(i, j)) == 0)
            f_comb(i, j) = f3_filter(i, j);
        elseif (isnan(f51_filter(i, j)) == 0)
            f_comb(i, j) = f15_filter(i, j);
        else
            f_comb(i, j) = f15_filter(i, j);
        end
    end
end

% Make dark green

tmp_PINGING = f_comb';
tmp_PINGING_dark = tmp_PINGING;

[rows, cols] = size(tmp_PINGING);
for i=1:1:rows
    for j=1:1:cols
        if (tmp_PING(i, j) > tmp_ING(i, j))
            tmp_PINGING(i, j) = NaN;
        else
            tmp_PINGING_dark(i, j) = NaN;
        end
    end
end

hSurface = surf(OneoverThetaPhiE_filter, OneoverThetaPhiI_filter, tmp_PINGING);
set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);

hSurface = surf(OneoverThetaPhiE_filter, OneoverThetaPhiI_filter, tmp_PINGING_dark);
set(hSurface,'FaceColor',[101/255 163/255 92/255],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [101/255 163/255 92/255]);

xx = 0.628;

f2_filter(f2_filter < xx) = NaN;
f2_filter(3,25) = (f2_filter(2,25) + f2_filter(4,25))/2;
f2_filter(13,25) = (f2_filter(12,25) + f2_filter(14,25))/2;
f2_filter(7,25) = (f2_filter(6,25) + f2_filter(8,25))/2;
f2_filter(3,24) = (f2_filter(2,24) + f2_filter(4,24))/2;
f2_filter(5,23) = (f2_filter(4,23) + f2_filter(6,23))/2;
f2_filter(9,23) = (f2_filter(8,23) + f2_filter(10,23))/2;
f2_filter(5,22) = (f2_filter(4,22) + f2_filter(6,22))/2;
f2_filter(9,22) = (f2_filter(8,22) + f2_filter(10,22))/2;
f2_filter(2,21) = (f2_filter(1,21) + f2_filter(3,21))/2;

hSurface = surf(OneoverThetaPhiE_filter, OneoverThetaPhiI_filter, f2_filter');
set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', [0 1 0]);



box on
axis tight
axis square
view([-37.500000000000000, 30])

make_me_pretty(gcf, ...
    gca, 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'XTick',[0.71 0.73 0.75 0.77],'XTickLabel',{'';'';'';''});
set(gca,'YTick',[0.48 0.49 0.50 0.51 0.52],'YTickLabel',{'';'';'';'';''});
set(gca,'ZTick',[0.60 0.62 0.64 0.66],'ZTickLabel',{'';'';'';''});
xlabel('');ylabel('');

maximize_a_fig(gcf);
m_savefig('Types12Tau04BifDiagFullVaryPhiEPhiI_v2', 'eps');

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