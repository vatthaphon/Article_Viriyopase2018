function plotTypes11BifDiagFullVaryPhiE
clc;
clear all;
close all;

N_OneoverThetaPhiE = 100;
OneoverThetaPhiE_min = 0.42;
OneoverThetaPhiE_max = 0.52;

N_OneoverThetaPhiI = 100;
OneoverThetaPhiI_min = 0.46;
OneoverThetaPhiI_max = 0.582;

PINGING_FaceAlpha = 0.5 + 0.0;
PING_FaceAlpha = 0.7 - 0.0;
ING_FaceAlpha = 0.7 - 0.0;

figure(1); hold on;

OneoverThetaPhiE_lin = linspace(OneoverThetaPhiE_min, OneoverThetaPhiE_max, N_OneoverThetaPhiE);
OneoverThetaPhiI_lin = linspace(OneoverThetaPhiI_min, OneoverThetaPhiI_max, N_OneoverThetaPhiI);

m_tau = 0.4;
epsilonEI = -0.5;       % From I to E
epsilonII = -1.00;      % From I to I
for i = 1:1:N_OneoverThetaPhiE
    for j = 1:1:N_OneoverThetaPhiI
        thetaPhiE = 1/OneoverThetaPhiE_lin(1, i);
        thetaPhiI = 1/OneoverThetaPhiI_lin(1, j);
        
        %% PINGING
        epsilonIE = 0.1;       % From E to I
        
        val = freq_sce2(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        PINGING_sce2_f1(i, j) = val(1, 1);
        PINGING_sce2_f2(i, j) = val(1, 2);
        
        PINGING_f2(i, j) = val(1, 2);

        val = freq_sce3(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        PINGING_sce3_f1(i, j) = val(1, 1);
        PINGING_sce3_f2(i, j) = val(1, 2);
        
        if (isnan(PINGING_f2(i, j)) == 1)
            PINGING_f2(i, j) = PINGING_sce3_f2(i, j);
        end
        
        PINGING_sce4_f1(i, j) = freq_sce4(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        
        val = freq_sce15(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        PINGING_sce15_f1(i, j) = val(1, 1);
        PINGING_sce15_f2(i, j) = val(1, 2);

        %% ING
        epsilonIE = 0;       % From E to I
        
        val = freq_sce2(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        ING_sce2_f1(i, j) = val(1, 1);
        ING_sce2_f2(i, j) = val(1, 2);
        
        ING_f2(i, j) = val(1, 2);

        val = freq_sce3(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        ING_sce3_f1(i, j) = val(1, 1);
        ING_sce3_f2(i, j) = val(1, 2);
        
        if (isnan(ING_f2(i, j)) == 1)
            ING_f2(i, j) = ING_sce3_f2(i, j);
        end        
        
        val = freq_sce15(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        if (isnan(ING_sce2_f2(i, j)) == 1) && (isnan(ING_sce3_f2(i, j)) == 1)
            ING_f1(i, j) = val(1, 1);
            ING_f2(i, j) = val(1, 2);
        end        
        val = freq_sce15(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        ING_sce15_f1(i, j) = val(1, 1);
        ING_sce15_f2(i, j) = val(1, 2);
        
        %% PING
        val = freq_PING(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII);
        PING(i, j) = val;
    end
end

% Kill PING when PING overlap with PINGING
for i = 1:1:N_OneoverThetaPhiE
    for j = 1:1:N_OneoverThetaPhiI
        if (isnan(PINGING_sce4_f1(i, j)) == 0)
            PING(i, j) = NaN;
        end
    end
end

%% PING
% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PING);set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', 'none');
hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, PING');set(hSurface,'FaceColor',[1 0 0],'FaceAlpha', PING_FaceAlpha, 'EdgeColor', 'none');

%% PINGING
% % hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PINGING_sce2_f1);set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PINGING_sce2_f2);set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');
% 
% % hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PINGING_sce3_f1);set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PINGING_sce3_f2);set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');

% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PINGING_sce4_f1);set(hSurface,'FaceColor',[101/255 163/255 92/255],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, PINGING_sce4_f1');set(hSurface,'FaceColor',[101/255 163/255 92/255],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');
hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, PINGING_sce4_f1');set(hSurface,'FaceColor',[0/255 102/255 92/255],'FaceAlpha', 0.95, 'EdgeColor', 'none');

% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PINGING_sce15_f1);set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PINGING_sce15_f2);set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');

% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, PINGING_f2);set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, PINGING_f2');set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', PINGING_FaceAlpha, 'EdgeColor', 'none');
hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, PINGING_f2');set(hSurface,'FaceColor',[0 1 0],'FaceAlpha', 0.7, 'EdgeColor', 'none');


%% ING
% % hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, ING_sce2_f1);set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, ING_sce2_f2);set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
% 
% % hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, ING_sce3_f1);set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, ING_sce3_f2);set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
%         
% % hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, ING_sce15_f1);set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, ING_sce15_f2);set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');

% hSurface = surf(OneoverThetaPhiI_lin, OneoverThetaPhiE_lin, ING_f2);set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
% hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, ING_f2');set(hSurface,'FaceColor',[0 0 1],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');
hSurface = surf(OneoverThetaPhiE_lin, OneoverThetaPhiI_lin, ING_f2');set(hSurface,'FaceColor', [102/255 102/255 255/255],'FaceAlpha', ING_FaceAlpha, 'EdgeColor', 'none');

box on
axis tight
axis square
view([-37.500000000000000, 30])

make_me_pretty(gcf, ...
    gca, 40, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12, ...
    [], 12)

set(gca,'YTick',[0.48 0.52 0.56],'YTickLabel',{'';'';''});
set(gca,'XTick',[0.43 0.45 0.47 0.49 0.51],'XTickLabel',{'';'';'';'';''});
set(gca,'ZTick',[0.33 0.35 0.37 0.39 0.41],'ZTickLabel',{'';'';'';'';''});

% xlabel('E');
% ylabel('I');
maximize_a_fig(gcf);
% m_savefig('Types11BifDiagFullVaryPhiEPhiI', 'eps');
m_savefig('Types11BifDiagFullVaryPhiEPhiI_v2', 'eps');

end

function val = freq_PING(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII)
    T_I = thetaPhiI + m_tau - H_LIF(m_tau, -epsilonII, thetaPhiI);
    T_E = 2*m_tau + thetaPhiE - H_LIF(2*m_tau, -epsilonEI, thetaPhiE);
    
    if (T_E <= T_I)
        A = (1 - exp(-thetaPhiE))*epsilonEI;
        val = 1/(2*m_tau + thetaPhiE + log(exp(-2*m_tau) - A));        
    else
        val = NaN;
    end
end

function val = GAMMA(thetaPhi, epsi)
val = (1 - exp(-thetaPhi))*epsi;
end

function val = freq_sce4(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII)
DeltaTheta = thetaPhiE - thetaPhiI;
H_I = H_LIF(thetaPhiI, -epsilonIE, thetaPhiI);
DeltaPhi = log((exp(-m_tau) - GAMMA(thetaPhiI, epsilonII))/(exp(-2*m_tau) - GAMMA(thetaPhiE, epsilonEI)));

if (((m_tau + DeltaTheta) <= DeltaPhi) && (DeltaPhi <= (DeltaTheta + thetaPhiI + m_tau - H_I)))
    A = (1 - exp(-thetaPhiE))*epsilonEI;
    val = 1/(2*m_tau + thetaPhiE + log(exp(-2*m_tau) - A));    
else
    val = NaN;
end
end

function val = freq_sce3(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII)
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

if (DeltaTheta <= DeltaPhi1) && (DeltaPhi1 <= (DeltaTheta + m_tau))
    T1 = thetaPhiE + log(exp(-(m_tau + DeltaPhi1 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI) + m_tau + DeltaPhi1 - DeltaTheta;

    val(1, 1) = 1/T1;
else
    val(1, 1) = NaN;
end

if (DeltaTheta <= DeltaPhi2) && (DeltaPhi2 <= (DeltaTheta + m_tau))
    T2 = thetaPhiE + log(exp(-(m_tau + DeltaPhi2 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI) + m_tau + DeltaPhi2 - DeltaTheta;

    val(1, 2) = 1/T2;
else
    val(1, 2) = NaN;
end

end

function val = freq_sce2(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII)
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

if ((DeltaTheta - m_tau ) <= DeltaPhi1) && (DeltaPhi1 <= DeltaTheta)
    T1 = m_tau + PhiE1 + log(exp(-(m_tau + DeltaPhi1 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI);

    val(1, 1) = 1/T1;
else
    val(1, 1) = NaN;
end

if ((DeltaTheta - m_tau ) <= DeltaPhi2) && (DeltaPhi2 <= DeltaTheta)
    T2 = m_tau + PhiE2 + log(exp(-(m_tau + DeltaPhi2 - DeltaTheta)) - (1 - exp(-thetaPhiE)).*epsilonEI);

    val(1, 2) = 1/T2;
else
    val(1, 2) = NaN;
end

end

function val = freq_sce15(m_tau, thetaPhiI, thetaPhiE, epsilonEI, epsilonIE, epsilonII)
DeltaTheta = thetaPhiE - thetaPhiI;

% % sce. 1
% phiE = 0.977;
% phiI = thetaPhiI;
% deltaphi = phiE - phiI;
% delta_xi = deltaphi - DeltaTheta
% delta_xi = H_LIF(thetaPhiE + delta_xi + m_tau, epsilonEI, thetaPhiE) - (H_LIF(m_tau, epsilonII, thetaPhiI) + DeltaTheta);
% 
% delta_xi + DeltaTheta + (-0.435);
% % sce. 5
% phiE = thetaPhiE;
% phiI = 1.09;
% deltaphi = phiE - phiI;
% delta_xi = deltaphi - DeltaTheta;
% delta_xi = m_tau - H_LIF(thetaPhiI - delta_xi + m_tau, epsilonIE, thetaPhiI) - DeltaTheta
% 
% delta_xi + DeltaTheta + thetaPhiI

A = exp(DeltaTheta - m_tau);
B = exp(-thetaPhiI - m_tau);
C = H_LIF(m_tau, epsilonII, thetaPhiI) + DeltaTheta;
D = exp(-thetaPhiE - m_tau);
E = B*exp(-C);
GAMMA1 = GAMMA(thetaPhiI, epsilonIE);
GAMMA2 = GAMMA(thetaPhiE, epsilonEI);

X = A*GAMMA2;
Y = E-A*D+GAMMA1*GAMMA2;
Z = GAMMA1*D;

DeltaXi1 = log((-Y + sqrt(Y*Y + 4*X*Z))/(2*X));
DeltaXi2 = log((-Y - sqrt(Y*Y + 4*X*Z))/(2*X));

% Check sce. 1
if (DeltaXi1 <= (-m_tau))
   
    DeltaXi1 = H_LIF(thetaPhiE + DeltaXi1 + m_tau, epsilonEI, thetaPhiE) - (H_LIF(m_tau, epsilonII, thetaPhiI) + DeltaTheta);
    
    % Check sce. 5
    if ((thetaPhiI + m_tau - H_LIF(thetaPhiI, -epsilonIE, thetaPhiI)) < DeltaXi1)
        DeltaPhi1 = DeltaXi1 + DeltaTheta;
        
        phiI = thetaPhiE - DeltaPhi1;
        A = phiI - H_LIF(m_tau, epsilonII, thetaPhiI);
        B = thetaPhiI - H_LIF(phiI + m_tau, epsilonIE, thetaPhiI);
        T = m_tau + A + m_tau + B;
        
        val(1, 1) = 1/T;    
    else
        val(1, 1) = NaN;    
    end
else
    val(1, 1) = NaN;
end

if (DeltaXi2 <= (-m_tau))
    
    DeltaXi1 = DeltaXi2;    
    DeltaXi1 = H_LIF(thetaPhiE + DeltaXi1 + m_tau, epsilonEI, thetaPhiE) - (H_LIF(m_tau, epsilonII, thetaPhiI) + DeltaTheta);    
    
    % Check sce. 5
    if ((thetaPhiI + m_tau - H_LIF(thetaPhiI, -epsilonIE, thetaPhiI)) < DeltaXi1)
        DeltaPhi1 = DeltaXi1 + DeltaTheta;
        
        phiI = thetaPhiE - DeltaPhi1;
        A = phiI - H_LIF(m_tau, epsilonII, thetaPhiI);
        B = thetaPhiI - H_LIF(phiI + m_tau, epsilonIE, thetaPhiI);
        T = m_tau + A + m_tau + B;
        
        val(1, 2) = 1/T;    
    else
        val(1, 2) = NaN;    
    end   
else
    val(1, 2) = NaN;
end

end

function val = H_LIF(phi, epsi, Theta)
val = -log(exp(-phi) - (1 - exp(-Theta))*epsi);
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