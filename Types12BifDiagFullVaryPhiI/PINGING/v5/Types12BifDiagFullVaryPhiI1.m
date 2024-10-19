function Types12BifDiagFullVaryPhiE1
global epsilonEI epsilonIE epsilonII tau
global gammaE thetaVE thetaPhiE I_E
global gammaI thetaVI thetaPhiI I_I

cur_fileID = 1;

%% Type of cells
mode = [1 2];
N_deltaPsi = 10001;

ThetaPhiI_min = 1/0.54;
ThetaPhiI_max = 1/0.46;

ThetaPhiE_min = 1/0.85;
ThetaPhiE_max = 1/0.63;

%% Connections
tau = 0.4;

epsilonEI = -0.2;       % From I to E
epsilonIE = 0.5;       % From E to I
epsilonII = -0.41514;  % From I to I

%% E-cell parameters
phiE = NaN;
gammaE = 1;
thetaVE = 1; % Voltage threshold
thetaPhiE = 1/0.74;
I_E = (gammaE*thetaVE)./(1 - exp(-gammaE*thetaPhiE));   % Current to drive the voltage of the LIF from 0 to thetaVE within thetaPhiE

%% E-cell parameters
phiI = NaN;
gammaI = 1;
thetaVI = 1; % Voltage threshold

tmp_I_E = (gammaE*thetaVE)./(1 - exp(-gammaE*ThetaPhiE_max));   % Current to drive the voltage of the LIF from 0 to thetaVE within thetaPhiE
Tmp_ThetaE_min = getInvUm(tmp_I_E, gammaE, epsilonEI, ThetaPhiE_max, 1, NaN, thetaVE);

lbb_DeltaPsi = Tmp_ThetaE_min - ThetaPhiE_max;
ubb_DeltaPsi = ThetaPhiI_max;

OneoverThetaPhiI_wholelin = linspace(0.46, 0.54, 11);
N_OneoverThetaPhiI_lin = 100;
OneoverThetaPhiI_lin = linspace(OneoverThetaPhiI_wholelin(1, cur_fileID), OneoverThetaPhiI_wholelin(1, cur_fileID + 1), N_OneoverThetaPhiI_lin);

N_fixepoint = 5;    % Max number of fixed point. If f(1, i, N_fixepoint) is not NaN, Beware!!! the results are crap.

firstIte_FP_DeltaPsi = NaN(5, N_OneoverThetaPhiI_lin, N_fixepoint);
secondIte_FP_DeltaPsi = NaN(5, 5, N_OneoverThetaPhiI_lin, N_fixepoint);

f_i = 1;
for OneoverThetaPhiI = OneoverThetaPhiI_lin
    OneoverThetaPhiI
    
    thetaPhiI = 1/OneoverThetaPhiI;  % Phase threshold
    I_I = (gammaI*thetaVI)./(1 - exp(-gammaI*thetaPhiI));
    
    for which_sce = 1:1:5
        FP_DeltaPsi = Sces_firstiter(tau, ...
            phiE, gammaE, I_E, thetaVE, thetaPhiE, ...
            phiI, gammaI, I_I, thetaVI, thetaPhiI, ...
            epsilonEI, epsilonIE, epsilonII, ...
            mode, ...
            which_sce, ...
            lbb_DeltaPsi, ubb_DeltaPsi, N_deltaPsi);
        
        for i = 1:1:size(FP_DeltaPsi, 2)
            firstIte_FP_DeltaPsi(which_sce, f_i, i) = FP_DeltaPsi(1, i);
        end
    end
    
    for first_sce = 1:1:5
        for second_sce = 1:1:5
            FP_DeltaPsi = Sces_twotiter(tau, ...
                phiE, gammaE, I_E, thetaVE, thetaPhiE, ...
                phiI, gammaI, I_I, thetaVI, thetaPhiI, ...
                epsilonEI, epsilonIE, epsilonII, ...
                mode, ...
                first_sce, ...
                second_sce, ...
                lbb_DeltaPsi, ubb_DeltaPsi, N_deltaPsi);
            
            for i = 1:1:size(FP_DeltaPsi, 2)
                secondIte_FP_DeltaPsi(first_sce, second_sce, f_i, i) = FP_DeltaPsi(1, i);
            end
        end
    end
    
    f_i = f_i + 1;
    
    FN = strcat('PINGING_Types12Tau04BifDiagFullVaryPhiI', num2str(cur_fileID), '.mat');
    save(FN, 'phiI', 'phiE', ...
        'firstIte_FP_DeltaPsi', 'secondIte_FP_DeltaPsi', ...
        'lbb_DeltaPsi', 'ubb_DeltaPsi', ...
        'OneoverThetaPhiI_lin', ...
        'thetaVI', 'thetaPhiI', 'gammaI', 'I_I', ...
        'thetaVE', 'thetaPhiE', 'gammaE', 'I_E', ...
        'epsilonII', 'epsilonIE', 'epsilonEI', 'tau');
end

exit
end