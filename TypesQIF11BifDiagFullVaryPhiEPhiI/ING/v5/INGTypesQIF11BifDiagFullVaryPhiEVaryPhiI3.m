function TypesQIF11BifDiagFullVaryPhiEVaryPhiI1
global epsilonEI epsilonIE epsilonII tau
global gammaE thetaVE thetaPhiE I_E
global gammaI thetaVI thetaPhiI I_I

cur_fileID = 3;

%% Type of cells
mode = [3 3];
N_deltaPsi = 10001;

ThetaPhiI_min = 1/(0.4401 + 0.0018);
ThetaPhiI_max = 1/(0.4401 - 0.0018);

ThetaPhiE_min = 1/(0.4744 + 0.0018);
ThetaPhiE_max = 1/(0.4744 - 0.0018);

N_OneoverThetaPhiE_lin = 10;
N_OneoverThetaPhiI_lin = 100;

% N_OneoverThetaPhiE_lin = 2;
% N_OneoverThetaPhiI_lin = 10;

%% Connections
tau = 0.4;

epsilonEI = -0.9;       % From I to E
epsilonIE = 0.0;       % From E to I
epsilonII = -1.0051;  % From I to I

%% E-cell parameters
phiE = NaN;
gammaE = 1;
thetaVE = 1; % Voltage threshold

%% I-cell parameters
phiI = NaN;
gammaI = 1;
thetaVI = 1; % Voltage threshold

lbb_DeltaPsi = -ThetaPhiE_max;
ubb_DeltaPsi = ThetaPhiI_max;

OneoverThetaPhiE_wholelin = linspace(1/ThetaPhiE_max, 1/ThetaPhiE_min, 11);
OneoverThetaPhiE_lin = linspace(OneoverThetaPhiE_wholelin(1, cur_fileID), OneoverThetaPhiE_wholelin(1, cur_fileID + 1), N_OneoverThetaPhiE_lin);

N_fixepoint = 5;    % Thsi can be dec. % Max number of fixed point. If f(1, i, N_fixepoint) is not NaN, Beware!!! the results are crap.

FN = strcat('ING_TypesQIF11Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(cur_fileID), '.mat');
if (exist(FN, 'file') == 2)
    load(FN);
    Begin_E_f_i = E_f_i + 1;
else
    firstIte_FP_DeltaPsi = NaN(5, N_OneoverThetaPhiE_lin, N_OneoverThetaPhiI_lin, N_fixepoint);
    secondIte_FP_DeltaPsi = NaN(5, 5, N_OneoverThetaPhiE_lin, N_OneoverThetaPhiI_lin, N_fixepoint);

    Begin_E_f_i = 1;    
end

for E_f_i = Begin_E_f_i:1:N_OneoverThetaPhiE_lin
    OneoverThetaPhiE = OneoverThetaPhiE_lin(1, E_f_i);
    thetaPhiE = 1/OneoverThetaPhiE;  % Phase threshold
    I_E = (pi*pi)/(thetaPhiE*thetaPhiE);
        
    OneoverThetaPhiI_lin = linspace(1/ThetaPhiI_max, 1/ThetaPhiI_min, N_OneoverThetaPhiI_lin);    
    
    f_i = 1;
    for OneoverThetaPhiI = OneoverThetaPhiI_lin
        [status, result] = system(strcat('top -n 1 -b -u viriyopa > ING_TypesQIF11Tau04BifDiagFullVaryPhiEVaryPhiI', num2str(cur_fileID), '.log'));
        display(strcat(num2str(OneoverThetaPhiE),',',num2str(OneoverThetaPhiI)));
        
        thetaPhiI = 1/OneoverThetaPhiI;  % Phase threshold
        I_I = (pi*pi)/(thetaPhiI*thetaPhiI);
        
        for which_sce = 1:1:5
            FP_DeltaPsi = Sces_firstiter(tau, ...
                phiE, gammaE, I_E, thetaVE, thetaPhiE, ...
                phiI, gammaI, I_I, thetaVI, thetaPhiI, ...
                epsilonEI, epsilonIE, epsilonII, ...
                mode, ...
                which_sce, ...
                lbb_DeltaPsi, ubb_DeltaPsi, N_deltaPsi);
            
            for i = 1:1:size(FP_DeltaPsi, 2)
                firstIte_FP_DeltaPsi(which_sce, E_f_i, f_i, i) = FP_DeltaPsi(1, i); %% We can throw exceptions if there are more fixed points.
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
                    secondIte_FP_DeltaPsi(first_sce, second_sce, E_f_i, f_i, i) = FP_DeltaPsi(1, i);
                end
            end
        end
        
        f_i = f_i + 1;
        
    end
    
    save(FN, 'phiI', 'phiE', ...
        'firstIte_FP_DeltaPsi', 'secondIte_FP_DeltaPsi', ...
        'lbb_DeltaPsi', 'ubb_DeltaPsi', ...
        'OneoverThetaPhiE_lin', 'OneoverThetaPhiI_lin', ...
        'thetaVI', 'thetaPhiI', 'gammaI', 'I_I', ...
        'thetaVE', 'thetaPhiE', 'gammaE', 'I_E', ...
        'epsilonII', 'epsilonIE', 'epsilonEI', 'tau', ...
        'E_f_i');
end

exit
end