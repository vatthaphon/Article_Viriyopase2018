function Types12BifDiagFullVaryPhiI1
global epsilonEI epsilonIE epsilonII tau
global gammaE thetaVE thetaPhiE I_E
global gammaI thetaVI thetaPhiI I_I

%% Type of cells
mode = [1 2];
N_deltaPsi = 10001;

%% Connections
tau = 0.4;

epsilonEI = -0.2;       % From I to E

epsilonIE = 0;       % From E to I
% epsilonIE = 0.1;       % From E to I
% epsilonIE = Inf;
% epsilonIE = 5;

% epsilonII = -0.41514;  % From I to I
N_epsilonII_lin = 100;
epsilonII_lin = linspace(-0.6, -0.5, N_epsilonII_lin);

%% E-cell parameters
phiE = NaN;
gammaE = 1;
thetaVE = 1; % Voltage threshold
thetaPhiE = 1/0.7455;
I_E = (gammaE*thetaVE)./(1 - exp(-gammaE*thetaPhiE));   % Current to drive the voltage of the LIF from 0 to thetaVE within thetaPhiE

%% I-cell parameters
phiI = NaN;
gammaI = 1;
thetaVI = 1; % Voltage threshold
thetaPhiI = 1/0.5;  % Phase threshold
I_I = (gammaI*thetaVI)./(1 - exp(-gammaI*thetaPhiI));

f1 = NaN(1, N_epsilonII_lin);
f2 = NaN(1, N_epsilonII_lin);
f3 = NaN(1, N_epsilonII_lin);
f4 = NaN(1, N_epsilonII_lin);
f5 = NaN(1, N_epsilonII_lin);

f_i = 0;
for epsilonII = epsilonII_lin
    epsilonII
        
    min_deltaPsi_Sce1 = 1000;   min_deltaPsi_Sce1_exist = 0; min_deltaPsi_Sce1_val = NaN;
    max_deltaPsi_Sce1 = -1000;  max_deltaPsi_Sce1_exist = 0; max_deltaPsi_Sce1_val = NaN;
    min_deltaPsi_Sce2 = 1000;   min_deltaPsi_Sce2_exist = 0; min_deltaPsi_Sce2_val = NaN;
    max_deltaPsi_Sce2 = -1000;  max_deltaPsi_Sce2_exist = 0; max_deltaPsi_Sce2_val = NaN;
    min_deltaPsi_Sce3 = 1000;   min_deltaPsi_Sce3_exist = 0; min_deltaPsi_Sce3_val = NaN;
    max_deltaPsi_Sce3 = -1000;  max_deltaPsi_Sce3_exist = 0; max_deltaPsi_Sce3_val = NaN;
    min_deltaPsi_Sce4 = 1000;   min_deltaPsi_Sce4_exist = 0; min_deltaPsi_Sce4_val = NaN;
    max_deltaPsi_Sce4 = -1000;  max_deltaPsi_Sce4_exist = 0; max_deltaPsi_Sce4_val = NaN;
    min_deltaPsi_Sce5 = 1000;   min_deltaPsi_Sce5_exist = 0; min_deltaPsi_Sce5_val = NaN;
    max_deltaPsi_Sce5 = -1000;  max_deltaPsi_Sce5_exist = 0; max_deltaPsi_Sce5_val = NaN;

    min_deltaPsi_Sce51 = 1000;   min_deltaPsi_Sce51_exist = 0; min_deltaPsi_Sce51_val = NaN;
    max_deltaPsi_Sce51 = -1000;  max_deltaPsi_Sce51_exist = 0; max_deltaPsi_Sce51_val = NaN;
    min_deltaPsi_Sce15 = 1000;   min_deltaPsi_Sce15_exist = 0; min_deltaPsi_Sce15_val = NaN;
    max_deltaPsi_Sce15 = -1000;  max_deltaPsi_Sce15_exist = 0; max_deltaPsi_Sce15_val = NaN;    
    
    
    best_deltaPsi_Sce1 = NaN; delta_deltaPsi_Sce1 = 1000;
    best_deltaPsi_Sce2 = NaN; delta_deltaPsi_Sce2 = 1000;
    best_deltaPsi_Sce3 = NaN; delta_deltaPsi_Sce3 = 1000;
    best_deltaPsi_Sce4 = NaN; delta_deltaPsi_Sce4 = 1000;
    best_deltaPsi_Sce5 = NaN; delta_deltaPsi_Sce5 = 1000;
    
    best_deltaPsi_Sce51 = NaN; delta_deltaPsi_Sce51 = 1000;
    best_deltaPsi_Sce15 = NaN; delta_deltaPsi_Sce15 = 1000;
       
    deltaPsi_lin = linspace(-2, 2, N_deltaPsi);
    
    for deltaPsi = deltaPsi_lin
        
        ret = Sce1(tau, deltaPsi, phiE, gammaE, I_E, thetaVE, thetaPhiE, phiI, gammaI, I_I, thetaVI, thetaPhiI, epsilonEI, epsilonIE, epsilonII, mode);
        if (isnan(ret) == 0)
            if (deltaPsi < min_deltaPsi_Sce1)
                min_deltaPsi_Sce1 = deltaPsi;
                min_deltaPsi_Sce1_exist = 1;
                min_deltaPsi_Sce1_val = ret;
            end
            
            if (deltaPsi > max_deltaPsi_Sce1)
                max_deltaPsi_Sce1 = deltaPsi;
                max_deltaPsi_Sce1_exist = 1;
                max_deltaPsi_Sce1_val = ret;
            end
            
            if (abs(ret - deltaPsi) < delta_deltaPsi_Sce1)
                delta_deltaPsi_Sce1 = abs(ret - deltaPsi);
                best_deltaPsi_Sce1 = deltaPsi;
            end
        end
        
        ret = Sce5(tau, ret, phiE, gammaE, I_E, thetaVE, thetaPhiE, phiI, gammaI, I_I, thetaVI, thetaPhiI, epsilonEI, epsilonIE, epsilonII, mode);
        if (isnan(ret) == 0)
            if (deltaPsi < min_deltaPsi_Sce15)
                min_deltaPsi_Sce15 = deltaPsi;
                min_deltaPsi_Sce15_exist = 1;
                min_deltaPsi_Sce15_val = ret;
            end
            
            if (deltaPsi > max_deltaPsi_Sce15)
                max_deltaPsi_Sce15 = deltaPsi;
                max_deltaPsi_Sce15_exist = 1;
                max_deltaPsi_Sce15_val = ret;
            end
            
            if (abs(ret - deltaPsi) < delta_deltaPsi_Sce15)
                delta_deltaPsi_Sce15 = abs(ret - deltaPsi);
                best_deltaPsi_Sce15 = deltaPsi;
            end
        end
        
        
        ret = Sce2(tau, deltaPsi, phiE, gammaE, I_E, thetaVE, thetaPhiE, phiI, gammaI, I_I, thetaVI, thetaPhiI, epsilonEI, epsilonIE, epsilonII, mode);
        if (isnan(ret) == 0)
            if (deltaPsi < min_deltaPsi_Sce2)
                min_deltaPsi_Sce2 = deltaPsi;
                min_deltaPsi_Sce2_exist = 1;
                min_deltaPsi_Sce2_val = ret;
            end
            
            if (deltaPsi > max_deltaPsi_Sce2)
                max_deltaPsi_Sce2 = deltaPsi;
                max_deltaPsi_Sce2_exist = 1;
                max_deltaPsi_Sce2_val = ret;
            end
            
            if (abs(ret - deltaPsi) < delta_deltaPsi_Sce2)
                delta_deltaPsi_Sce2 = abs(ret - deltaPsi);
                best_deltaPsi_Sce2 = deltaPsi;
            end
        end
        
        ret = Sce3(tau, deltaPsi, phiE, gammaE, I_E, thetaVE, thetaPhiE, phiI, gammaI, I_I, thetaVI, thetaPhiI, epsilonEI, epsilonIE, epsilonII, mode);
        if (isnan(ret) == 0)
            if (deltaPsi < min_deltaPsi_Sce3)
                min_deltaPsi_Sce3 = deltaPsi;
                min_deltaPsi_Sce3_exist = 1;
                min_deltaPsi_Sce3_val = ret;
            end
            
            if (deltaPsi > max_deltaPsi_Sce3)
                max_deltaPsi_Sce3 = deltaPsi;
                max_deltaPsi_Sce3_exist = 1;
                max_deltaPsi_Sce3_val = ret;
            end
            
            if (abs(ret - deltaPsi) < delta_deltaPsi_Sce3)
                delta_deltaPsi_Sce3 = abs(ret - deltaPsi);
                best_deltaPsi_Sce3 = deltaPsi;
            end
        end
        
        ret = Sce4(tau, deltaPsi, phiE, gammaE, I_E, thetaVE, thetaPhiE, phiI, gammaI, I_I, thetaVI, thetaPhiI, epsilonEI, epsilonIE, epsilonII, mode);
        if (isnan(ret) == 0)
            if (deltaPsi < min_deltaPsi_Sce4)
                min_deltaPsi_Sce4 = deltaPsi;
                min_deltaPsi_Sce4_exist = 1;
                min_deltaPsi_Sce4_val = ret;
            end
            
            if (deltaPsi > max_deltaPsi_Sce4)
                max_deltaPsi_Sce4 = deltaPsi;
                max_deltaPsi_Sce4_exist = 1;
                max_deltaPsi_Sce4_val = ret;
            end
            
            if (abs(ret - deltaPsi) < delta_deltaPsi_Sce4)
                delta_deltaPsi_Sce4 = abs(ret - deltaPsi);
                best_deltaPsi_Sce4 = deltaPsi;
            end
        end
        
        ret = Sce5(tau, deltaPsi, phiE, gammaE, I_E, thetaVE, thetaPhiE, phiI, gammaI, I_I, thetaVI, thetaPhiI, epsilonEI, epsilonIE, epsilonII, mode);
        if (isnan(ret) == 0)
            if (deltaPsi < min_deltaPsi_Sce5)
                min_deltaPsi_Sce5 = deltaPsi;
                min_deltaPsi_Sce5_exist = 1;
                min_deltaPsi_Sce5_val = ret;
            end
            
            if (deltaPsi > max_deltaPsi_Sce5)
                max_deltaPsi_Sce5 = deltaPsi;
                max_deltaPsi_Sce5_exist = 1;
                max_deltaPsi_Sce5_val = ret;
            end
            
            if (abs(ret - deltaPsi) < delta_deltaPsi_Sce5)
                delta_deltaPsi_Sce5 = abs(ret - deltaPsi);
                best_deltaPsi_Sce5 = deltaPsi;
            end
        end
        
        ret = Sce1(tau, ret, phiE, gammaE, I_E, thetaVE, thetaPhiE, phiI, gammaI, I_I, thetaVI, thetaPhiI, epsilonEI, epsilonIE, epsilonII, mode);
        if (isnan(ret) == 0)
            if (deltaPsi < min_deltaPsi_Sce51)
                min_deltaPsi_Sce51 = deltaPsi;
                min_deltaPsi_Sce51_exist = 1;
                min_deltaPsi_Sce51_val = ret;
            end
            
            if (deltaPsi > max_deltaPsi_Sce51)
                max_deltaPsi_Sce51 = deltaPsi;
                max_deltaPsi_Sce51_exist = 1;
                max_deltaPsi_Sce51_val = ret;
            end
            
            if (abs(ret - deltaPsi) < delta_deltaPsi_Sce51)
                delta_deltaPsi_Sce51 = abs(ret - deltaPsi);
                best_deltaPsi_Sce51 = deltaPsi;
            end
        end        
    end
    
    % Determine range
    f_i = f_i + 1;
    if ((((min_deltaPsi_Sce1 <= min_deltaPsi_Sce1_val) && (max_deltaPsi_Sce1 >= max_deltaPsi_Sce1_val)) || ...
            ((min_deltaPsi_Sce1 >= min_deltaPsi_Sce1_val) && (max_deltaPsi_Sce1 <= max_deltaPsi_Sce1_val))) && ...
            (min_deltaPsi_Sce1_exist) && (max_deltaPsi_Sce1_exist))
        f1(1, f_i) = calTTypes12(best_deltaPsi_Sce1);
    else
        f1(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce15 <= min_deltaPsi_Sce15_val) && (max_deltaPsi_Sce15 >= max_deltaPsi_Sce15_val)) || ...
            ((min_deltaPsi_Sce15 >= min_deltaPsi_Sce15_val) && (max_deltaPsi_Sce15 <= max_deltaPsi_Sce15_val))) && ...
            (min_deltaPsi_Sce15_exist) && (max_deltaPsi_Sce15_exist))
        f15(1, f_i) = calTTypes12(best_deltaPsi_Sce15);
    else
        f15(1, f_i) = NaN;
    end    
    
    if ((((min_deltaPsi_Sce2 <= min_deltaPsi_Sce2_val) && (max_deltaPsi_Sce2 >= max_deltaPsi_Sce2_val)) || ...
            ((min_deltaPsi_Sce2 >= min_deltaPsi_Sce2_val) && (max_deltaPsi_Sce2 <= max_deltaPsi_Sce2_val))) && ...
            (min_deltaPsi_Sce2_exist) && (max_deltaPsi_Sce2_exist))
        f2(1, f_i) = calTTypes12(best_deltaPsi_Sce2);
    else
        f2(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce3 <= min_deltaPsi_Sce3_val) && (max_deltaPsi_Sce3 >= max_deltaPsi_Sce3_val)) || ...
            ((min_deltaPsi_Sce3 >= min_deltaPsi_Sce3_val) && (max_deltaPsi_Sce3 <= max_deltaPsi_Sce3_val))) && ...
            (min_deltaPsi_Sce3_exist) && (max_deltaPsi_Sce3_exist))
        f3(1, f_i) = calTTypes12(best_deltaPsi_Sce3);
    else
        f3(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce4 <= min_deltaPsi_Sce4_val) && (max_deltaPsi_Sce4 >= max_deltaPsi_Sce4_val)) || ...
            ((min_deltaPsi_Sce4 >= min_deltaPsi_Sce4_val) && (max_deltaPsi_Sce4 <= max_deltaPsi_Sce4_val))) && ...
            (min_deltaPsi_Sce4_exist) && (max_deltaPsi_Sce4_exist))
        f4(1, f_i) = calTTypes12(best_deltaPsi_Sce4);
    else
        f4(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce5 <= min_deltaPsi_Sce5_val) && (max_deltaPsi_Sce5 >= max_deltaPsi_Sce5_val)) || ...
            ((min_deltaPsi_Sce5 >= min_deltaPsi_Sce5_val) && (max_deltaPsi_Sce5 <= max_deltaPsi_Sce5_val))) && ...
            (min_deltaPsi_Sce5_exist) && (max_deltaPsi_Sce5_exist))
        f5(1, f_i) = calTTypes12(best_deltaPsi_Sce5);
    else
        f5(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce51 <= min_deltaPsi_Sce51_val) && (max_deltaPsi_Sce51 >= max_deltaPsi_Sce51_val)) || ...
            ((min_deltaPsi_Sce51 >= min_deltaPsi_Sce51_val) && (max_deltaPsi_Sce51 <= max_deltaPsi_Sce51_val))) && ...
            (min_deltaPsi_Sce51_exist) && (max_deltaPsi_Sce51_exist))
        f51(1, f_i) = calTTypes12(best_deltaPsi_Sce51);
    else
        f51(1, f_i) = NaN;
    end     
    
    save('ING_Types12Tau04BifDiagFullVarygI2I4.mat', 'phiI', 'phiE', 'f51', 'f15', 'f5', 'f4', 'f3', 'f2', 'f1', 'epsilonII_lin', 'thetaVI', 'gammaI', 'I_E', 'thetaPhiE', 'thetaVE', 'gammaE', 'epsilonIE', 'epsilonEI', 'tau');
end

exit
end

function f = calTTypes12(deltaPsi)
global NoConn neuronModel
global epsilonEI epsilonIE epsilonII tau
global gammaE thetaPhiE I_E
global gammaI thetaPhiI I_I
global thetaVE thetaVI
NoConn=-1;

% 0: the Mirollo Strogatz neuron
% 1: the Leaky Integrate-and-Fire neuron
% 2: the SINE neuron
neuronModel = [1 2];

tol_eq = 1e-6;
tEnd=100;
dt = 0.0001;

% Column = from, Row = to. 1 = E, 2 = I.
epsiMat=[
    0.0  epsilonEI
    epsilonIE  epsilonII
    ];

deltaTheta = thetaPhiE - thetaPhiI;
deltaPhi = deltaPsi + deltaTheta;
phiE = 0.0;
phiI = phiE - deltaPhi;
initPhase = [phiE phiI];
tauMat=[
    NaN   tau
    tau   tau
    ];

tBegin=0;

I_i=[I_E I_I];
osc_types=[1 -1];   % -1=inh., 1=exc.

%% Pre-processing
k_i = [1 1];
thetaPhi_i = [thetaPhiE thetaPhiI];
gamma_i=[gammaE gammaI];
thetaV_i = [thetaVE thetaVI];

%% Begin processing
[rec_phases spikeTimes discrete_t] = simNetworkInPhase(osc_types, epsiMat, tauMat, initPhase, tBegin, tEnd, dt, gamma_i, I_i, tol_eq, thetaPhi_i, k_i, thetaV_i);

E_spks = spikeTimes(isnan(spikeTimes(:, 1)) == 0, 1);
I_spks = spikeTimes(isnan(spikeTimes(:, 2)) == 0, 1);

if (size(E_spks, 1) >= 2)
    T = E_spks(end, 1) - E_spks(end - 1, 1);
else
    if (size(I_spks, 1) >= 2)
        T = I_spks(end, 1) - I_spks(end - 1, 1);
    else
        T = NaN;
    end
end


f = 1/T;
end