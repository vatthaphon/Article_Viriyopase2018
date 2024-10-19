function Types11BifDiagFullVaryPhiI6
global epsilonEI epsilonIE epsilonII tau
global gammaE thetaVE thetaPhiE I_E
global gammaI thetaVI thetaPhiI I_I

%% Type of cells
mode = [1 1];

%% Connections
tau = 0.4;
epsilonEI = -0.5;       % From I to E
% epsilonIE = 0;       % From E to I
epsilonIE = 0.1;       % From E to I
% epsilonIE = Inf;       % From E to I
% epsilonIE = 5;       % From E to I
epsilonII = -1.00;  % From I to I

%% E-cell parameters
phiE = NaN;
gammaE = 1;
thetaVE = 1; % Voltage threshold
thetaPhiE = 1/0.495;   % Phase threshold
I_E = (gammaE*thetaVE)./(1 - exp(-gammaE*thetaPhiE));   % Current to drive the voltage of the LIF from 0 to thetaVE within thetaPhiE

%% I-cell parameters
phiI = NaN;
gammaI = 1;
thetaVI = 1; % Voltage threshold

N_OneoverThetaPhiI_lin = 100;
OneoverThetaPhiI_lin = linspace(0.500, 0.510, N_OneoverThetaPhiI_lin);

f1 = NaN(1, N_OneoverThetaPhiI_lin);
f2 = NaN(1, N_OneoverThetaPhiI_lin);
f3 = NaN(1, N_OneoverThetaPhiI_lin);
f4 = NaN(1, N_OneoverThetaPhiI_lin);
f5 = NaN(1, N_OneoverThetaPhiI_lin);

f_i = 0;
for OneoverThetaPhiI = OneoverThetaPhiI_lin
    OneoverThetaPhiI
    
    thetaPhiI = 1/OneoverThetaPhiI;  % Phase threshold
    I_I = (gammaI*thetaVI)./(1 - exp(-gammaI*thetaPhiI));    
    
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
    
    best_deltaPsi_Sce1 = NaN; delta_deltaPsi_Sce1 = 1000;
    best_deltaPsi_Sce2 = NaN; delta_deltaPsi_Sce2 = 1000;
    best_deltaPsi_Sce3 = NaN; delta_deltaPsi_Sce3 = 1000;
    best_deltaPsi_Sce4 = NaN; delta_deltaPsi_Sce4 = 1000;
    best_deltaPsi_Sce5 = NaN; delta_deltaPsi_Sce5 = 1000;
    
    N_deltaPsi = 10001;
    deltaPsi_lin = linspace(-20, 20, N_deltaPsi);
    
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
    end
    
    % Determine range
    f_i = f_i + 1;
    if ((((min_deltaPsi_Sce1 <= min_deltaPsi_Sce1_val) && (max_deltaPsi_Sce1 >= max_deltaPsi_Sce1_val)) || ...
            ((min_deltaPsi_Sce1 >= min_deltaPsi_Sce1_val) && (max_deltaPsi_Sce1 <= max_deltaPsi_Sce1_val))) && ...
            (min_deltaPsi_Sce1_exist) && (max_deltaPsi_Sce1_exist))
        f1(1, f_i) = calTTypes11(best_deltaPsi_Sce1);
    else
        f1(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce2 <= min_deltaPsi_Sce2_val) && (max_deltaPsi_Sce2 >= max_deltaPsi_Sce2_val)) || ...
            ((min_deltaPsi_Sce2 >= min_deltaPsi_Sce2_val) && (max_deltaPsi_Sce2 <= max_deltaPsi_Sce2_val))) && ...
            (min_deltaPsi_Sce2_exist) && (max_deltaPsi_Sce2_exist))
        f2(1, f_i) = calTTypes11(best_deltaPsi_Sce2);
    else
        f2(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce3 <= min_deltaPsi_Sce3_val) && (max_deltaPsi_Sce3 >= max_deltaPsi_Sce3_val)) || ...
            ((min_deltaPsi_Sce3 >= min_deltaPsi_Sce3_val) && (max_deltaPsi_Sce3 <= max_deltaPsi_Sce3_val))) && ...
            (min_deltaPsi_Sce3_exist) && (max_deltaPsi_Sce3_exist))
        f3(1, f_i) = calTTypes11(best_deltaPsi_Sce3);
    else
        f3(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce4 <= min_deltaPsi_Sce4_val) && (max_deltaPsi_Sce4 >= max_deltaPsi_Sce4_val)) || ...
            ((min_deltaPsi_Sce4 >= min_deltaPsi_Sce4_val) && (max_deltaPsi_Sce4 <= max_deltaPsi_Sce4_val))) && ...
            (min_deltaPsi_Sce4_exist) && (max_deltaPsi_Sce4_exist))
        f4(1, f_i) = calTTypes11(best_deltaPsi_Sce4);
    else
        f4(1, f_i) = NaN;
    end
    
    if ((((min_deltaPsi_Sce5 <= min_deltaPsi_Sce5_val) && (max_deltaPsi_Sce5 >= max_deltaPsi_Sce5_val)) || ...
            ((min_deltaPsi_Sce5 >= min_deltaPsi_Sce5_val) && (max_deltaPsi_Sce5 <= max_deltaPsi_Sce5_val))) && ...
            (min_deltaPsi_Sce5_exist) && (max_deltaPsi_Sce5_exist))
        f5(1, f_i) = calTTypes11(best_deltaPsi_Sce5);
    else
        f5(1, f_i) = NaN;
    end
    
    save('PINGING_Types11BifDiagFullVaryPhiI6.mat', 'phiI', 'phiE', 'f5', 'f4', 'f3', 'f2', 'f1', 'OneoverThetaPhiI_lin', 'thetaVI', 'gammaI', 'I_E', 'thetaPhiE', 'thetaVE', 'gammaE', 'epsilonII', 'epsilonIE', 'epsilonEI', 'tau');
end

exit
end

function f = calTTypes11(deltaPsi)
global NoConn neuronModel
global epsilonEI epsilonIE epsilonII tau
global gammaE thetaVE thetaPhiE I_E
global gammaI thetaVI thetaPhiI I_I
NoConn=-1;

% 0: the Mirollo Strogatz neuron
% 1: the Leaky Integrate-and-Fire neuron
% 2: the SINE neuron
neuronModel = [1 1];

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
thetaV_i = [thetaVE thetaVI];
thetaPhi_i = [thetaPhiE thetaPhiI];
gamma_i=[gammaE gammaI];

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

function ret = phis2vs(phis, I_i, gamma_i)
global neuronModel

N_oscs = size(phis, 2);

if (neuronModel == 0)
elseif (neuronModel == 1)
    ret = NaN(1, N_oscs);
    for i = 1:1:N_oscs
        ret(1, i) = phi2vLIF(phis(1, i), I_i(1, i), gamma_i(1, i));
    end
    
elseif (neuronModel == 2)
end

end

function ret = v2phiss(vs, I_i, gamma_i)
global neuronModel

N_oscs = size(vs, 2);

if (neuronModel == 0)
elseif (neuronModel == 1)
    ret = NaN(1, N_oscs);
    for i = 1:1:N_oscs
        ret(1, i) = v2phiLIF(vs(1, i), I_i(1, i), gamma_i(1, i));
    end
elseif (neuronModel == 2)
end

end

function ret = phi2vLIF(phi, I, gamma)

ret = I/gamma*(1 - exp(-gamma*phi));

end

function ret = v2phiLIF(v, I, gamma)

ret = 1/gamma*log(I/(I - gamma*v));

end