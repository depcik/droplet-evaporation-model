%% Droplet Evaporation Model
% Authors: Amir Madani, Dr. Christopher Depcik
% Version: 2023a
% Date: 08/17/2023


%% Properties

M_F = 0.10021; % [kg / mol] n-Heptane molecular mass
M_N2 = 0.0280134; % [kg / mol] Nitrogen molecular mass
R_u = 8.3144626; % [J / mol K] Universal Gas Constant 
R =  R_u / M_F; % [J / kg K] fuel Gas Constant
rho_l = 684; % [kg / m3] density of the liquid fuel 
L_bar = 36000; % [J / mol] molar vaporization enthalpy
L = L_bar / M_F; % [J / kg] vaporization enthalpy
cp = 2381.4; % [J / kg . K] gas specific heat from analytical solution interpolation
lambda = 0.0504; % [J / s m K] thermal condictivity from analytical solution interpolation

T_amb = 471; % [K] ambient temperature
p_amb = 1e+5; % [Pa] ambient pressure
T_F_ref = 371.5; % [K] liquid fuel boiling temperature from NIST
A = 4.02832;  B = 1268.636; C = -56.199; % Antoine constants
p_F_ref = 100000 * (10 ^ (A - B / (T_F_ref + C))); % [Pa] saturation pressure at the boiling temperature
Y_F_inf = 0; % fuel mass fraction far from droplet
T_inf = T_amb; % temperature far from droplet

d_s = 0.0008; % [m] droplet diameter
r_s = d_s / 2; % [m] droplet radius
r_end = 200 * d_s; % gas domain is 200 times the initial diameter of the droplet
N = 200; % constant number of meshes
del_r = (r_end - r_s) / (N-1); % mesh size
r = r_s : del_r : r_end; % radial direction

%% Initial Guesses

mdot = 5.85e-8; % guess for mass evaporation rate from analytical solution
y = zeros(1, N); y(1:N-1) = 0.001; y(N) = Y_F_inf; % guess for mass fractions to initiate the solver
T = zeros(1, N); T(1:N-1) = 638; T(N) = T_inf; % guess for temperature to initiate the solver
y1 = 1; T1 = 1; % dummy variables to begin the solution
% Errors
Error_m = inf; Error_Req = 1e-8; Error_ReqT = 1e-8; % errors required for each of the unknowns

%% Solution
% Here, mass fraction and temperature of the gas phase are solved until a
% certain error is reached and the solution is complete once the error for
% mass evporation rates reaches the treshold. 

while Error_m > Error_Req
    Z = mdot * cp / (4 * pi * lambda); % a variabel defined to shorten the code
    y = solveSpeciesEquation(N, r, del_r, Z, y, y1, Error_Req); % calling mass fraction function solution
    T = solveTemperatureEquation(N, r, del_r, Z, T, T1, cp, L, T_amb, Error_ReqT); % calling temperature function solution
    [mdot , Error_m] = solveMassEvaporationRate(r_s, R, lambda, cp, T_amb, L, mdot, p_amb, p_F_ref, T_F_ref, M_F, M_N2, y); % calling mass evaporation function solution
end

[xx, yy] = solveEvaporationRateConstant(cp, T,T_inf, rho_l, y, lambda, L, Y_F_inf); % evaporation rate results
[xxx, yyy] = calculateAnalyticalSolution(Y_F_inf, L, cp,lambda, T_inf, p_amb, Error_Req, M_F, M_N2, p_F_ref, R, T_F_ref, rho_l); % analytical solution results to compare with numerical model

