% Analytical Solution Function 

function [xxx, yyy] = calculateAnalyticalSolution(Y_F_inf, L, cp_an,lambda_an, T_inf, p_amb, Error_Req, M_F, M_N2, p_F_ref, R, T_F_ref, rho_l)
    T_s_An = 350; % [K] guess for temperature at the droplet surface
    err = 1; % defining error to initiate the solver
    
    % Here, mass fraction at the droplet surface, temperature at the
    % droplet surface and the temperature given by Clausius-Clapeyron
    % equation are solved together until a certain error is reached. Before
    % that, we need to find the specific heat and thermal conductivity of
    % the mixture based with respect to the droplet surface mass fractions 
    % and ambient temperatures. 

    while err > Error_Req
        Y_F_s_An = (L * Y_F_inf + cp_an * (T_inf - T_s_An)) / (L + cp_an * (T_inf - T_s_An)) ; 
        p_F_s_An = p_amb * (Y_F_s_An / M_F) / ((Y_F_s_An / M_F) + ((1-Y_F_s_An) / M_N2));
        T_s_CC_An = -(T_F_ref * L) / (T_F_ref * R * log (p_F_s_An / p_F_ref) - L);
        err = abs (T_s_CC_An - T_s_An);
        T_s_An = T_s_CC_An;
    end

    B_T_An = cp_an * (T_inf - T_s_An) / L; % Spalding heat trasnfer number
    B_M_An = (Y_F_inf - Y_F_s_An) / (Y_F_s_An - 1); % Spalding mass trasnfer number

    K_An = ((8 * lambda_an / cp_an / rho_l) * log (1 + B_T_An)); % [m^2 s^-1] evaporation rate
    K_mm_An = K_An * 1000000; % [mm^2 s^-1] evaporation rate
 
    % Here, evaporation rate is calculated against distance from the
    % droplet surface
    for ii = 1 : 61
        xxx (ii) = (ii - 1) * 0.1;
        yyy (ii) = 1 - K_mm_An * xxx (ii) ;
    end
end