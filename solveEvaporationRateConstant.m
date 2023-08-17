% Function to solve the evaporation rate constant
function [xx, yy] = solveEvaporationRateConstant(cp, T,T_inf, rho_l, y, lambda, L, Y_F_inf)
    B_T = cp * (T_inf - T(1)) / L; 
    B_M = (Y_F_inf - y(1)) / (y(1) - 1);

    K = ((8 * lambda / cp / rho_l) * log (1 + B_T));  % evaporation rate constant 
    K_mm = K * 1000000;

    for i = 1 : 61 % domain setup for plot
        xx (i) = (i - 1) * 0.1;
        yy (i) = 1 - K_mm * xx (i) ;
    end
end