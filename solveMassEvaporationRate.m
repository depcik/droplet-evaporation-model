% Function to solve the mass evaporation rate

function [mdot , Error_m] = solveMassEvaporationRate(r_s, R, lambda, cp, T_amb, L, mdot, p_amb, p_F_ref, T_F_ref, M_F, M_N2, y);

    p_F_s = p_amb * (y(1) / M_F) / ((y(1) / M_F) + ((1-y(1)) / M_N2));
    T_s = -(T_F_ref * L) / (T_F_ref * R * log(p_F_s / p_F_ref) - L);
    mdot1 = 4 * pi * r_s * lambda / cp * log((cp * (T_amb - T_s) + L) / L);
    Error_m = abs(mdot1 - mdot);
    mdot = mdot1;
end
