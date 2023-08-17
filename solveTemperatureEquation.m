% Function to solve the energy equation
function T = solveTemperatureEquation(N, r, del_r, Z, T, T1, cp, L, T_amb, Error_ReqT)
    ErrorT = inf; % defining error to initiate the solver

    % Here, energy equation is solved until a certain error is reached.

    while ErrorT > Error_ReqT
        ErrorT = 0;
        for i = 2 : N-1
            Told(i) = T(i);
            T(i) = (((T(i+1)+T(i-1))*r(i)^2) + (2*T(i+1)*del_r*r(i)) - T(i+1)*del_r*Z) / (2*r(i)^2 + 2*del_r*r(i) - del_r*Z);
            ErrorT = max(ErrorT, abs(Told(i)-T(i)));
        end
        T(1) = ((2*T(2)*cp*r(1)^2) + (2*T(2)*cp*del_r*r(1)) - ((2*del_r*L+T(2)*del_r*cp )*Z)) / (2*cp*r(1)^2 + 2*cp*del_r*r(1) - cp*del_r*Z);
        ErrorT = max(ErrorT, abs(T(1) - T1));
        T1 = T(1);
    end
end