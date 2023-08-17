% Function to solve the species equation
function y = solveSpeciesEquation(N, r, del_r, Z, y, y1, Error_Req)
    Errory = inf; % defining error to initiate the solver

    % Here, species equation is solved until a certain error is reached.
   
    while Errory > Error_Req
        Errory = 0;
        for i = 2 : N-1
            yold(i) = y(i);
            y(i) = (((y(i+1)+y(i-1))*r(i)^2) + (2*y(i+1)*del_r*r(i)) - y(i+1)*del_r*Z) / (2*r(i)^2 + 2*del_r*r(i) - del_r*Z);
            Errory = max(Errory,abs(yold(i)-y(i)));
        end
        y(1) = ((2*y(2)*r(1)^2) + (2*y(2)*del_r*r(1)) + ((2-y(2))*del_r*Z )) / (2*r(1)^2 + 2*del_r*r(1) + del_r*Z);
        Errory = max(Errory, abs(y(1) - y1));
        y1 = y(1);
    end
end