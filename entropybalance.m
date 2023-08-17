% Each of the terms in the entropy balance equation are solved separately
% in this code, based on the results from the solved model. This file
% should be run after a solution is found.

y_N2 = 1 - y; % nitrogen mass fraction
X_N2 = y_N2./M_N2 ./ (y_N2./M_N2+y./M_F); % nitrogen mole fraction
X = 1 - X_N2; % fuel mole fraction
m_bar = M_N2 .* X_N2 + M_F .* X; % average molecular mass of the mixture 

% Advection of entropy

I1 = zeros (1,N);
for i=1:N-1
    I1(i) = mdot*cp/4/pi/r(i)^2/T(i) * (T(i+1)-T(i))/del_r;  
end
I1(1) = mdot*cp/4/pi/r(1)^2/T(1) * (T(2)-T(1))/del_r; 
I2 = zeros (1,N);
for i=2:N-1
    I2(i) = lambda / r(i)^2 *(2*r(i)/T(i)*(T(i+1)-T(i))/del_r ...
        + r(i)^2*(T(i+1)-T(i))/del_r*(1/T(i+1)-1/T(i))/del_r ...
        + r(i)^2/T(i)*(T(i+1)-2*T(i)+T(i-1))/del_r^2);
end

% Conduction of entropy

I2(1) = lambda / r(1)^2 *(2*r(1)/T(1)*(T(2)-T(1))/del_r ...
        + r(1)^2*(T(2)-T(1))/del_r*(1/T(2)-1/T(1))/del_r ...
        + r(1)^2/T(1)*(T(3)-2*T(2)+T(1))/del_r^2);

% Entropy generated due to conduction

I4 = zeros (1,N);
for i=1:N-1
    I4(i) = lambda / T(i)^2 * ((T(i+1)-T(i))/del_r)^2;
end

% Entropy generated due to mass transfer

for i=1:N-1
I5(i) = lambda / cp * R_u / m_bar(i) * ((((X(i+1)-X(i))/del_r))*(((y(i+1)-y(i))/del_r)) ...
    / y(i)+(((X_N2(i+1)-X_N2(i))/del_r))*((y_N2(i+1)-y_N2(i))/del_r) / y_N2(i));
end

% Diffusion of entropy

I3 = I1-I2-I5-I4;
I_gen = I4+I5;
