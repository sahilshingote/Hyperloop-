%parameter - defining

syms v; %Linear Speed
syms is; %Inductor current space vector in the inductor reference frame
syms ir_dash; % Induced-part current space vector in the inductor reference frame
syms psi_r_dash; %Induced-part flux space vector in the inductor reference frame
syms psi_s; %Inductor flux space vector in the inductor reference frame
syms us; %. Inductor voltage space vector in the inductor reference frame
syms L_m_hat;
syms R_r_hat;
syms L_s_hat;
syms R_S_hat;
syms T_r_hat; %induced part time machine constant
syms Q; %end effect factor
syms sigma_hat; %equivalent global leakage factor
syms t; %Time

tau_m = 0.3; %length of the inductor
tau_p = 1.075; %polar pitch
v = 6.85; %linear speed
wr = 60; %angular rotor speed
p = 3; %number of pole pairs
L_s = 637.6e-3; %inductor resistance
L_r = 757.8; %induced-part resistance
L_m = 517.5; %three-phase magnetizing inductance
R_r = 32.57 ; %d induced-part resistance -- ohms
Rs = 11; % inductor resistance
L_sigma_r = 1; %induced-part leakage inductance
L_sigma_s = 1; %inductor leakage resistance
I = eye(2);
C = [ I  zeros(2, 2)];
I = [1 0;0 1];
J = [0 -1;1 0];


Q = (tau_m*R_r)/((L_m + L_sigma_r) + v) ;

out = f_q(Q);

L_m_hat = L_m*(1-out);
R_r_hat = R_r*out;



psi_s = (L_sigma_s + L_m*(1-out))*is + L_m * (1-out)*ir_dash;
psi_r_dash = L_m*(1-out)*is + (L_sigma_r + L_m*(1-out)*ir_dash);

us = Rs*is + R_r*out * (is + ir_dash) + diff(psi_s, t);

sigma_hat = 1 - ((L_m * (1 - out))^2)/( (L_sigma_r + L_m * (1-out))*(L_sigma_s + L_m * (1-out)));

L_s_hat = L_sigma_s + L_m * (1-out);
L_r_hat = L_sigma_r + L_m * (1-out);

T_r_hat = (L_sigma_r + L_m * (1-out))/(R_r*(1+out));

%A-matrix

A11 = -(1/(sigma_hat*L_s_hat)*(Rs + R_r_hat*(1 - (L_m_hat/L_r_hat)) + (L_m_hat/L_r_hat)*((L_m_hat/T_r_hat)-R_r_hat)))*I;


a  = (1/T_r_hat) + (R_r_hat);
b = (p*pi/tau_p)*v;

A12 = (L_m_hat/(sigma_hat*L_s_hat*L_r_hat))*(a*I - b*J);
A21 = ((L_m_hat/T_r_hat) - R_r_hat)*I;
A22 = -((1/T_r_hat) - b*J);

B1 = (1/sigma_hat)*(1/L_s_hat)*I;

A  = [A11 A12; A21 A22];

B = [B1; zeros(2, 2)];


D = [zeros(1, 2); zeros(1, 2 )];

f = linspace(10,250,800);
w = 2*pi*f;
sys = ss(A, B ,C, D);
Fresp = freqresp(sys,w);
Fresp = Fresp(4:4:3200);
%plot(f,transpose(abs(Fresp)))
Fe = 3/2*pi*p/tau_p*L_m_hat/L_r_hat*abs(psi_r_dash)*is;
%Fresp = lsim(sys,u);

ode1 = diff(is, t) == A11*is + A12*psi_r_dash + B1*us;
ode2 = diff(psi_r_dash, t) == A21*is + A22*psi_r_dash;

odes = [ode1; ode2];

%[isSol, psi_r_dashSol] = dsolve(odes);

%isSol = simplify(isSol, t);
%psi_r_dashSol = simplify(psi_r_dashSol, t);



function out = f_q(Q)

out = (1-exp(-Q))/Q;

end