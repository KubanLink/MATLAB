%
%   ARO 3220
%   Control and Feedback Systems
%   Dr. Marco Maggia
%   12/15/23
%
%   N. Della Rocco, T. Roaf, C. Young
%
%                           Final Project: 
%   Longitudinal Dynamics Controller Design for a Boeing 747-100
%   

K=1
%% User Defined Variables

%Physical Data:
g = 9.807; %Gravitational Acceleration at Sea Level (meters/s^2)
rho_ssl = 1.225; %Air Density at SSL (kg/meters^3)
h_s = 8500; %Atmospheric Scale Height of Earth (meters)

%For a Boeing 747-100

m = 636636 * 0.453592; %lbs to kg
b_wing = 195.7 * 0.3048; %ft to meters
c_wing = 27.3 * 0.3048; %ft to meters
s_wing = 5500 * 0.092903; %ft^2 to meters^2

Ixx=18.2*(10^6) * 1.35581795; %slug-ft^2 to kg*meters^2
Iyy=33.1*(10^6) * 1.35581795; %slug-ft^2 to kg*meters^2
Izz=49.7*(10^6) * 1.35581795; %slug-ft^2 to kg*meters^2

Ixy = 0;
Iyz = 0;
Ixz = -1.56*(10^6) * 1.35581795; %slug-ft^2 to kg*meters^2

T_max = 827.2*(10^3); %N %kg-m/s^2

%Drag Coefficients:
C_D0 = 0.0275; C_Dalpha = 0.425; C_DdelE = 0.000;

%Lift Coefficients:
C_L0 = 0.2050; C_Lalpha = 4.920; C_LdelE = 0.367; C_LQ = 6.000;

%Pitching Moment Coefficients:
C_m0 = 0.1660; C_malpha = -1.033; C_mdelE = -1.450; C_mQ = -24.00;

%Side Force Coefficients:
C_Ybeta = -0.880; C_YdelA = 0.000; C_YdelR = 0.1157;

%Rolling Moment Coefficients:
C_lbeta = -0.277; C_ldelA = 0.0137; C_ldelR = 0.007; C_lP = -0.334;
C_lR = -0.3;

%Yawing Moment Coefficients:
C_nbeta = 0.195; C_ndelA = 0.0002; C_ndelR = -0.1256; C_nP = -0.0415;
C_nR = -0.327;


%%


%Given:
h = 38000 *0.3048; %feet to meters
V_inf_0 = 450 * 0.514444; %knots to meters/second


% Part One: Steady Level Flight (Quasy-Steady-State or Trim Condition) 
% 
% Find the angle of attack (in deg), the elevator deflection
% (in deg) and the throttle level (in %) necessary for a B747-100 to
% maintain a constant altitude of 38,000 ft at a constant speed of 450
% knots. 

rho_inf_0 = rho_ssl * exp(-h / h_s); % Free Stream Density (kg/m^3)
q_inf_0 = 0.5 * rho_inf_0*(V_inf_0^2);   % Free Stream Dynamic Pressure (Pascals, kg/(m*s^2))


eqn = @(alpha) ((C_D0 - (C_m0 / C_mdelE) * C_DdelE) + ...
    (C_Dalpha - (C_malpha / C_mdelE) * C_DdelE) * alpha) * sin(alpha) + ...
    ((C_L0 - (C_m0 / C_mdelE) * C_LdelE - ((m * g) / (q_inf_0*s_wing))) + ...
    (C_Lalpha - (C_malpha / C_mdelE) * C_LdelE) * alpha) * cos(alpha);

alpha_0 = fsolve(eqn, 0)

alpha_deg = (alpha_0)*(180/pi);

delE_0 = (-1 / C_mdelE) * (C_m0 + C_malpha * alpha_0)
delE_deg = delE_0*(180/pi);

D_0= (C_D0 + C_Dalpha * alpha_0 + C_DdelE * delE_0) * q_inf_0 * s_wing;

L_0 = (C_L0 + C_Lalpha * alpha_0 + C_LdelE * delE_0) * q_inf_0 * s_wing;

pi_T_0 = (1 / T_max) * (D_0 * cos(alpha_0) + ((m * g) - L_0) * sin(alpha_0))

T_0 = pi_T_0 * T_max

par=[m,c_wing, s_wing, T_max, Iyy, ...
    g, h_s, rho_ssl, C_D0, C_Dalpha, C_DdelE, ...
    C_L0, C_Lalpha, C_LQ, C_LdelE, C_malpha, C_mQ, C_mdelE, ...
    D_0, L_0, T_0, V_inf_0, q_inf_0, alpha_0];


A(1,:) = [0,0,0,0,1];
A(2,:) = [V_inf_0, 0, 0, -V_inf_0, 0];
A(3,:) = [ -g, D_0/(m*h_s), (-2*D_0)/(m*V_inf_0), g-((T_0*sin(alpha_0))/m)-((q_inf_0*s_wing*C_Dalpha)/m),0];
A(4,:) = [0, L_0/(m*V_inf_0*h_s), (-g/(V_inf_0^2))-(T_0*sin(alpha_0)/(m*V_inf_0^2))-(L_0/(m*V_inf_0^2)), ...
    (-T_0*cos(alpha_0)/(m*V_inf_0))-(q_inf_0*s_wing*C_Lalpha)/(m*V_inf_0), 1 - ((q_inf_0*s_wing*C_LQ)/(m*V_inf_0)) ];
A(5,:) = [0, 0, 0, (q_inf_0*s_wing*c_wing*C_malpha)/Iyy, (q_inf_0*s_wing*c_wing*C_mQ)/Iyy];

eig_A = eig(A)
eig_realA = real(eig_A)

wn = abs(eig_A)

zeta = -eig_realA./wn

Dh0=2000*0.3048; %2000ft to meters
DVinf0 = 25* 0.514444;

%%

%Dynamics Block:

% dtheta = 0;
% dh = h;
% dVinf = V_inf_0;
% dalpha = alpha_0;
% dQ = 0;
% dpi_T = pi_T_0;
% ddelE = delE_0;
% 
% 
% m = par(1); c_wing = par(2); s_wing=par(3); T_max=par(4);
% Iyy = par(5); g = par(6); h_s = par(7); rho_ssl = par(8);
% C_D0 = par(9); C_Dalpha = par(10); C_DdelE = par(11);
% C_L0 = par(12); C_Lalpha = par(13); C_LQ = par(14);
% C_LdelE = par(15); C_malpha = par(16); C_mQ = par(17);
% C_mdelE = par(18); D_0 = par(19); L_0 = par(20);
% T_0 = par(21); V_inf_0 = par(22); q_inf_0 = par(23);
% alpha_0 = par(24);
% 
% 
% x = [dtheta, dh, dVinf, dalpha, dQ]';
% u = [dpi_T, ddelE]';
% 
% A = zeros(5,5);
% B = zeros(5,2);
% 
% A(1,:) = [0,0,0,0,1];
% A(2,:) = [V_inf_0, 0, 0, -V_inf_0, 0];
% A(3,:) = [ -g, D_0/(m*h_s), (-2*D_0)/(m*V_inf_0), g-(T_0*sin(alpha_0)/m)-(q_inf_0*s_wing*C_Dalpha)/m,0];
% A(4,:) = [0, L_0/(m*V_inf_0*h_s), (-g+(T_0*sin(alpha_0)/m)-(L_0/m))/V_inf_0^2, ...
%     (-T_0*cos(alpha_0)/(m*V_inf_0))-(q_inf_0*s_wing*C_Lalpha)/(m*V_inf_0), 1 - ((q_inf_0*s_wing*C_LQ)/(m*V_inf_0)) ];
% A(5,:) = [0, 0, 0, (q_inf_0*s_wing*c_wing*C_malpha)/(Iyy), (q_inf_0*s_wing*c_wing)/(Iyy*C_mQ)];
% 
% 
% B(1,:) = [0,0];
% B(2,:) = [0,0];
% B(3,:) = [(T_max*cos(alpha_0))/(m*V_inf_0), (-q_inf_0*s_wing*C_LdelE)/(m*V_inf_0)];
% B(4,:) = [(-T_max*cos(alpha_0))/(m*V_inf_0), (-q_inf_0*s_wing*C_DdelE)/(m*V_inf_0)];
% B(5,:) = [0, (q_inf_0*s_wing*c_wing*C_mdelE)/Iyy];
% 
% x_dot = A*x + B*u;
% 
% 
% 
% x_dot_T = transpose(x_dot)

    
