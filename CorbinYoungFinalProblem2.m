
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% PROBLEM 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the parameters for m, b, and k
m = 0.5; % mass in kg
b = 1.5; % damping coefficient in N*s/m
k = 10.0; % spring constant in N/m

% Initial conditions
x0 = 0.1; % initial position in m
v0 = 0.05; % initial velocity in m/s
initial_conditions = [x0, v0];

% Time span
tspan = [0 5]; % from 0 to 5 seconds

% Define the ODE as a function handle
odefun = @(t, x) [x(2); -(b/m)*x(2) - (k/m)*x(1)];

% Use ode45 to solve the ODE
[t, x] = ode45(odefun, tspan, initial_conditions);

% Part B: Plotting
tiledlayout(2,1)

% Position vs Time
nexttile
plot(t, x(:,1))
xlabel('Time (s)')
ylabel('Position (m)')
title('Position vs Time')

% Velocity vs Time
nexttile
plot(t, x(:,2))
xlabel('Time (s)')
ylabel('Velocity (m/s)')
title('Velocity vs Time')

% Part C: Position and velocity at 5 seconds
fprintf('The position at 5 seconds is: %.4f m\n', x(end,1));
fprintf('The velocity at 5 seconds is: %.4f m/s\n', x(end,2));

% Part D: Minimum and maximum position values
fprintf('The minimum position value is: %.4f m\n', min(x(:,1)));
fprintf('The maximum position value is: %.4f m\n', max(x(:,1)));



function mf = finalMass(mi, deltaV, iSP)
    g0 = 9.8; % standard gravity in m/s^2
    mf = mi / exp((deltaV * 1000) / (g0 * iSP)); % Rocket equation, converting km/s to m/s
end

function getDvDisplay(e, rp, mu, v_inf)
    dV2 = getDv(e, rp, mu, v_inf); % Call the getDv function to calculate delta V
    fprintf("Part B: To capture this orbit, the cost, delta V is: %f km/s\n", dV2);
end

function dV2 = getDv(e, rp, mu, v_inf)
    a = rp ./ (1 - e); % Semi-major axis for each e and rp
    vp = sqrt(mu .* (2 ./ rp - 1 ./ a)); % Periapsis velocity for each e and rp
    vc = sqrt(2 .* (v_inf.^2 ./ 2 + mu ./ rp)); % Circular orbit velocity for each rp
    dV2 = vc - vp; % Total delta V for each e and rp, element-wise subtraction
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%    EXTRA CREDIT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% I believe the most interesting assignment I was able to work on 
% was the Map Searching algorithm, for UAV, this lab was very challenging
% for me to understand because I would struggle writing for loops and
% conditions.
