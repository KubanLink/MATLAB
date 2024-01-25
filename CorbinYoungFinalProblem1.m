% Part A
mi = 7000; % initial mass kg
iSP = 300; % burn out time seconds
dV1 = 0.5; % deltaV in km/s
finalMassValue = finalMass(mi, dV1, iSP);
fprintf("The initial final mass for Part A is: %.2f kg\n", finalMassValue);


% Part B
rMars = 3396.19; % radius of Mars in kilometers
alt_P = 300; % periapsis altitude in km
e = 0.95; % eccentricity
mu = 42828.31; % gravitational parameter in km^3/s^2
v_inf = 2; % v infinity in km/s
getDvDisplay(e, alt_P + rMars, mu, v_inf);

% Part C
[ecc, alt_P_grid] = meshgrid(0.00:0.01:0.99, 100:50:1000); % Eccentricity and altitude grid
rP_grid = rMars + alt_P_grid; % Calculate rP for each altitude in alt_P grid
valueDV = getDv(ecc, rP_grid, mu, v_inf); % Calculate delta V for each combination of e and rP

% Part D Calculate final masses (mf) for the contour plot using the function from Part A
final_masses = arrayfun(@(dV) finalMass(mi, dV, iSP), valueDV)

% Plotting the contour plot with labels
figure()
contourf(ecc, alt_P_grid, final_masses, 'ShowText', 'on');
c = colorbar; % Add a colorbar
c.Label.String = 'Spacecraft Final Mass (m_f) (kg)'; % Label the colorbar
c.Label.FontSize = 16; % Set the font size for the colorbar label
xlabel('Eccentricity (e)'); % Label for the x-axis
ylabel('Periapsis Altitude (km)'); % Label for the y-axis
title('Contour Plot of Spacecraft Final Mass'); % Title of the plot


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
