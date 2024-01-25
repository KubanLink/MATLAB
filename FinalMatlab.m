%Part A
mi = 7000; % initial mass kg
iSP = 300; % burn out time seconds
dV1 = 0.5; %km/s
G = 6.67*10^-11; 
g = 9.8; %gravity
finalMass(mi,dV1,iSP);


%PART B

rMars = 3396.19; %radius of mars in kilometers
alt_P = 300; %km
rP = rMars + alt_P; %Argument of Periapsis
e = 0.95; %eccentricity of problem
a = rP/(1-e); %Arguement of Apoapsis
eEllipse = 0.99; %eccentricity for highly eliptical
eCircle = 0; %for circular orbit
mu = 42828.31; %gravitivational parameter km^3/s^2
v_inf = 2; %km/s
vp = sqrt(mu*(2/rP - 1/a)); %periapsis velocity
vc = sqrt(2*(v_inf^2/2 + mu/rP)); %velocity of capture orbit
getDv(e,rP,mu,v_inf);

%Part C
[ecc,alt_P] = meshgrid(0.00:0.01:0.99,100:50:1000);
valueDV = getDv(ecc,rP,mu,v_inf);
contourf(ecc,alt_P,valueDV)






function finalMass(mi, deltaV, iSP)
g = 9.8; %gravity
mf = (mi)/exp(deltaV/(g*iSP));
fprintf("PART A: The value final mass of the propellant is: %f\n",mf);
end


function getDv(e,rp,mu,v_inf)
a = rp./(1-e); %Arguement of Apoapsis
vp = sqrt(mu*(2./rp - 1./a)); %periapsis velocity
vc = sqrt(2*(v_inf^2/2 + mu./rp)); %velocity of capture orbit
dV2 = vc-vp;
fprintf("Part B: To capture this orbit, the cost, delta V is: %f\n",dV2)
end
