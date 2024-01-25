clc;clear;close all
s = tf('s'); %Defines general TF
%% %%%%%%%%%%%  USER INPUTS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Uncompensated OL System TF
G = 1/(s*(s^2*(s+1))); 

% Performance Requirements
PO_req = 35;  % Percent Overshoot
Ts_req = 4;   % Settling Time
Tp_req = Inf; % Peak Time

% Desired Dominant Poles     --> 1st "Arbitrary" Design Step
Re_DP = -1; % Real Part
Im_DP =  2; % Imaginary Part

% Compensator ZERO Placement --> 2nd "Arbitrary" Design Step
z_LEAD = abs(Re_DP);

%% %%%%%%%%%% NO NEED TO EDIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wn   = norm([Re_DP,Im_DP]);
zeta = abs(Re_DP)/wn;

fprintf('Damping Ratio: zeta = %.3f\n',zeta)
fprintf('Natural Freq.:   wn = %.3f rad/s\n\n',wn)

s_DP1 = Re_DP + Im_DP*1i;
s_DP2 = Re_DP - Im_DP*1i;
fprintf('Dominant Poles:  s_DP = %.2f +/- %.2fi\n\n',Re_DP,Im_DP)

% CL Sys Performance Estimates (based on DP)
PO_est = 100*exp(-pi*zeta/(sqrt(1-zeta^2)));
Ts_est = 4/(zeta*wn);
Tp_est = pi/(wn*sqrt(1-zeta^2));

if PO_est<=PO_req
    fprintf('Estimate PO = %.1f <= %.1f %%\n',PO_est,PO_req)
else
    fprintf('Estimate PO = %.1f > %.1f %%\n Change Dominant Poles!\n',PO_est,PO_req)
end
if Ts_est<=Ts_req
    fprintf('Estimate Ts = %.2f <= %.2f sec\n',Ts_est,Ts_req)
else
    fprintf('Estimate Ts = %.2f > %.2f sec\n Change Dominant Poles!\n',Ts_est,Ts_req)
end
if Tp_est<=Tp_req
    fprintf('Estimate Tp = %.1f <= %.1f sec\n\n',Tp_est,Tp_req)
else
    fprintf('Estimate Tp = %.1f > %.1f sec\n Change Dominant Poles!\n\n',Tp_est,Tp_req)
end

%% Compensator ZERO Placement (see earlier)
fprintf('Compensator zero at  s = %.3f\n',-z_LEAD)

%% Compensator POLE Placement
% 1 + K*(s+z_LEAD)/(s+p_LEAD)/s^2 = 0
% K*(s+z_LEAD)/(s+p_LEAD)/s^2 = -1
% K*(Re_DP+z_LEAD+Im_DP*1i)/(Re_DP+p_LEAD+Im_DP*1i)/(Re_DP+Im_DP*1i)^2 = -1
% atan2d(Im_DP,Re_DP+z_LEAD) - atan2d(Im_DP,Re_DP+p_LEAD) - 2*atand2(Im_DP,Re_DP) = 180
% -atan2d(Im_DP,Re_DP+p_LEAD) = 180 + 2*atand2(Im_DP,Re_DP) - atan2d(Im_DP,Re_DP+z_LEAD)
% atan2d(Im_DP,Re_DP+p_LEAD) = -180 - 2*atand2(Im_DP,Re_DP) + atan2d(Im_DP,Re_DP+z_LEAD)
% Im_DP/(Re_DP+p_LEAD) = tand(-180 - 2*atan2d(Im_DP,Re_DP) + atan2d(Im_DP,Re_DP+z_LEAD))

a = tand(-180 - 2*atan2d(Im_DP,Re_DP) + atan2d(Im_DP,Re_DP+z_LEAD));

% a = Im_DP/(Re_DP+p_LEAD)

p_LEAD = Im_DP/a - Re_DP;
fprintf('Compensator pole at  s = %.3f\n',-p_LEAD)

alpha = p_LEAD/z_LEAD;
if alpha>1
    fprintf('Phase-LEAD Compensator: alpha = %.2f\n\n',alpha)
else
    fprintf('Review Your Design!\n\n')
    return
end

% Compensator TF
Gc = (s+z_LEAD)/(s+p_LEAD)

%% GAIN Determination
K = norm(s_DP1+p_LEAD)*norm(s_DP1)^2/norm(s_DP1+z_LEAD);
fprintf('GAIN: K = %.2f\n\n',K)

% Loop TF
L = K*Gc*G;

% CL Uncompensated (T_unc) and Compensated (T_com) TFs.
T_unc = minreal(G/(1+G));
T_com = minreal(L/(1+L));

figure(2),hold on
set(gcf,'units','normalized','position',[0.5 0 0.5 0.8])
% Uncompensated System
subplot(2,4,[1,2]),hold on
    rlocus(G)
    l1 = plot(Re_DP, Im_DP,'Marker','square','MarkerEdgeColor','r','MarkerFaceColor','none','MarkerSize',12,'LineStyle','none');
    plot(Re_DP,-Im_DP,'Marker','square','MarkerEdgeColor','r','MarkerFaceColor','none','MarkerSize',12,'LineStyle','none')
    l2 = plot(real(pole(T_unc)),imag(pole(T_unc)),'Marker','^','MarkerEdgeColor','m','MarkerFaceColor','m','MarkerSize',6,'LineStyle','none');
    axis('equal')
    legend([l1,l2],'Desired Dominant Poles','Uncompensated CL Poles')
subplot(2,4,[3,4]),hold on
    step(T_unc,50)

% Compensated System
subplot(2,4,[5,6]),hold on
    rlocus(L)
    l3 = plot(Re_DP, Im_DP,'Marker','square','MarkerEdgeColor','r','MarkerFaceColor','none','MarkerSize',12,'LineStyle','none');
    plot(Re_DP,-Im_DP,'Marker','square','MarkerEdgeColor','r','MarkerFaceColor','none','MarkerSize',12,'LineStyle','none')
    l4 = plot(real(pole(T_com)),imag(pole(T_com)),'Marker','>','MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',6,'LineStyle','none');
    axis('equal')
    legend([l3,l4],'Desired Dominant Poles','Compensated CL Poles')
subplot(2,4,[7,8]),hold on
    step(T_com)