function [posEst,linVelEst,oriEst,windEst,driftEst,...
    posVar,linVelVar,oriVar,windVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,windEst,driftEst,...
%    posVar,linVelVar,oriVar,windVar,driftVar,estState] =
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time t_k, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   windEst         wind direction estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   windVar         variance of wind direction estimate(time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!

    % initial state mean
    posEst = [0 0]; % 1x2 matrix
    linVelEst = [0 0]; % 1x2 matrix
    oriEst = 0; % 1x1 matrix
    windEst = 0; % 1x1 matrix
    driftEst = 0; % 1x1 matrix

    % initial state variance
    posVar = [estConst.StartRadiusBound^2/4 estConst.StartRadiusBound^2/4];
    % 1x2 matrix %uniform disk distri. variance, see draft
    linVelVar = [0 0]; % 1x2 matrix
    oriVar = estConst.RotationStartBound^2/3; % 1x1 matrix
    windVar = estConst.WindAngleStartBound^2/3; % 1x1 matrix
    driftVar = 0; % 1x1 matrix

    %initial random values
    r_ini = random('Uniform',0,estConst.StartRadiusBound);
    angle_ini = random('Uniform',0,2*pi);
    px_ini = sqrt(r_ini)*cos(angle_ini);
    py_ini = sqrt(r_ini)*sin(angle_ini);
    phi_ini = random('Uniform',-estConst.RotationStartBound,-estConst.RotationStartBound);
    rou_ini = random('Uniform',-estConst.WindAngleStartBound,-estConst.WindAngleStartBound);

    estState = struct;
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar,linVelVar,oriVar,windVar,driftVar]);
    % estimator state
    estState.xm = [px_ini, py_ini, 0, 0, phi_ini, rou_ini, 0]';%follow the orders above
    %disp(estState.xm);
    % time of last update
    estState.tm = tm;
end

%% Estimator iteration.
% get time since last estimator update
syms px py sx sy phi rou b

syms ut ur
syms vd vr vp vb

syms za zb zc zg zn
syms wa wb wc wg wn

%% state vectors input
if tm==0
    [ut, ur] = [0,0];
else
    [ut, ur] = deal(actuate(1),actuate(2));
end
%[za, zb, zc, zg, zn] = deal(sense(1),sense(2),sense(3),sense(4),sense(5));
V = [vd vr vp vb];
X = [px py sx sy phi rou b];
U = [ut ur];
Z = [za zb zc zg zn];
W = [wa wb wc wg wn];

%% Process equations
Dpx = sx;
Dpy = sy;
Dsx = cos(phi)*(tanh(ut)-const.dragCoefficientHydr*(sx^2+sy^2)*(1+vd))...
    - const.dragCoefficientAir*(sx-const.windVel*cos(rou))...
    *sqrt((sx-const.windVel*cos(rou))^2+(sy-const.windVel*sin(rou))^2);
Dsy = sin(phi)*(tanh(ut)-const.dragCoefficientHydr*(sx^2+sy^2)*(1+vd))...
    - const.dragCoefficientAir*(sy-const.windVel*sin(rou))...
    *sqrt((sx-const.windVel*cos(rou))^2+(sy-const.windVel*sin(rou))^2);
Dphi = const.rudderCoefficient*ur*(1+vr);
Drou = vp;
Db = vb;

%% measurement equations
za = sqrt((px-const.pos_radioA(1))^2+(py-const.pos_radioA(2))^2)+wa;
zb = sqrt((px-const.pos_radioB(1))^2+(py-const.pos_radioB(2))^2)+wb;
zc = sqrt((px-const.pos_radioC(1))^2+(py-const.pos_radioC(2))^2)+wc;
zg = phi+b+wg;
zn = phi+wn;

%% jacobian matrices
At = jacobian([Dpx,Dpy,Dsx,Dsy,Dphi,Drou,Db], X);
% JacoU = jacobian([Dpx,Dpy,Dsx,Dsy,Dphi,Drou,Db], U)
Lt = jacobian([Dpx,Dpy,Dsx,Dsy,Dphi,Drou,Db], V);
Hk = jacobian([za zb zc zg zn], X);
Mk = jacobian([za zb zc zg zn], W);
% measurement noise variance
R = diag([estConst.DistNoiseA,estConst.DistNoiseB,estConst.DistNoiseC,estConst.GyroNoise,estConst.CompassNoise]);
% process noise variance
Qsys= diag([const.DragNoise,const.RudderNoise,const.WindAngleNoise,const.GyroDriftNoise]);

%dt = ...
% update measurement update time
estState.tm = tm; 

%% prior update
syms px(t) py(t) sx(t) sy(t) phi(t) rou(t) b(t)
eqs = [diff(px(t),t) == sx,...
       diff(py(t),t) == sy,...
       diff(sx(t),t) == cos(phi)*(tanh(ut)-const.dragCoefficientHydr*(sx^2+sy^2)*(1+vd))...
       - const.dragCoefficientAir*(sx-const.windVel*cos(rou))...
       *sqrt((sx-const.windVel*cos(rou))^2+(sy-const.windVel*sin(rou))^2),...
       diff(sx,t) == sin(phi)*(tanh(ut)-const.dragCoefficientHydr*(sx^2+sy^2)*(1+vd))...
       - const.dragCoefficientAir*(sy-const.windVel*sin(rou))...
       *sqrt((sx-const.windVel*cos(rou))^2+(sy-const.windVel*sin(rou))^2),...
       diff(phi(t),t) == const.rudderCoefficient*ur*(1+vr),...
       diff(rou(t),t) == vp,...
       diff(b(t),t) == vb];
vars = [px(t) py(t) sx(t) sy(t) phi(t) rou(t) b(t)];

[M,F] = massMatrixForm(eqs,vars);
M = odeFunction(M,vars);
F = odeFunction(F,vars,vd,vr,vp,vb);
vd=0;vr=0;vp=0;vb=0;
F = @(t,Y) F(t,Y,vd,vr,vp,vb);
[t_list,xp_list] = ode45(F,[tm*0.1 (tm+1)*0.1],estState.xm);

xp=xp_list(end,:)';

% dot Psv = At*Psv+Psv*At.'+Lt*Qsys*Lt.';
% assign values to system states
[px, py, sx, sy, phi, rou, b] = deal(estState.xm(1),estState.xm(2),estState.xm(3),estState.xm(4),estState.xm(5),estState.xm(6),estState.xm(7));
% obtain double from symbolic matrix
Akt = double(subs(At));
Qkt = double(subs(Lt*Qsys*Lt.'));
% feed to ode45
[Tlist, Plist] = ode45(@(t,Ps)VarianceDE(t, Ps, Akt, Qkt), [tm*0.1 (tm+1)*0.1],estState.Pm);
Pp=reshape(Plist(end,:),[7,7]);

%% measurement update
% assign xp(k) to symbolic
[px, py, sx, sy, phi, rou, b] = deal(xp(1),xp(2),xp(3),xp(4),xp(5),xp(6),xp(7));
Hk = double(subs(Hk));
Mk = double(subs(Mk));
wa=0;wb=0;wc=0;wg=0;wn=0;
% our nonlinear prediction:hk(xp,0)
hkxp = [double(subs(za)),double(subs(zb)),double(subs(zc)),double(subs(zg)),double(subs(zn))]';
% Kalman gain:
Kk = Pp*Hk'/(Hk*Pp*Hk' + Mk*R*Mk');

% give more trust to zc if possible
if isinf(sense(3))
    sense(3)=double(subs(zc));
end

%update mean
%no measurement at initial step, use xp
if tm==0
    estState.xm = xp;
else
    estState.xm = xp + Kk*(sense' - hkxp);
end

estState.Pm = (eye(7)-Kk*Hk)*Pp;
P_diag= diag(estState.Pm);

%% Get resulting estimates and variances

% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(3:4);
oriEst = estState.xm(5);
windEst = estState.xm(6);
driftEst = estState.xm(7);
% 
posVar = P_diag(1:2);    % 1x2 matrix
linVelVar = P_diag(3:4); % 1x2 matrix
oriVar = P_diag(5); 
windVar = P_diag(6);
driftVar = P_diag(7); 

end

%% ode function handle
function dPdt = VarianceDE(t, X, A, Q)
X = reshape(X, size(A)); %Convert from "n^2"-by-1 to "n"-by-"n"
dPdt = A*X + X*A.' + Q; %Determine derivative
dPdt = dPdt(:); %Convert from "n"-by-"n" to "n^2"-by-1
end
