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
    % 1x2 matrix %uniform disk distri. variance, see draft
    posVar = [estConst.StartRadiusBound^2/4 estConst.StartRadiusBound^2/4];
    linVelVar = [0 0]; % 1x2 matrix
    oriVar = estConst.RotationStartBound^2/3; % 1x1 matrix
    windVar = estConst.WindAngleStartBound^2/3; % 1x1 matrix
    driftVar = 0; % 1x1 matrix

    estState = struct;
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar,linVelVar,oriVar,windVar,driftVar]);
    % estimator state
    estState.xm = [posEst, linVelEst, oriEst, windEst, driftEst]';%follow the orders above
    % time of last update
    estState.tm = -0.1;
end

%% Estimator iteration.
% get time since last estimator update
%px=0;py=0;sx=0;sy=0;phi=0;rou=0;b=0;
wa=0;wb=0;wc=0;wg=0;wn=0;
vd=0;%vr=0;vp=0;vb=0;

if tm==0
    ut=0;ur=0; %no input
else
    [ut, ur] = deal(actuate(1),actuate(2));
end

%% jacobian matrices
% measurement noise variance
R = diag([estConst.DistNoiseA,estConst.DistNoiseB,estConst.DistNoiseC,estConst.GyroNoise,estConst.CompassNoise]);
% process noise variance
Qsys= diag([estConst.DragNoise,estConst.RudderNoise,estConst.WindAngleNoise,estConst.GyroDriftNoise]);

% update time
dt = tm - estState.tm;
estState.tm = tm; 

%% prior update
% mean
% offline calculated function handle

ht = @(t,Y)[Y(3);Y(4);-cos(Y(5)).*(-tanh(ut)+Y(3).^2./1.0e+1+Y(4).^2./1.0e+1)+sqrt((cos(Y(6)).*(3.0./4.0)-Y(3)).^2+(sin(Y(6)).*(3.0./4.0)-Y(4)).^2).*(cos(Y(6)).*(9.0./2.0e+2)-Y(3).*(3.0./5.0e+1));-sin(Y(5)).*(-tanh(ut)+Y(3).^2./1.0e+1+Y(4).^2./1.0e+1)+(sin(Y(6)).*(9.0./2.0e+2)-Y(4).*(3.0./5.0e+1)).*sqrt((cos(Y(6)).*(3.0./4.0)-Y(3)).^2+(sin(Y(6)).*(3.0./4.0)-Y(4)).^2);ur.*2.0;0.0;0.0];

[~,xp_list] = ode45(ht,[tm-dt tm],estState.xm);
xp=xp_list(end,:)';

% variance
% dot Psv = At*Psv+Psv*At.'+Lt*Qsys*Lt.';
% assign values to system states
[px, py, sx, sy, phi, rou, b] = deal(estState.xm(1),estState.xm(2),estState.xm(3),estState.xm(4),estState.xm(5),estState.xm(6),estState.xm(7));

At = [0, 0,                                                                                                                                                                                                                      1,                                                                                                                                                                                                                      0,                                                   0,                                                                                                                                                                                                                                                           0, 0;
0, 0,                                                                                                                                                                                                                      0,                                                                                                                                                                                                                      1,                                                   0,                                                                                                                                                                                                                                                           0, 0;
0, 0, - (3*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2))/50 - (sx*cos(phi)*(vd + 1))/5 - ((2*sx - (3*cos(rou))/2)*((3*sx)/50 - (9*cos(rou))/200))/(2*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2)),                                                                    - (sy*cos(phi)*(vd + 1))/5 - (((3*sx)/50 - (9*cos(rou))/200)*(2*sy - (3*sin(rou))/2))/(2*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2)), -sin(phi)*(tanh(ut) - (sx^2/10 + sy^2/10)*(vd + 1)), - (9*sin(rou)*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2))/200 - (((3*sx)/50 - (9*cos(rou))/200)*((3*sin(rou)*(sx - (3*cos(rou))/4))/2 - (3*cos(rou)*(sy - (3*sin(rou))/4))/2))/(2*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2)), 0;
0, 0,                                                                    - (sx*sin(phi)*(vd + 1))/5 - ((2*sx - (3*cos(rou))/2)*((3*sy)/50 - (9*sin(rou))/200))/(2*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2)), - (3*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2))/50 - (sy*sin(phi)*(vd + 1))/5 - ((2*sy - (3*sin(rou))/2)*((3*sy)/50 - (9*sin(rou))/200))/(2*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2)),  cos(phi)*(tanh(ut) - (sx^2/10 + sy^2/10)*(vd + 1)),   (9*cos(rou)*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2))/200 - (((3*sy)/50 - (9*sin(rou))/200)*((3*sin(rou)*(sx - (3*cos(rou))/4))/2 - (3*cos(rou)*(sy - (3*sin(rou))/4))/2))/(2*((sx - (3*cos(rou))/4)^2 + (sy - (3*sin(rou))/4)^2)^(1/2)), 0;
0, 0,                                                                                                                                                                                                                      0,                                                                                                                                                                                                                      0,                                                   0,                                                                                                                                                                                                                                                           0, 0;
0, 0,                                                                                                                                                                                                                      0,                                                                                                                                                                                                                      0,                                                   0,                                                                                                                                                                                                                                                           0, 0;
0, 0,                                                                                                                                                                                                                      0,                                                                                                                                                                                                                      0,                                                   0,                                                                                                                                                                                                                                                           0, 0];

Lt =[                       0,    0, 0, 0;
                            0,    0, 0, 0;
-cos(phi)*(sx^2/10 + sy^2/10),    0, 0, 0;
-sin(phi)*(sx^2/10 + sy^2/10),    0, 0, 0;
                            0, 2*ur, 0, 0;
                            0,    0, 1, 0;
                            0,    0, 0, 1];
Qt = Lt*Qsys*Lt.';

% feed to ode45
% ode function handle for matrix input
function dPdt = VarianceDE(X, A, Q)
    X = reshape(X, size(A));  % Convert from n^2-by-1 to n-by-n
    dPdt = A*X + X*A.' + Q;   % Determine derivative
    dPdt = dPdt(:);           % Convert from n-by-n to n^2-by-1
end
[~, Plist] = ode45(@(t,Ps)VarianceDE(Ps, At, Qt), [tm-dt tm], estState.Pm);
Pp=reshape(Plist(end,:),[7,7]);

%% measurement update
% assign values to system states
[px, py, ~, ~, phi, ~, b] = deal(xp(1),xp(2),xp(3),xp(4),xp(5),xp(6),xp(7));

Hk =[(2*px + 2000)/(2*((px + 1000)^2 + (py - 1000)^2)^(1/2)), (2*py - 2000)/(2*((px + 1000)^2 + (py - 1000)^2)^(1/2)), 0, 0, 0, 0, 0;
         (2*px - 4000)/(2*((px - 2000)^2 + py^2)^(1/2)),                         py/((px - 2000)^2 + py^2)^(1/2), 0, 0, 0, 0, 0;
                        px/((py - 2000)^2 + px^2)^(1/2),          (2*py - 4000)/(2*((py - 2000)^2 + px^2)^(1/2)), 0, 0, 0, 0, 0;
                                                      0,                                                       0, 0, 0, 1, 0, 1;
                                                      0,                                                       0, 0, 0, 1, 0, 0];

Mk = eye(5);

% measurement equations
za = sqrt((px-estConst.pos_radioA(1))^2+(py-estConst.pos_radioA(2))^2)+wa;
zb = sqrt((px-estConst.pos_radioB(1))^2+(py-estConst.pos_radioB(2))^2)+wb;
zc = sqrt((px-estConst.pos_radioC(1))^2+(py-estConst.pos_radioC(2))^2)+wc;
zg = phi+b+wg;
zn = phi+wn;

% our nonlinear prediction:hk(xp,0)
hkxp = [za zb zc zg zn]';

% Kalman gain:
Kk = Pp*Hk'/(Hk*Pp*Hk' + Mk*R*Mk');

% update mean
% no measurement at initial step, use xp
if tm==0
    estState.xm = xp;
else
    % give more trust to zc if possible
    if isinf(sense(3))
        %sense(3)=zc; %test here
        Hk(3,:) = [];
        Mk(3,:) = [];
        Mk(:,3) = [];
        R(3,:) = [];
        R(:,3) = [];
        sense(3) = [];
        hkxp(3) = [];
        Kk = Pp*Hk'/(Hk*Pp*Hk' + Mk*R*Mk');
    end
    estState.xm = xp + Kk*(sense' - hkxp);
end

estState.Pm = (eye(7)-Kk*Hk)*Pp;
P_diag= diag(estState.Pm);

%% Get resulting estimates and variances

% Output quantities
posEst = estState.xm(1:2)';    % 1x2 matrix
linVelEst = estState.xm(3:4)';
oriEst = estState.xm(5);
windEst = estState.xm(6);
driftEst = estState.xm(7);
% 
posVar = P_diag(1:2)';    % 1x2 matrix
linVelVar = P_diag(3:4)'; % 1x2 matrix
oriVar = P_diag(5); 
windVar = P_diag(6);
driftVar = P_diag(7); 

end
