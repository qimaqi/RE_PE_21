function [postParticles] = Estimator(prevPostParticles, sens, act, estConst, km)
% The estimator function. The function will be called in two different
% modes: If km==0, the estimator is initialized. If km > 0, the
% estimator does an iteration for a single sample time interval using the 
% previous posterior particles passed to the estimator in
% prevPostParticles and the sensor measurement and control inputs.
%
% Inputs:
%   prevPostParticles   previous posterior particles at time step k-1
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%                           
%   sens                Sensor measurement z(k), scalar
%
%   act                 Control inputs u(k-1), [1x2]-vector
%                       act(1): u_f, forward control input
%                       act(2): u_phi, angular control input
%
%   estConst            estimator constants (as in EstimatorConst.m)
%
%   km                  time index k, scalar
%                       corresponds to continous time t = k*Ts
%                       If km==0 initialization, otherwise estimator
%                       iteration step.
%
% Outputs:
%   postParticles       Posterior particles at time step k
%                       The fields of the struct are [1xN_particles]-vector
%                       (N_particles is the number of particles) 
%                       corresponding to: 
%                       .x_r: x-locations of the robot [m]
%                       .y_r: y-locations of the robot [m]
%                       .phi: headings of the robot [rad]
%                       .kappa: wall offset [m]
%
%
% Class:
% Recursive Estimation
% Spring 2021
% Programming Exercise 2
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch

% Set number of particles:
N_particles = 50; % obviously, you will need more particles than 10.

%% Mode 1: Initialization
if (km == 0)
    % Do the initialization of your estimator here!
    pA = estConst.pA;
    pB = estConst.pB;
    d = estConst.d;
    phi_0 = estConst.phi_0;
    l = estConst.l;
    
    r = d * sqrt(rand(1,N_particles));
    theta = rand(1,N_particles) * 2 * pi;

    center_seed = rand(1,N_particles);
    centerX = zeros(1,10);
    centerY = zeros(1,10);
    
    for particle_index = 1:N_particles
        if center_seed(particle_index)<=0.5
            centerX = pA(1) ;
            centerY = pA(2) ;
        else
            centerX = pB(1);
            centerY = pB(2);
            %disp('using pB');
            %disp('using r');
            %disp(r);
        end
    end
    
    

    xr = centerX + r .* cos(theta);
    yr = centerY + r .* sin(theta);

    phi_r = -phi_0 + 2*phi_0*rand(1,N_particles);
    l_r = -l + 2*l*rand(1,N_particles);
    
    postParticles.x_r = xr;  % 1xN_particles matrix
    postParticles.y_r = yr; % 1xN_particles matrix
    postParticles.phi = phi_r ; % 1xN_particles matrix
    postParticles.kappa = l_r; % 1xN_particles matrix
    
    % and leave the function
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement your estimator here!
sigma_phi = estConst.sigma_phi;
sigma_f = estConst.sigma_f;

v_f = -sigma_f/2 + sigma_f*rand(1,N_particles);
v_phi = -sigma_phi/2 + sigma_phi*rand(1,N_particles);

x_r_prev = prevPostParticles.x_r;
y_r_prev = prevPostParticles.y_r;
phi_r_prev = prevPostParticles.phi;
kappa_r_prev = prevPostParticles.kappa;

contour = estConst.contour;
map_x =contour(:,1);
map_y =contour(:,2);
epsilon = estConst.epsilon;

% Prior Update:
x_r_p_hat = x_r_prev + (act(1)+v_f).*cos(phi_r_prev);
y_r_p_hat = y_r_prev + (act(1)+v_f).*sin(phi_r_prev);
phi_r_p_hat = phi_r_prev + act(2) + v_phi;
kappa_r_p_hat = kappa_r_prev;

% Posterior Update:

lineseg_x_start = x_r_p_hat;
lineseg_y_start = y_r_p_hat;
lineseg_x_end = lineseg_x_start + 10*cos(phi_r_p_hat);
lineseg_y_end = lineseg_y_start + 10*sin(phi_r_p_hat);
intersect_points_x = zeros(1,N_particles);
intersect_points_y = zeros(1,N_particles);
for particle_index = 1:N_particles
    map_x(8) = kappa_r_p_hat(particle_index);
    map_x(9) = kappa_r_p_hat(particle_index);
    map_poly = polyshape(map_x,map_y);
    lineseg = [lineseg_x_start(particle_index),lineseg_y_start(particle_index);
             lineseg_x_end(particle_index),lineseg_y_end(particle_index)];
    [in,~] = intersect(map_poly,lineseg);
    if size(in,1)>=2
        intersect_points_x(particle_index) = in(2,1);
        intersect_points_y(particle_index) = in(2,2);
    else
        warning('Particle outside the map')
    end

    
predict_distance = sqrt((x_r_p_hat-intersect_points_x).^2+(y_r_p_hat-intersect_points_y).^2);

sense_repeat = repmat(sens,1,N_particles);
error = sense_repeat - predict_distance;
p_post = zeros(1,N_particles);
% one possible bug
for particle_index = 1:N_particles
    if abs(error(particle_index))>=2*epsilon && abs(error(particle_index))<=3*epsilon
        p_post(particle_index) = 1/(5*epsilon) - 2/(5*epsilon^2) * abs(2.5*epsilon - abs(error(particle_index)));
    elseif abs(error(particle_index))<=2*epsilon
        p_post(particle_index) =  2/(5*epsilon) -  1/(5*epsilon^2) * abs(error(particle_index));
    end
end

if(sum(p_post) > 0)
    p_post = p_post/sum(p_post);
else
    p_post = ones(1,N_particles)*1/N_particles;
    % warning('Particle weights were all zero')
end

postCumSum = cumsum(p_post);
x_m_hat = zeros(1,N_particles);
y_m_hat = zeros(1,N_particles);
phi_m_hat = zeros(1,N_particles);
kappa_m_hat = zeros(1,N_particles);

for i = 1:N_particles
    randNumber = rand;
    ind = find(postCumSum >= randNumber,1,'first');
    x_m_hat(i) = x_r_p_hat(ind);
    y_m_hat(i) = y_r_p_hat(ind);
    phi_m_hat(i) = phi_r_p_hat(ind);
    kappa_m_hat(i) = kappa_r_p_hat(ind);
end

postParticles.x_r = x_m_hat;
postParticles.y_r = y_m_hat;
postParticles.phi = phi_m_hat;
postParticles.kappa = kappa_m_hat;

end % end estimator