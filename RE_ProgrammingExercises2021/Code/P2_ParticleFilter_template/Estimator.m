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
N_particles = 2000;

%% Mode 1: Initialization
% robustness evaluation will vary d, l and epsilon
if (km == 0)
    % initialization of estimator
    pA = estConst.pA;
    pB = estConst.pB;
    d = estConst.d;
    phi_0 = estConst.phi_0;
    l = estConst.l;
    
    % uniform distri. inside a circle
    r = d * sqrt(rand(1,N_particles)); 
    theta = rand(1,N_particles) * 2 * pi;
    
    % choose a initial circle
    center_seed = rand(1,N_particles);
    centerX = zeros(1,N_particles);
    centerY = zeros(1,N_particles);

    for particle_index = 1:N_particles
        if center_seed(particle_index)<=0.5
            centerX(particle_index) = pA(1);
            centerY(particle_index) = pA(2);
        else
            centerX(particle_index) = pB(1);
            centerY(particle_index) = pB(2);
        end
    end

    xr = centerX + r .* cos(theta);
    yr = centerY + r .* sin(theta);

    phi_r = -phi_0 + 2*phi_0*rand(1,N_particles);
    l_r = -l + 2*l*rand(1,N_particles);
    
    postParticles.x_r = xr;     % 1xN_particles matrix
    postParticles.y_r = yr;     % 1xN_particles matrix
    postParticles.phi = phi_r ; % 1xN_particles matrix
    postParticles.kappa = l_r;  % 1xN_particles matrix
    
    return;
end % end init

%% Mode 2: Estimator iteration.
% If km > 0, we perform a regular update of the estimator.

% Implement estimator
sigma_phi = estConst.sigma_phi;
sigma_f = estConst.sigma_f;

rng shuffle

v_f = -sigma_f/2 + sigma_f*rand(1,N_particles);
v_phi = -sigma_phi/2 + sigma_phi*rand(1,N_particles);

x_r_prev = prevPostParticles.x_r;
y_r_prev = prevPostParticles.y_r;
phi_r_prev = prevPostParticles.phi;
kappa_r_prev = prevPostParticles.kappa;

contour = estConst.contour;
epsilon = estConst.epsilon;

% Prior Update:
% according to process model
x_r_p_hat = x_r_prev + (act(1)+v_f).*cos(phi_r_prev);
y_r_p_hat = y_r_prev + (act(1)+v_f).*sin(phi_r_prev);
phi_r_p_hat = phi_r_prev + act(2) + v_phi;
kappa_r_p_hat = kappa_r_prev;

% Posterior Update:
% compute the distance to the intersection
distance_from_wall = zeros(1,N_particles);

for i = 1:N_particles
    contour(8,1) = kappa_r_p_hat(i);
    contour(9,1) = kappa_r_p_hat(i);
    distance_from_wall(i) = compute_distance_from_wall(x_r_p_hat(i),y_r_p_hat(i),phi_r_p_hat(i),contour);
end

sense_repeat = repmat(sens,1,N_particles);
error = sense_repeat - distance_from_wall;
p_post = zeros(1,N_particles);
p_post_gaussian = zeros(1,N_particles);

for particle_index = 1:N_particles
    if abs(error(particle_index))>=2*epsilon && abs(error(particle_index))<=3*epsilon
        p_post(particle_index) = 1/(5*epsilon) - 2/(5*epsilon^2) * abs(2.5*epsilon - abs(error(particle_index)));
    elseif abs(error(particle_index))<=2*epsilon
        p_post(particle_index) =  2/(5*epsilon) -  1/(5*epsilon^2) * abs(error(particle_index));
    %elseif abs(error(particle_index))>3*epsilon
    %    
    %    p_post(particle_index) = normpdf(error(particle_index),0,normal_var);  %1/(2*pi*normal_var)*exp(-(abs(error(particle_index))-0).^2/(2*normal_var^2));
        %p_post(particle_index)=0;
    end
    normal_var=0.10;
    p_post_gaussian(particle_index) = 1/(2*pi*normal_var)*exp(-(abs(error(particle_index))-0).^2/(2*normal_var^2));
end
% if>3*epsilon, strategy 0: 0.2752, 7/50 outliers>0.5
%                           0.2052, 5/50 outliers>0.5 (without update phi)
%        strategy gaussian: 0.1802  4/50 outliers>0.5 (without update phi)  

% allocate weights
if(sum(p_post) > 0)
    p_post = p_post./sum(p_post);
else
    %normal_var=0.10;
    %p_post_gaussian = normpdf(error,0,normal_var);
    p_post = p_post_gaussian./sum(p_post_gaussian);
    %ones(1,N_particles)*1/N_particles;
%   warning('Particle weights were all zero')
end

% resampling
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

K_tune=0.0005;

% add roughening
% not much info from estimate to update kappa, causes outliers 
[postParticles.x_r,postParticles.y_r,postParticles.phi,postParticles.kappa] = roughening_samples(x_m_hat,y_m_hat,phi_m_hat,kappa_m_hat,K_tune,N_particles);

end % end estimator

function distance = compute_distance_from_wall(x0,y0,phi,contour)
    % parameters
    length = size(contour,1);
    alpha = zeros(1,length);
    beta  = zeros(1,length);
    % loop to the start
    contour_aug = [contour; contour(1,:)];
    % calculate by parameter equations
    for i = 1:length
        x1 = contour_aug(i,1);
        y1 = contour_aug(i,2);
        x2 = contour_aug(i+1,1);
        y2 = contour_aug(i+1,2);
        den = sin(phi)*(x2-x1)+cos(phi)*(y1-y2);
        alpha(i) = (sin(phi)*(x2-x0)+cos(phi)*(y0-y2))/den;
        beta(i)  = ((y2-y1)*(x0-x2)+(x1-x2)*(y0-y2))/den;
    end
    % determine the wall intersected
    intersect_pos = alpha >= 0 & alpha <= 1 & beta >= 0;
    distance_list = beta.*intersect_pos;
    distance = min(distance_list(distance_list~=0));
    if isempty(distance)
        distance = inf;
    end
end

function [x_r,y_r,phi,kappa]=roughening_samples(x_m_hat,y_m_hat,phi_m_hat,kappa_m_hat,K_tune,N_particles)
    x_r = x_m_hat + sqrt(K_tune*(max(x_m_hat)-min(x_m_hat))*N_particles^(-1/4))*randn(1,N_particles);
    y_r = y_m_hat + sqrt(K_tune*(max(y_m_hat)-min(y_m_hat))*N_particles^(-1/4))*randn(1,N_particles);
    phi = phi_m_hat ;%+ sqrt(K_tune*(max(phi_m_hat)-min(phi_m_hat))*N_particles^(-1/4))*randn(1,N_particles);
    kappa = kappa_m_hat + sqrt(K_tune*(max(kappa_m_hat)-min(kappa_m_hat))*N_particles^(-1/4))*randn(1,N_particles);
end