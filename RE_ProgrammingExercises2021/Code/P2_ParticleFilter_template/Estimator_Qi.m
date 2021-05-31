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
    p_post(particle_index) = calculate_post(error(particle_index),epsilon);
    normal_var=0.3 ;
    p_post_gaussian(particle_index) = 1/(2*pi*normal_var)*exp(-(abs(error(particle_index))-0).^2/(2*normal_var^2));
end

if(sum(p_post) > 0)
    p_post = p_post./sum(p_post);
else  
    string = "Zero Probability after Measurement Update at " + int2str(km) + ". Resampling ...";
    disp(string);
    % one time random resample kappa
    resample_kappa = -estConst.l + 2*estConst.l*rand(1,N_particles); 
    kappa_r_p_hat = resample_kappa;
    distance_from_wall = zeros(1,N_particles);
    for i = 1:N_particles
        contour(8,1) = kappa_r_p_hat(i);
        contour(9,1) = kappa_r_p_hat(i);
        distance_from_wall(i) = compute_distance_from_wall(x_r_p_hat(i),y_r_p_hat(i),phi_r_p_hat(i),contour);
    end
    error = sense_repeat - distance_from_wall;
    for particle_index = 1:N_particles
        p_post(particle_index) = calculate_post(error(particle_index),epsilon);
    end
    if(sum(p_post) > 0)
        p_post = p_post./sum(p_post);
    else
         p_post = p_post_gaussian./sum(p_post_gaussian);
         disp('kappa resample do not work');
    end 
 
    
    %%% 10 times resample M %%%
%     M = 10;
%     for i = 1:N_particles
%         resample_kappa = -estConst.l + 2*estConst.l*rand(1,M);
%         temp_p_post = zeros(1,M);
%         for j = 1:M
%            contour(8,1) = resample_kappa(j);
%            contour(9,1) = resample_kappa(j);
%            distance_from_wall = compute_distance_from_wall(x_r_p_hat(i),y_r_p_hat(i),phi_r_p_hat(i),contour);
%            error = sens - distance_from_wall;
%            temp_p_post(j) = calculate_post(error,epsilon);    
%         end
%         [p_post(i),index] = max(temp_p_post);
%         kappa_r_p_hat(i) = kappa_r_p_hat(index);
%     end
%     if(sum(p_post) > 0)
%         p_post = p_post./sum(p_post);
%     else
%          p_post = p_post_gaussian./sum(p_post_gaussian);
%          %disp('kappa resample do not work');
%     end
    %%% 10 times resample M %%%
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
    phi = phi_m_hat;% + sqrt(K_tune*(max(phi_m_hat)-min(phi_m_hat))*N_particles^(-1/4))*randn(1,N_particles);
    kappa = kappa_m_hat + 0.8*sqrt(K_tune*(max(kappa_m_hat)-min(kappa_m_hat))*N_particles^(-1/4))*randn(1,N_particles);
end

function p_post = calculate_post(error,epsilon)
    error = abs(error);
    if error>=2*epsilon && error<=3*epsilon
        p_post = 1/(5*epsilon) - 2/(5*epsilon^2) * abs(2.5*epsilon - error);
    elseif error <= 2*epsilon
        p_post = 2/(5*epsilon) -  1/(5*epsilon^2) * error;
    else
        p_post = 0;
    end
end
