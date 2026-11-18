function path = simulate_tc_path(lambda, params)

    % Initialize TC state
    R = params.R_init; V = params.V_init;
    c = params.c_begin; heading = params.h_begin;
    
    lat = params.lat_min + ...
        (params.lat_max - params.lat_min) * rand;
    lon = params.lon_min + ...
        (params.lon_max - params.lon_min) * rand;
    
    % Initialize simulation
    t = 0; n = 1;
    hit = false; terminated = false;
    
    % Store initial state
    R_history = zeros(params.max_steps, 1); V_history = zeros(params.max_steps, 1);
    lat_history = zeros(params.max_steps, 1); lon_history = zeros(params.max_steps, 1);
    c_history = zeros(params.max_steps, 1); h_history = zeros(params.max_steps, 1);
    chi_history = zeros(params.max_steps, 1);
    time_history = zeros(params.max_steps, 1);
    
    R_history(1) = R; V_history(1) = V;
    lat_history(1) = lat; lon_history(1) = lon;
    c_history(1) = c; h_history(1) = heading;

    CHI = calculate_chi(R, V, params);
    chi_history(1) = CHI;

    % Initialize previous heading
    heading_prev = heading;
    
    % Simulate TC evolution
    while t < params.timelength && n < params.max_steps
    
        % Sample next information arrival
        dt = sample_information_arrival(lambda);
    
        % Stop at simulation horizon
        if t + dt > params.timelength
            break;
        end
    
        % Update TC intensity and CHI
        [R, V] = update_tc_intensity(R, V, dt, params);
        CHI = calculate_chi(R, V, params);

        heading_old = heading;
    
        % Update TC track
        [lat, lon, c, heading] = update_tc_track( ...
            lat, lon, c, heading, heading_prev, dt, params);
    
        heading_prev = heading_old;
    
        % Advance simulation
        t = t + dt; n = n + 1;
    
        % Store state
        R_history(n) = R; V_history(n) = V;
        lat_history(n) = lat; lon_history(n) = lon;
        c_history(n) = c; h_history(n) = heading;
        chi_history(n) = CHI;
        time_history(n) = t;
    
        % Check landfall or geographic exit
        [hit, terminated] = check_tc_termination(lat, lon, params);
    
        if terminated
            break;
        end
    
    end
    
    % Return simulation results
    path.R_end = R; path.V_end = V;
    path.lat_end = lat; path.lon_end = lon;
    path.time = t; path.steps = n;
    path.CHI = CHI;
    path.hit = hit; path.terminated = terminated;
    
    path.R_history = R_history(1:n); path.V_history = V_history(1:n);
    path.lat_history = lat_history(1:n); path.lon_history = lon_history(1:n);
    path.c_history = c_history(1:n); path.h_history = h_history(1:n);
    path.chi_history = chi_history(1:n);
    path.time_history = time_history(1:n);

end

function [lat, lon, c, heading] = update_tc_track( ...
    lat, lon, c, heading, heading_prev, dt, params)

    % Update TC translation speed
    d_logc = params.a1 + params.a2*lat + params.a3*lon ...
        + params.a4*log(c) + params.a5*heading ...
        + params.a6*randn;
    
    c_new = exp(log(c) + d_logc);
    
    % Update TC heading
    d_heading = params.b1 + params.b2*lat + params.b3*lon ...
        + params.b4*c + params.b5*heading ...
        + params.b6*heading_prev + params.b7*randn;
    
    heading_new = heading + d_heading;
    
    % Update geographic position using previous heading
    distance = c_new * 3600 * dt * 720 / 1000;
    
    lon = lon - rad2deg(distance / params.earthr) ...
        * sin(deg2rad(heading));
    
    lat = lat + rad2deg(distance / params.earthr) ...
        * cos(deg2rad(heading));
    
    % Return updated speed and heading
    c = c_new; heading = heading_new;

end

function [R, V] = update_tc_intensity(R, V, dt, params)

    % Generate random innovation (shared by R and V)
    z = randn;
    
    % Update R using geometric Brownian motion
    R = R * exp((params.mu_R - 0.5*params.sigma_R^2)*dt ...
        + params.sigma_R*sqrt(dt)*z);
    
    % Update V using geometric Brownian motion
    V = V * exp((params.mu_V - 0.5*params.sigma_V^2)*dt ...
        + params.sigma_V*sqrt(dt)*z);

end

function CHI = calculate_chi(R, V, params)

    % Calculate CME Hurricane Index
    CHI = (V/params.V0)^3 + 1.5*(R/params.R0)*(V/params.V0)^2;

end

function [hit, terminated] = check_tc_termination(lat, lon, params)

    % Initialize outputs
    hit = false;
    terminated = false;
    
    % Check geographic domain exit
    if lon > params.lon_exit_max || lon < params.lon_exit_min || ...
            lat > params.lat_exit_max || lat < params.lat_exit_min
    
        terminated = true;
        return;
    end

    % Check landfall
    atlantic_hit = lon <= params.atlantic_lon && ...
    lat >= params.atlantic_lat_min && lat <= params.atlantic_lat_max;

    gulf_hit = lat >= params.gulf_lat && ...
        lon >= params.gulf_lon_min && lon <= params.gulf_lon_max;
    
    if atlantic_hit || gulf_hit
        hit = true; terminated = true;
    end

end

function dt = sample_information_arrival(lambda)

    % Sample exponential waiting time between information arrivals
    dt = -log(rand) / lambda;

end
