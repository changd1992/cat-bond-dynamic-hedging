function params = initialize_parameters(nsample, c, heading, R, V)

    % Monte Carlo simulation
    params.nsample = nsample;
    params.max_steps = 3000;
    params.timelength = 2;
    
    % Initial TC position
    params.lat_min = 10;
    params.lat_max = 20;
    params.lon_min = -60;
    params.lon_max = -40;
    params.earthr = 6356;
    
    % Initial TC characteristics
    params.R0 = 60;
    params.V0 = 74;
    params.R_init = R;
    params.V_init = V;
    params.c_begin = c;
    params.h_begin = heading;
    
    % TC intensity dynamics (GBM)
    params.mu_R = 0.01;
    params.sigma_R = 0.02;
    params.mu_V = 0.01;
    params.sigma_V = 0.02;
    
    % TC translation speed coefficients
    params.a1 = 2.5e-4;
    params.a2 = 5e-5;
    params.a3 = 5e-5;
    params.a4 = 5e-5;
    params.a5 = 5e-5;
    params.a6 = 5e-5; % DOUBLE CHECK
    
    % TC heading coefficients
    params.b1 = 2.5e-3;
    params.b2 = 5e-4;
    params.b3 = 5e-4;
    params.b4 = 5e-4;
    params.b5 = 5e-4;
    params.b6 = 5e-4;
    params.b7 = 5; % DOUBLE CHECK
    
    % CAT bond contract
    params.K = 10;
    params.Face = 8;
    params.short_X = params.K;
    params.long_X = params.K + params.Face;
    
    % Discounting
    params.spread = 0.08925;
    params.rate = 0.0101;
    params.rate1 = 0.0101;
    params.rate2 = params.rate1 + params.spread;
    
    % Geographic termination boundaries
    params.lon_exit_min = -100;
    params.lon_exit_max = -10;
    params.lat_exit_min = 0;
    params.lat_exit_max = 42;
    
    % Landfall boundaries
    params.atlantic_lon = -75;
    params.atlantic_lat_min = 30;
    params.atlantic_lat_max = 42;
    
    params.gulf_lat = 30;
    params.gulf_lon_min = -100;
    params.gulf_lon_max = -75;

end