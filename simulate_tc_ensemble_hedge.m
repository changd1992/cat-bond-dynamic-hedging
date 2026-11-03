function [results, path_results] = ...
    simulate_tc_ensemble_hedge(lambda, params)

    nsample = params.nsample;

    for i = 1:nsample

        % Simulate TC path
        path = simulate_tc_path(lambda, params);

        % Calculate CAT bond values
        path.bond = calculate_catbond_path(path, params);

        % Calculate delta hedge
        path.delta_hedge = calculate_delta_hedge(path, params);

        % Calculate static minimum variance hedge
        path.smv_hedge = calculate_mvhr(path, params);

        % Store entire path including histories
        if i == 1
            path_results = repmat(path, nsample, 1);
        end

        path_results(i) = path;

    end

    % Aggregate results
    results = aggregate_hedging_results(path_results);

end

function hedge = calculate_delta_hedge(path, params)

    R = path.R_history;
    V = path.V_history;
    t = path.time_history;
    n = numel(t);

    % Baseline bond values
    bond = path.bond;

    % Reconstruct common GBM innovations from baseline path
    dt = diff(t);
    z = zeros(n-1, 1);

    for k = 1:n-1
        z(k) = (log(R(k+1)/R(k)) ...
            - (params.mu_R - 0.5*params.sigma_R^2)*dt(k)) ...
            / (params.sigma_R*sqrt(dt(k)));
    end

    % Perturb initial R and V by +/- 1%
    eps = 0.01;

    R_plus = zeros(n,1); R_minus = zeros(n,1);
    V_plus = zeros(n,1); V_minus = zeros(n,1);

    R_plus(1) = R(1)*(1+eps); R_minus(1) = R(1)*(1-eps);
    V_plus(1) = V(1)*(1+eps); V_minus(1) = V(1)*(1-eps);

    for k = 2:n
        growth_R = exp( ...
            (params.mu_R - 0.5*params.sigma_R^2)*dt(k-1) ...
            + params.sigma_R*sqrt(dt(k-1))*z(k-1));

        growth_V = exp( ...
            (params.mu_V - 0.5*params.sigma_V^2)*dt(k-1) ...
            + params.sigma_V*sqrt(dt(k-1))*z(k-1));

        R_plus(k) = R_plus(k-1)*growth_R; R_minus(k) = R_minus(k-1)*growth_R;
        V_plus(k) = V_plus(k-1)*growth_V; V_minus(k) = V_minus(k-1)*growth_V;
    end

    % CHI trajectories
    CHI_plus = calculate_chi(R_plus, V_plus, params);
    CHI_minus = calculate_chi(R_minus, V_minus, params);

    % Perturbed terminal payoffs
    if path.hit
        payoff_plus = params.Face ...
            - max(0, CHI_plus(end)-params.short_X) ...
            + max(0, CHI_plus(end)-params.long_X);

        payoff_minus = params.Face ...
            - max(0, CHI_minus(end)-params.short_X) ...
            + max(0, CHI_minus(end)-params.long_X);
    else
        payoff_plus = params.Face;
        payoff_minus = params.Face;
    end

    % Discount perturbed payoffs
    discount = exp(-params.rate*(t(end)-t));

    P_plus = payoff_plus .* discount; P_minus = payoff_minus .* discount;

    % Finite-difference delta
    delta = (P_plus - P_minus) ./ (CHI_plus - CHI_minus);

    % Replicate original portfolio definition
    portfolio = path.chi_history - delta .* bond.P_bond;

    % Return results
    hedge.CHI = path.chi_history;
    hedge.CHI_plus = CHI_plus;
    hedge.CHI_minus = CHI_minus;

    hedge.P_bond = bond.P_bond;
    hedge.P_plus = P_plus;
    hedge.P_minus = P_minus;

    hedge.delta = delta;
    hedge.portfolio = portfolio;
    hedge.time = t;

end

function hedge = calculate_mvhr(path_results, params)

    nsample = numel(path_results);
    
    % Terminal CHI and bond values
    CHI_T = zeros(nsample,1);
    P_bond = zeros(nsample,1);
    
    for i = 1:nsample
        path = path_results(i);
    
        CHI_T(i) = path.chi_history(end);
        P_bond(i) = path.bond.P_bond(end);
    end
    
    % Static minimum-variance hedge ratio
    a = cov(CHI_T, P_bond);
    minhedge2 = a(1,2) / a(2,2);
    
    % Terminal hedged portfolio
    Portfolio2 = CHI_T - minhedge2 .* P_bond;
    
    % Variance reduction
    varianceCHI1 = var(CHI_T);
    varianceP2 = var(Portfolio2);
    
    Ratio2 = (varianceCHI1 - varianceP2) / varianceCHI1;
    
    % Return results
    hedge.CHI_T = CHI_T;
    hedge.P_bond = P_bond;
    hedge.delta = minhedge2;
    hedge.portfolio = Portfolio2;
    hedge.varianceCHI = varianceCHI1;
    hedge.varianceP = varianceP2;
    hedge.ratio = Ratio2;

end

function bond = calculate_catbond_path(path, params)

    CHI = path.chi_history;
    t = path.time_history;

    % Terminal CAT bond payoff
    if path.hit
        short_payoff = max(0, CHI(end) - params.short_X);
        long_payoff = max(0, CHI(end) - params.long_X);
    else
        short_payoff = 0;
        long_payoff = 0;
    end

    terminal_payoff = params.Face - short_payoff + long_payoff;

    % Discount from termination time to each observation
    discount = exp(-params.rate * (t(end) - t));

    bond.P_short = short_payoff .* discount;
    bond.P_long = long_payoff .* discount;
    bond.P_bond = terminal_payoff .* discount;

    bond.terminal_payoff = terminal_payoff;
    bond.time = t;

end

function results = aggregate_hedging_results(path_results)

    nsample = numel(path_results);
    
    variance_CHI = NaN(nsample,1);
    variance_portfolio = NaN(nsample,1);
    performance = NaN(nsample,1);
    
    for i = 1:nsample
    
        path = path_results(i);
    
        if ~isfield(path, 'hedge') || isempty(path.hedge)
            continue;
        end
    
        CHI = path.hedge.CHI;
        portfolio = path.hedge.portfolio;
    
        if numel(CHI) < 2
            continue;
        end
    
        variance_CHI(i) = var(CHI);
        variance_portfolio(i) = var(portfolio);
    
        if variance_CHI(i) > 0
            performance(i) = 1 ...
                - variance_portfolio(i)/variance_CHI(i);
        end
    
    end
    
    results.variance_CHI = variance_CHI;
    results.variance_portfolio = variance_portfolio;
    results.performance = performance;
    
    results.mean_effectiveness = mean(performance,'omitnan');
    results.median_effectiveness = median(performance,'omitnan');

end

function CHI = calculate_chi(R, V, params)

    % Calculate CME Hurricane Index
    CHI = (V/params.V0).^3 + 1.5*(R/params.R0).*(V/params.V0).^2;

end