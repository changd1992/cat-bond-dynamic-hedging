function [results, path_results] = simulate_tc_ensemble(lambda, params)

    % Initialize results structure
    nsample = params.nsample;

    % Monte Carlo simulation
    for i = 1:nsample

        % Simulate one TC trajectory
        path = simulate_tc_path(lambda, params);

        % Calculate CAT bond payoff
        path.payoff = calculate_catbond_payoff( ...
            path.CHI, path.hit, params);

        % Discount payoff
        path.discounted1 = discount_payoff( ...
            path.payoff, path.time, params.rate1);

        path.discounted2 = discount_payoff( ...
            path.payoff, path.time, params.rate2);

        % Store entire path including histories
        if i == 1
            path_results = repmat(path, nsample, 1);
        end

        path_results(i) = path;

    end

    % Aggregate Monte Carlo results
    results = aggregate_mc_results(path_results);

end

function results = aggregate_mc_results(path_results)

    % Extract path-level results
    discounted1 = [path_results.discounted1];
    discounted2 = [path_results.discounted2];
    hits = [path_results.hit];
    steps = [path_results.steps];
    times = [path_results.time];
    
    % Monte Carlo sample size
    results.nsample = numel(path_results);
    
    % CAT bond prices
    results.bondPrice1 = mean(discounted1);
    results.bondPrice2 = mean(discounted2);
    
    % Monte Carlo standard errors
    results.SE1 = std(discounted1) / sqrt(results.nsample);
    results.SE2 = std(discounted2) / sqrt(results.nsample);
    
    % Landfall statistics
    results.hit_count = sum(hits);
    results.hit_probability = mean(hits);
    
    % Conditional landfall time statistics
    if any(hits)
        results.mean_landfall_time = mean(times(hits));
    else
        results.mean_landfall_time = NaN;
    end
    
    % Simulation diagnostics
    results.mean_steps = mean(steps);
    results.mean_stopping_time = mean(times);

end

function payoff = calculate_catbond_payoff(CHI, hit, params)

    % Full principal repayment if no landfall
    if ~hit
        payoff = params.Face;
        return;
    end
    
    % Calculate option payoffs
    P_short = max(0, CHI - params.short_X);
    P_long = max(0, CHI - params.long_X);
    
    % Calculate CAT bond payoff
    payoff = params.Face - P_short + P_long;

end

function discounted_payoff = discount_payoff(payoff, time, rate)

    % Discount payoff to present value
    discounted_payoff = payoff * exp(-rate * time);

end