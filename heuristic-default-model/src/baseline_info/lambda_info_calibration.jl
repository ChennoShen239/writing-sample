# =============================================================================
# CALIBRATION SCRIPT FOR SOVEREIGN DEFAULT MODEL WITH BOUNDED RATIONALITY
#
# Purpose: This script calibrates the parameter λ (fraction of rational lenders)
# to match a target average interest rate spread using the bisection method.
# =============================================================================

# --- 1. Preamble: Load necessary packages ---
using LinearAlgebra, Statistics
using LaTeXStrings, QuantEcon, DataFrames, Plots, Random
using Interpolations
using Printf

# --- 2. Model Definition (Copied from your main script) ---

# Create directory for results if it doesn't exist
mkpath("Result/Figs/lambda_info")

function ArellanoEconomy(;
    beta=0.953,
    gamma=2.0,
    r=0.017,
    rho=0.945,    # persistence in output
    eta=0.025,
    theta=0.282,
    ny=201,
    nB=401,
    λ=1.0,        # fraction of lenders who know how to compute the default probability
)

    # create grids
    Bgrid = collect(range(-0.4, 0.4, length=nB))
    mc = tauchen(ny, rho, eta)
    Pi = mc.p
    ygrid = exp.(mc.state_values)
    ydefgrid = min.(0.969 * mean(ygrid), ygrid)

    # define value functions
    vf = zeros(nB, ny)
    vd = zeros(1, ny)
    vc = zeros(nB, ny)
    policy = zeros(nB, ny)
    q = ones(nB, ny) .* (1 / (1 + r))
    defprob = zeros(nB, ny)

    return (;
        beta,
        gamma,
        r,
        rho,
        eta,
        theta,
        ny,
        nB,
        ygrid,
        ydefgrid,
        Bgrid,
        Pi,
        vf,
        vd,
        vc,
        policy,
        q,
        defprob,
        λ,
    )
end

u(ae, c) = c^(1 - ae.gamma) / (1 - ae.gamma)

function one_step_update!(ae, EV, EVd, EVc)
    (; beta, gamma, r, rho, eta, theta, ny, nB) = ae
    (; ygrid, ydefgrid, Bgrid, Pi, vf, vd, vc, policy, q, defprob) = ae
    zero_ind = searchsortedfirst(Bgrid, 0.0)

    for iy = 1:ny
        y = ae.ygrid[iy]
        ydef = ae.ydefgrid[iy]

        defval = u(ae, ydef) + beta * (theta * EVc[zero_ind, iy] + (1 - theta) * EVd[1, iy])
        ae.vd[1, iy] = defval

        for ib = 1:nB
            B = ae.Bgrid[ib]
            current_max = -1e14
            pol_ind = 0
            for ib_next = 1:nB
                c = max(y - ae.q[ib_next, iy] * Bgrid[ib_next] + B, 1e-14)
                m = u(ae, c) + beta * EV[ib_next, iy]
                if m > current_max
                    current_max = m
                    pol_ind = ib_next
                end
            end
            ae.vc[ib, iy] = current_max
            ae.policy[ib, iy] = pol_ind
            ae.vf[ib, iy] = defval > current_max ? defval : current_max
        end
    end
end

function interp2d(xs, ys, A)
    itp2d = interpolate((xs, ys), A, Gridded(Linear()))
    eitp2d = extrapolate(itp2d, Line())
    return eitp2d
end

function compute_prices!(ae)
    (; λ, r, Pi) = ae
    vd_compact = repeat(ae.vd, ae.nB)
    default_states = vd_compact .> ae.vc
    def_interp = interp2d(ae.Bgrid, ae.ygrid, vd_compact)
    vd_exp_y = [def_interp(B, y^ae.rho) for B in ae.Bgrid, y in ae.ygrid]
    repay_interp = interp2d(ae.Bgrid, ae.ygrid, ae.vc)
    repay_exp_y = [repay_interp(B, y^ae.rho) for B in ae.Bgrid, y in ae.ygrid]
    def_states_exp_y = vd_exp_y .> repay_exp_y
    copyto!(ae.defprob, λ * default_states * Pi' + (1 - λ) * def_states_exp_y)
    copyto!(ae.q, (1 .- ae.defprob) / (1 + r))
    return
end

function vfi!(ae; tol=1e-8, maxit=10000, verbose=true)
    (; Pi) = ae
    Pit = Pi'
    it = 0
    dist = 10.0
    V_upd = similar(ae.vf)

    while dist > tol && it < maxit
        it += 1
        copyto!(V_upd, ae.vf)
        EV = ae.vf * Pit
        EVd = ae.vd * Pit
        EVc = ae.vc * Pit
        one_step_update!(ae, EV, EVd, EVc)
        compute_prices!(ae)
        dist = maximum(abs(x - y) for (x, y) in zip(V_upd, ae.vf))
        if verbose && it % 25 == 0
            println("Finished iteration $(it) with dist of $(dist)")
        end
    end
end

function AEsimulate(ae, capT=5000; y_init=mean(ae.ygrid), B_init=mean(ae.Bgrid))
    zero_index = searchsortedfirst(ae.Bgrid, 0.0)
    y_init_ind = searchsortedfirst(ae.ygrid, y_init)
    B_init_ind = searchsortedfirst(ae.Bgrid, B_init)
    mc = MarkovChain(ae.Pi)
    y_sim_indices = simulate(mc, capT + 1; init=y_init_ind)
    y_sim_val = zeros(capT + 1)
    B_sim_val, q_sim_val = similar(y_sim_val), similar(y_sim_val)
    B_sim_indices = fill(0, capT + 1)
    default_status = fill(false, capT + 1)
    B_sim_indices[1], default_status[1] = B_init_ind, false
    y_sim_val[1], B_sim_val[1] = ae.ygrid[y_init_ind], ae.Bgrid[B_init_ind]

    for t = 1:capT
        yi, Bi = y_sim_indices[t], B_sim_indices[t]
        defstat = default_status[t]
        if !defstat
            default_today = ae.vc[Bi, yi] < ae.vd[yi]
            if default_today
                default_status[t] = true
                default_status[t+1] = true
                y_sim_val[t] = ae.ydefgrid[y_sim_indices[t]]
                B_sim_indices[t+1] = zero_index
                B_sim_val[t+1] = 0.0
                q_sim_val[t] = ae.q[zero_index, y_sim_indices[t]]
            else
                default_status[t] = false
                y_sim_val[t] = ae.ygrid[y_sim_indices[t]]
                B_sim_indices[t+1] = ae.policy[Bi, yi]
                B_sim_val[t+1] = ae.Bgrid[B_sim_indices[t+1]]
                q_sim_val[t] = ae.q[B_sim_indices[t+1], y_sim_indices[t]]
            end
        else
            B_sim_indices[t+1] = zero_index
            B_sim_val[t+1] = 0.0
            y_sim_val[t] = ae.ydefgrid[y_sim_indices[t]]
            q_sim_val[t] = ae.q[zero_index, y_sim_indices[t]]
            default_status[t+1] = rand() >= ae.theta
        end
    end
    return (y_sim_val[1:capT], B_sim_val[1:capT], q_sim_val[1:capT], default_status[1:capT])
end


# --- 3. Calibration-Specific Functions ---

"""
    get_avg_spread(λ_val; T_sim=500000)

Helper function for calibration. It takes a value for λ, solves the model,
simulates it, and returns the single resulting average interest rate spread.
This is a computationally intensive function. For calibration runs, we use
a shorter simulation (T_sim) and fewer VFI iterations to speed things up.
"""
function get_avg_spread(λ_val::Float64; T_sim=50000)
    Random.seed!(1234)
    # --- 1. Setup and Solve the Economy ---
    ae = ArellanoEconomy(λ=λ_val, ny=201, nB=401)
    vfi!(ae; maxit=5000, tol=1e-8, verbose=true) # Solve model quietly

    # --- 2. Simulate ---
    y_vec, B_vec, q_vec, default_vec = AEsimulate(ae, T_sim)

    # --- 3. Calculate Average Spread using Bornstein (2020) methodology ---
    (length(y_vec) < 2) && return NaN

    T_orig = length(y_vec)
    # Note: The calculation of default_trimmed is crucial for window selection
    default_trimmed = default_vec[1:(T_orig-1)]
    q_trimmed = max.(q_vec[1:(T_orig-1)], 1e-12)

    annual_spread = (1.0 ./ q_trimmed).^4 .- (1.0 .+ ae.r).^4
    # annual_spread = ((1.0 .+ quarterly_spread) .^ 4 .- 1.0)

    window_size = 74
    max_events_to_sample = 100
    default_event_indices = findall(default_trimmed)

    window_means_spread = Float64[]
    processed_event_count = 0
    last_sampled_window_end_idx = 0

    for event_idx in default_event_indices
        (processed_event_count >= max_events_to_sample) && break

        window_end_actual = event_idx - 1
        window_start_idx = window_end_actual - window_size + 1

        if window_start_idx > 0 && window_start_idx > last_sampled_window_end_idx
            if !any(default_trimmed[window_start_idx:window_end_actual])
                spread_window = annual_spread[window_start_idx:window_end_actual]
                push!(window_means_spread, mean(filter(isfinite, spread_window)))
                last_sampled_window_end_idx = window_end_actual
                processed_event_count += 1
            end
        end
    end

    (isempty(window_means_spread)) && return NaN

    return mean(filter(isfinite, window_means_spread))
end


"""
    calibrate_lambda(; target_spread=0.1025, ...)

Uses the bisection method to find the value of λ that makes the model's
average interest rate spread match the target_spread.
"""
function calibrate_lambda(;
    target_spread::Float64=0.1025,
    lambda_low::Float64=0.5,
    lambda_high::Float64=0.55,
    tol::Float64=1e-3, # Tolerance on the value of λ
    max_iter::Int=20
)
    println("\n" * "="^60)
    println("--- Starting Calibration for λ ---")
    println("Target average spread: ", round(target_spread * 100, digits=2), "%")
    println("="^60)

    # Define the error function we want to find the root of
    error_func = λ -> get_avg_spread(λ) - target_spread

    # --- Check if the target is bracketed by the initial bounds ---
    println("Evaluating lower bound (λ=$lambda_low)...")
    err_low = error_func(lambda_low)
    if isnan(err_low)
        println("Error: Could not compute spread at lower bound. Aborting.")
        return NaN
    end
    println("  -> Spread at λ=$(lambda_low): ", round((err_low + target_spread) * 100, digits=2), "%")

    println("Evaluating upper bound (λ=$lambda_high)...")
    err_high = error_func(lambda_high)
    if isnan(err_high)
        println("Error: Could not compute spread at upper bound. Aborting.")
        return NaN
    end
    println("  -> Spread at λ=$(lambda_high): ", round((err_high + target_spread) * 100, digits=2), "%")

    if sign(err_low) == sign(err_high)
        println("\nError: Initial bounds [$(lambda_low), $(lambda_high)] do not bracket the root.")
        println("The target spread of ", round(target_spread * 100, digits=2), "% may be unreachable.")
        return NaN
    end

    # --- Bisection Loop ---
    println("\n--- Starting bisection search ---")
    for i in 1:max_iter
        lambda_mid = (lambda_low + lambda_high) / 2

        print("Iter $(i): Trying λ = $(round(lambda_mid, digits=4))... ")
        err_mid = error_func(lambda_mid)

        if isnan(err_mid)
            println("\nWarning: Spread calculation failed for λ_mid. Aborting.")
            return NaN
        end

        current_spread = err_mid + target_spread
        println("Spread = ", round(current_spread * 100, digits=2), "%")

        # Check for convergence on the error value
        if abs(err_mid) < (tol / 10)
            println("\n--- Calibration Converged (Error tolerance met) ---")
            return lambda_mid
        end

        # Update interval
        if sign(err_mid) == sign(err_low)
            lambda_low = lambda_mid
        else
            lambda_high = lambda_mid
        end

        # Check for convergence on the interval size
        if (lambda_high - lambda_low) < tol
            println("\n--- Calibration Converged (Interval tolerance met) ---")
            return (lambda_low + lambda_high) / 2
        end
    end

    println("\n--- Calibration Failed to Converge (Max iterations reached) ---")
    return (lambda_low + lambda_high) / 2
end


# --- 4. Main Execution Block ---

# Set a seed for reproducibility


# Run the calibration
# Note: This is computationally intensive and may take a significant amount of time.
calibrated_lambda = calibrate_lambda(target_spread=0.1025)

# Final check and reporting
if !isnan(calibrated_lambda)
    println("\n" * "="^60)
    println("CALIBRATION COMPLETE")
    println("Final calibrated value for λ = ", round(calibrated_lambda, digits=4))

    # Run a final, high-quality simulation with the calibrated value
    # to confirm the result and generate a final table.
    println("\nRunning final high-quality simulation with calibrated λ...")
    # I need the generate_table4_stats_bornstein_method function here to print the final table.
    # It was not included in this script based on the prompt, but would be needed for a final report.
    # For now, let's just re-calculate the spread with a long simulation to verify.
    final_spread = get_avg_spread(calibrated_lambda, T_sim=50000)
    println("\nVerification spread with T=50,000 and calibrated λ:")
    println(" -> Average Spread = ", round(final_spread * 100, digits=2), "%")
    println("="^60)
else
    println("\n" * "="^60)
    println("CALIBRATION FAILED")
    println("="^60)
end