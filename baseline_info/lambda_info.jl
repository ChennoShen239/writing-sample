using LinearAlgebra, Statistics
using LaTeXStrings, QuantEcon, DataFrames, Plots, Random
using Interpolations
using Printf


mkpath("Result/Figs/lambda_info")

function ArellanoEconomy(;
    beta=0.953,
    gamma=2.0,
    r=0.017,
    rho=0.945,    # persistence in output
    eta=0.025,
    theta=0.282,
    ny=201,
    nB=251,
    λ=1.0,        # fraction of lenders who know how to compute the default probability
)

    # create grids
    Bgrid = collect(range(-0.4, 0.4, length=nB))
    mc = tauchen(ny, rho, eta)
    Pi = mc.p
    ygrid = exp.(mc.state_values)
    ydefgrid = min.(0.969 * mean(ygrid), ygrid)

    # define value functions
    # notice ordered different than Python to take
    # advantage of column major layout of Julia
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

    # unpack stuff
    (; beta, gamma, r, rho, eta, theta, ny, nB) = ae
    (; ygrid, ydefgrid, Bgrid, Pi, vf, vd, vc, policy, q, defprob) = ae
    zero_ind = searchsortedfirst(Bgrid, 0.0)

    for iy = 1:ny
        y = ae.ygrid[iy]
        ydef = ae.ydefgrid[iy]

        # value of being in default with income y
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

            # update value and policy functions
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
    # unpack parameters
    (; beta, gamma, r, rho, eta, theta, ny, nB) = ae

    # create default values with a matching size
    vd_compact = repeat(ae.vd, nB)
    default_states = vd_compact .> ae.vc
    def_interp = interp2d(ae.Bgrid, ae.ygrid, vd_compact)
    vd_exp_y = [def_interp(B, y^ae.rho) for B in ae.Bgrid, y in ae.ygrid]

    repay_interp = interp2d(ae.Bgrid, ae.ygrid, ae.vc)
    repay_exp_y = [repay_interp(B, y^ae.rho) for B in ae.Bgrid, y in ae.ygrid]

    def_states_exp_y = vd_exp_y .> repay_exp_y

    # only λ fraction of lenders know how to compute the default probability
    # the other 1 - λ fraction of lenders will only look at the default state of the expected next period y 
    # update default probabilities and prices
    copyto!(ae.defprob, ae.λ * default_states * ae.Pi' + (1 - ae.λ) * def_states_exp_y)
    copyto!(ae.q, (1 .- ae.defprob) / (1 + r))
    return
end

function vfi!(ae; tol=1e-8, maxit=10000)

    # unpack stuff
    (; beta, gamma, r, rho, eta, theta, ny, nB) = ae
    (; ygrid, ydefgrid, Bgrid, Pi, vf, vd, vc, policy, q, defprob) = ae
    Pit = Pi'

    # Iteration stuff
    it = 0
    dist = 10.0

    # allocate memory for update
    V_upd = similar(ae.vf)

    while dist > tol && it < maxit
        it += 1

        # compute expectations for this iterations
        # (we need Pi' because of order value function dimensions)
        copyto!(V_upd, ae.vf)
        EV = ae.vf * Pit
        EVd = ae.vd * Pit
        EVc = ae.vc * Pit

        # update value function
        one_step_update!(ae, EV, EVd, EVc)

        # update prices
        compute_prices!(ae)

        dist = maximum(abs(x - y) for (x, y) in zip(V_upd, ae.vf))

        if it % 25 == 0
            println("Finished iteration $(it) with dist of $(dist)")
        end
    end
end

function AEsimulate(ae, capT=5000; y_init=mean(ae.ygrid), B_init=mean(ae.Bgrid))

    # get initial indices
    zero_index = searchsortedfirst(ae.Bgrid, 0.0)
    y_init_ind = searchsortedfirst(ae.ygrid, y_init)
    B_init_ind = searchsortedfirst(ae.Bgrid, B_init)

    # create a QE MarkovChain
    mc = MarkovChain(ae.Pi)
    y_sim_indices = simulate(mc, capT + 1; init=y_init_ind)

    # allocate and fill output
    y_sim_val = zeros(capT + 1)
    B_sim_val, q_sim_val = similar(y_sim_val), similar(y_sim_val)
    B_sim_indices = fill(0, capT + 1)
    default_status = fill(false, capT + 1)
    B_sim_indices[1], default_status[1] = B_init_ind, false
    y_sim_val[1], B_sim_val[1] = ae.ygrid[y_init_ind], ae.Bgrid[B_init_ind]

    for t = 1:capT
        # get today's indexes
        yi, Bi = y_sim_indices[t], B_sim_indices[t]
        defstat = default_status[t]

        # if you are not in default
        if !defstat
            default_today = ae.vc[Bi, yi] < ae.vd[yi]

            if default_today
                # default values
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

            # if you are in default
        else
            B_sim_indices[t+1] = zero_index
            B_sim_val[t+1] = 0.0
            y_sim_val[t] = ae.ydefgrid[y_sim_indices[t]]
            q_sim_val[t] = ae.q[zero_index, y_sim_indices[t]]

            # with probability theta exit default status
            default_status[t+1] = rand() >= ae.theta
        end
    end

    return (y_sim_val[1:capT], B_sim_val[1:capT], q_sim_val[1:capT], default_status[1:capT])
end

ae_λ = ArellanoEconomy(
    beta=0.953,      # time discount rate
    gamma=2.0,       # risk aversion
    r=0.017,        # international interest rate
    rho=0.945,       # persistence in output
    eta=0.025,      # st dev of output shock
    theta=0.282,    # prob of regaining access
    ny=201,          # number of points in y grid
    nB=401,
    λ=0.5179687500000001,
)         # number of points in B grid

# now solve the model on the grid.
@time vfi!(ae_λ)

ae_full_info = ArellanoEconomy(
    beta=0.953,      # time discount rate
    gamma=2.0,       # risk aversion
    r=0.017,        # international interest rate
    rho=0.945,       # persistence in output
    eta=0.025,      # st dev of output shock
    theta=0.282,    # prob of regaining access
    ny=201,          # number of points in y grid
    nB=401,
    λ=1.0,
)

@time vfi!(ae_full_info)

function plot_lambda_info()
    gr()
    default(fontfamily="Times New Roman")
    default(background_color=:transparent)
    default(background_color_legend=:transparent)
    default(grid=false)
    # create "Y High" and "Y Low" values as 5% devs from mean
    high, low = 1.05 * mean(ae_λ.ygrid), 0.95 * mean(ae_λ.ygrid)
    iy_high, iy_low = map(x -> searchsortedfirst(ae_λ.ygrid, x), (high, low))

    # extract a suitable plot grid
    x = zeros(0)
    q_low_λ = zeros(0)
    q_high_λ = zeros(0)
    for i = 1:(ae_λ.nB)
        b = ae_λ.Bgrid[i]
        if -0.35 <= b <= 0  # to match fig 3 of Arellano
            push!(x, b)
            push!(q_low_λ, ae_λ.q[i, iy_low])
            push!(q_high_λ, ae_λ.q[i, iy_high])
        end
    end

    q_low_full_info = zeros(0)
    q_high_full_info = zeros(0)
    for i = 1:(ae_full_info.nB)
        b = ae_full_info.Bgrid[i]
        if -0.35 <= b <= 0  # to match fig 3 of Arellano
            push!(q_low_full_info, ae_full_info.q[i, iy_low])
            push!(q_high_full_info, ae_full_info.q[i, iy_high])
        end
    end

    # generate plot and save to pdf
    p1 = plot(x, q_low_λ, label=L"y_{Low} (λ = 0.518)", linewidth=2)
    plot!(p1, x, q_high_λ, label=L"y_{High} (λ = 0.518)", linewidth=2)
    plot!(p1, x, q_low_full_info, label=L"y_{Low} (λ = 1.0)", linewidth=2)
    plot!(p1, x, q_high_full_info, label=L"y_{High} (λ = 1.0)", linewidth=2)
    plot!(
        p1,
        # title=L"Bond price schedule $q(y, B^\prime)$",
        xlabel=L"B^\prime",
        ylabel=L"q",
        # legend = :topleft,
        legend_frame=false,
    )

    savefig(p1, "Result/Figs/lambda_info/bond_price_schedule_$(ae_λ.λ).pdf")

    # Value functions plot
    p2 = plot(ae_λ.Bgrid, ae_λ.vf[:, iy_low], label=L"y_{Low} (λ = 0.518)", linewidth=2)
    plot!(p2, ae_λ.Bgrid, ae_λ.vf[:, iy_high], label=L"y_{High} (λ = 0.518)", linewidth=2)
    plot!(p2, ae_full_info.Bgrid, ae_full_info.vf[:, iy_low], label=L"y_{Low} (λ = 1.0)")
    plot!(
        p2,
        ae_full_info.Bgrid,
        ae_full_info.vf[:, iy_high],
        label=L"y_{High} (λ = 1.0)",
    )
    plot!(p2, xlabel=L"B", ylabel=L"V(y,B)")
    savefig(p2, "Result/Figs/lambda_info/value_functions_$(ae_λ.λ).pdf")

    # Ensure 'ae' is the solved ArellanoEconomy instance from your previous code.
    # Ensure iy_high and iy_low are defined as in your previous code:
    # high, low = 1.05 * mean(ae.ygrid), 0.95 * mean(ae.ygrid)
    # iy_high, iy_low = map(x -> searchsortedfirst(ae.ygrid, x), (high, low))

    # Calculate interest rates for low income
    rates_low_income_λ = zeros(ae_λ.nB)
    for ib = 1:ae_λ.nB
        # Get the optimal policy for B_prime (index)
        b_prime_idx = Int(ae_λ.policy[ib, iy_low])
        # Get the corresponding bond price q(B_prime, y_low)
        q_val = ae_λ.q[b_prime_idx, iy_low]
        # Calculate interest rate, handling potential q_val = 0
        if q_val > 1e-12 # Avoid division by zero or extremely small q
            rates_low_income_λ[ib] = (1 / q_val) - 1
        else
            rates_low_income_λ[ib] = NaN # Or some large number if preferred
        end
    end

    # Calculate interest rates for high income
    rates_high_income_λ = zeros(ae_λ.nB)
    for ib = 1:ae_λ.nB
        # Get the optimal policy for B_prime (index)
        b_prime_idx = Int(ae_λ.policy[ib, iy_high])
        # Get the corresponding bond price q(B_prime, y_high)
        q_val = ae_λ.q[b_prime_idx, iy_high]
        # Calculate interest rate
        if q_val > 1e-12
            rates_high_income_λ[ib] = (1 / q_val) - 1
        else
            rates_high_income_λ[ib] = NaN
        end
    end

    rates_low_income_full_info = zeros(ae_full_info.nB)
    for ib = 1:ae_full_info.nB
        b_prime_idx = Int(ae_full_info.policy[ib, iy_low])
        q_val = ae_full_info.q[b_prime_idx, iy_low]
        if q_val > 1e-12
            rates_low_income_full_info[ib] = (1 / q_val) - 1
        else
            rates_low_income_full_info[ib] = NaN
        end
    end

    rates_high_income_full_info = zeros(ae_full_info.nB)
    for ib = 1:ae_full_info.nB
        b_prime_idx = Int(ae_full_info.policy[ib, iy_high])
        q_val = ae_full_info.q[b_prime_idx, iy_high]
        if q_val > 1e-12
            rates_high_income_full_info[ib] = (1 / q_val) - 1
        else
            rates_high_income_full_info[ib] = NaN
        end
    end

    # Generate plot
    p_interest_rate =
        plot(ae_λ.Bgrid, rates_low_income_λ, label=L"y_{Low} (λ = 0.518)", linewidth=2)
    plot!(
        p_interest_rate,
        ae_λ.Bgrid,
        rates_high_income_λ,
        label=L"y_{High} (λ = 0.518)",
        linewidth=2,
    )
    plot!(
        p_interest_rate,
        ae_full_info.Bgrid,
        rates_low_income_full_info,
        label=L"y_{Low} (λ = 1.0)",
        linewidth=2,
    )
    plot!(
        p_interest_rate,
        ae_full_info.Bgrid,
        rates_high_income_full_info,
        label=L"y_{High} (λ = 1.0)",
        linewidth=2,
    )
    plot!(
        p_interest_rate,
        # title=L"Equilibrium Interest Rate $1/q(B'(B,y),y) - 1$",
        xlabel=L"$B$",
        ylabel=L"$1/q - 1$",
        # legend = :topright, # Adjust legend position as needed
        xlims=(-0.1, 0.03), # Set x-axis limits
        ylims=(0, 0.3), # Set y-axis limits
    )

    # Save the plot
    savefig(
        p_interest_rate,
        "Result/Figs/lambda_info/equilibrium_interest_rate_$(ae_λ.λ).pdf",
    )

    # Display the plot (optional, depending on your environment)
    # display(p_interest_rate)

    # Ensure 'ae' is the solved ArellanoEconomy instance from your previous code.
    # Ensure iy_high and iy_low are defined as in your previous code:
    # high, low = 1.05 * mean(ae.ygrid), 0.95 * mean(ae.ygrid)
    # iy_high, iy_low = map(x -> searchsortedfirst(ae.ygrid, x), (high, low))

    # Get savings (B') for low income
    savings_low_income_λ = zeros(ae_λ.nB)
    for ib = 1:ae_λ.nB
        # Get the optimal policy for B_prime (index)
        b_prime_idx = Int(ae_λ.policy[ib, iy_low])
        # Convert index to B_prime value
        savings_low_income_λ[ib] = ae_λ.Bgrid[b_prime_idx]
    end

    # Get savings (B') for high income
    savings_high_income_λ = zeros(ae_λ.nB)
    for ib = 1:ae_λ.nB
        # Get the optimal policy for B_prime (index)
        b_prime_idx = Int(ae_λ.policy[ib, iy_high])
        # Convert index to B_prime value
        savings_high_income_λ[ib] = ae_λ.Bgrid[b_prime_idx]
    end

    savings_low_income_full_info = zeros(ae_full_info.nB)
    for ib = 1:ae_full_info.nB
        b_prime_idx = Int(ae_full_info.policy[ib, iy_low])
        savings_low_income_full_info[ib] = ae_full_info.Bgrid[b_prime_idx]
    end

    savings_high_income_full_info = zeros(ae_full_info.nB)
    for ib = 1:ae_full_info.nB
        b_prime_idx = Int(ae_full_info.policy[ib, iy_high])
        savings_high_income_full_info[ib] = ae_full_info.Bgrid[b_prime_idx]
    end

    # Generate plot
    p_savings =
        plot(ae_λ.Bgrid, savings_low_income_λ, label=L"y_{Low} (λ = 0.518)", linewidth=2)
    plot!(
        p_savings,
        ae_λ.Bgrid,
        savings_high_income_λ,
        label=L"y_{High} (λ = 0.518)",
        linewidth=2,
    )
    plot!(
        p_savings,
        ae_full_info.Bgrid,
        savings_low_income_full_info,
        label=L"y_{Low} (λ = 1.0)",
        linewidth=2,
    )
    plot!(
        p_savings,
        ae_full_info.Bgrid,
        savings_high_income_full_info,
        label=L"y_{High} (λ = 1.0)",
        linewidth=2,
    )

    # Add a 45-degree line for reference (B' = B)
    # This helps visualize where assets are increasing/decreasing
    plot!(
        p_savings,
        ae_λ.Bgrid,
        ae_λ.Bgrid,
        label=L"B' = B",
        linestyle=:dash,
        color=:grey,
    )

    plot!(
        p_savings,
        #  title=L"Savings Function $B'(B,y)$",
        xlabel=L"$B$",
        ylabel=L"$B'$",
        # legend = :bottomright, # Adjust legend position as needed
        xlims=(-0.3, 0.15),
        ylims=(-0.2, 0.1),
    )

    # Save the plot
    savefig(p_savings, "Result/Figs/lambda_info/savings_function_$(ae_λ.λ).pdf")
end

plot_lambda_info()


function generate_table4_stats_bornstein_method(
    ae,
    y_vec_orig,
    B_vec_orig,
    q_vec_orig,
    default_vec_orig,
)

    # --- 0. Initial Checks and Burn-in Reminder ---
    if length(y_vec_orig) < 2
        println("Warning: Input vectors are too short for stats calculation. Aborting.")
        return
    end
    println("Reminder: Ensure input vectors are post-burn-in for accurate stats.")
    T_orig = length(y_vec_orig)

    # --- 1. Define Core Economic Variables (Length T_orig - 1 generally) ---
    # Ensure q_vec_orig is at least T_orig-1 for B_{t+1}
    if length(q_vec_orig) < T_orig - 1 || length(B_vec_orig) < T_orig
        println("Warning: q_vec or B_vec too short for consumption calculation. Aborting.")
        return
    end

    # Consumption: c_t = y_t + B_t - q_t * B_{t+1}
    # B_vec_orig[2:T_orig] is B_{t+1}
    # q_vec_orig[1:(T_orig-1)] is q_t
    # B_vec_orig[1:(T_orig-1)] is B_t
    # y_vec_orig[1:(T_orig-1)] is y_t
    c_vec =
        y_vec_orig[1:(T_orig-1)] .+ B_vec_orig[1:(T_orig-1)] .-
        q_vec_orig[1:(T_orig-1)] .* B_vec_orig[2:T_orig]

    y_trimmed = y_vec_orig[1:(T_orig-1)]
    B_trimmed = B_vec_orig[1:(T_orig-1)] # This is B_t
    q_trimmed = max.(q_vec_orig[1:(T_orig-1)], 1e-12) # q_t
    default_trimmed = default_vec_orig[1:(T_orig-1)]

    # Quarterly spread: s_q = 1/q_t - (1+r)
    # quarterly_spread_trimmed = (1.0 ./ q_trimmed) .- (1 .+ ae.r)
    # Annualized percentage spread: ((1+s_q)^4 - 1) * 100
    # annual_pct_spread_trimmed = ((1 .+ quarterly_spread_trimmed) .^ 4 .- 1.0) # Keep as fraction for now, scale at print
    annual_pct_spread_trimmed = (1.0 ./ q_trimmed) .^ 4 .- (1 .+ ae.r) .^ 4

    trade_balance_ratio_trimmed = zeros(T_orig - 1)
    for i = 1:length(y_trimmed)
        if abs(y_trimmed[i]) > 1e-9
            trade_balance_ratio_trimmed[i] = (y_trimmed[i] - c_vec[i]) / y_trimmed[i]
        else
            trade_balance_ratio_trimmed[i] = NaN
        end
    end

    # --- 2. Sample Windows & Calculate Per-Window Statistics ---
    window_size = 74 # Quarters preceding a default event
    max_events_to_sample = 100 # Or as desired
    default_event_indices = findall(default_trimmed .== true) # Indices where default_trimmed is true

    # Lists to store statistics from each valid window
    window_means_Spread_val = Float64[] # store annualized spread value (fraction)
    window_means_TB_ratio_val = Float64[]

    window_stds_Y_log = Float64[]
    window_stds_C_log = Float64[]
    window_stds_TB_ratio = Float64[]
    window_stds_Spread_val = Float64[] # store std of annualized spread value (fraction)

    window_corrs_YC_log = Float64[]
    window_corrs_Yspread_log_vs_raw = Float64[]
    window_corrs_Cspread_log_vs_raw = Float64[]
    window_corrs_TBspread_ratio_vs_raw = Float64[]
    window_corrs_YTB_log_vs_ratio = Float64[]

    window_mean_debt_to_output_ratio = Float64[] # For -B/Y

    processed_event_count = 0
    last_sampled_window_end_idx = 0

    for event_idx in default_event_indices
        if processed_event_count >= max_events_to_sample
            break
        end
        # Window is [event_idx - window_size + 1, event_idx]
        # Bornstein uses [def_period - 74, def_period - 1]
        # So, if event_idx is the first period OF default, window ends at event_idx - 1
        window_end_actual = event_idx - 1
        window_start_idx = window_end_actual - window_size + 1


        if window_start_idx > 0 && window_start_idx > last_sampled_window_end_idx
            is_window_internally_clean = true
            # Check for defaults within the [window_start_idx, window_end_actual] window
            if window_end_actual >= window_start_idx # ensure window has some length
                if any(default_trimmed[window_start_idx:window_end_actual])
                    is_window_internally_clean = false
                end
            else # window has zero or negative length
                is_window_internally_clean = false
            end


            if is_window_internally_clean
                # Extract data for this window
                y_w = y_trimmed[window_start_idx:window_end_actual]
                c_w = c_vec[window_start_idx:window_end_actual]
                tb_ratio_w = trade_balance_ratio_trimmed[window_start_idx:window_end_actual]
                spread_annual_val_w =
                    annual_pct_spread_trimmed[window_start_idx:window_end_actual] # fraction
                B_w = B_trimmed[window_start_idx:window_end_actual] # Assets B_t

                if length(y_w) < 2 # Need at least 2 data points for std/corr
                    continue
                end

                # Log transform Y and C for this window
                log_y_w = log.(max.(y_w, 1e-9))
                log_c_w = log.(max.(c_w, 1e-9))

                # Store means (for Spread and TB/Y only, as per Bornstein's table structure)
                push!(window_means_Spread_val, mean(filter(isfinite, spread_annual_val_w)))
                push!(window_means_TB_ratio_val, mean(filter(isfinite, tb_ratio_w)))

                # Store std devs (log for Y, C; raw for TB, Spread)
                push!(window_stds_Y_log, std(filter(isfinite, log_y_w)))
                push!(window_stds_C_log, std(filter(isfinite, log_c_w)))
                push!(window_stds_TB_ratio, std(filter(isfinite, tb_ratio_w)))
                push!(window_stds_Spread_val, std(filter(isfinite, spread_annual_val_w)))

                # Helper for safe correlation calculation
                function safe_cor_window(v1, v2)
                    f1 = filter(isfinite, v1)
                    f2 = filter(isfinite, v2)
                    min_len = min(length(f1), length(f2))
                    if min_len < 2
                        return NaN
                    end
                    xf1 = f1[1:min_len]
                    xf2 = f2[1:min_len]
                    # Check for zero standard deviation to prevent errors/NaNs from cor()
                    if std(xf1) < 1e-9 || std(xf2) < 1e-9
                        return NaN
                    end
                    return cor(xf1, xf2)
                end

                push!(window_corrs_YC_log, safe_cor_window(log_y_w, log_c_w))
                push!(
                    window_corrs_Yspread_log_vs_raw,
                    safe_cor_window(log_y_w, spread_annual_val_w),
                )
                push!(
                    window_corrs_Cspread_log_vs_raw,
                    safe_cor_window(log_c_w, spread_annual_val_w),
                )
                push!(
                    window_corrs_TBspread_ratio_vs_raw,
                    safe_cor_window(tb_ratio_w, spread_annual_val_w),
                )
                push!(window_corrs_YTB_log_vs_ratio, safe_cor_window(log_y_w, tb_ratio_w))

                # Mean Debt to Output for this window (-B_t / Y_t)
                current_window_debt_to_output_ratios = -B_w ./ y_w
                push!(
                    window_mean_debt_to_output_ratio,
                    mean(filter(isfinite, current_window_debt_to_output_ratios)),
                )

                last_sampled_window_end_idx = window_end_actual
                processed_event_count += 1
            end
        end
    end

    if processed_event_count == 0
        println(
            "Warning: No valid default windows sampled. Cannot compute average statistics.",
        )
        return
    end

    # --- 3. Calculate Final Table Statistics (Averages of Per-Window Stats) ---
    safe_mean(stat_list) =
        isempty(stat_list) || all(isnan.(stat_list)) ? NaN :
        mean(filter(isfinite, stat_list))

    # Bornstein's table structure for correlations:
    # corr(c,y) -> col 3 for Consumption
    # corr(TB/Y,y) -> col 3 for Trade Balance
    # corr(spread,y) -> col 3 for Interest Rate Spread (and col 4 for Output)
    # corr(spread,TB/Y) -> col 4 for Trade Balance
    # corr(spread,c) -> col 4 for Consumption

    # Values for the table (Mean and Std will be scaled by 100 later)
    mean_spread = safe_mean(window_means_Spread_val)
    std_spread = safe_mean(window_stds_Spread_val)
    corr_spread_y = safe_mean(window_corrs_Yspread_log_vs_raw) # corr(spread, log Y)

    mean_tb = safe_mean(window_means_TB_ratio_val)
    std_tb = safe_mean(window_stds_TB_ratio)
    corr_tb_y = safe_mean(window_corrs_YTB_log_vs_ratio)    # corr(TB/Y, log Y)
    corr_tb_spread = safe_mean(window_corrs_TBspread_ratio_vs_raw) # corr(TB/Y, Spread)

    # Mean for C and Y are NaN as per Bornstein's typical table output
    std_c = safe_mean(window_stds_C_log)
    corr_c_y = safe_mean(window_corrs_YC_log)              # corr(log C, log Y)
    corr_c_spread = safe_mean(window_corrs_Cspread_log_vs_raw) # corr(log C, Spread)

    std_y = safe_mean(window_stds_Y_log)
    # corr_y_spread is corr_spread_y

    final_stats = [
        ("Interest rate spread", mean_spread, std_spread, corr_spread_y, NaN),
        ("Trade balance (/Y)", mean_tb, std_tb, corr_tb_y, corr_tb_spread),
        ("Consumption", NaN, std_c, corr_c_y, corr_c_spread),
        ("Output", NaN, std_y, NaN, corr_spread_y),
    ]

    # --- 4. Calculate "Other Statistics" ---
    total_sim_periods_default_prob = length(default_vec_orig)
    overall_default_prob = sum(default_vec_orig .== true) / total_sim_periods_default_prob # as fraction

    mean_debt_over_y_all_windows = safe_mean(window_mean_debt_to_output_ratio) # as fraction

    # --- 5. Print Results ---
    println("\nTable: Business Cycle Statistics (Bornstein's Method)\n")
    @printf(
        "%-28s %-18s %-10s %-15s %-15s\n",
        "Variable",
        "Mean (Def.Epis.)",
        "Std(x)",
        "corr(x,y)",
        "corr(x,spread)"
    )
    for (name, val_mean, val_std, val_corr_xy, val_corr_xspread) in final_stats
        # Scale means and stds by 100 for printing, except for NaNs
        mean_print = isnan(val_mean) ? val_mean : val_mean * 100
        std_print = isnan(val_std) ? val_std : val_std * 100
        @printf(
            "%-28s %-18.2f %-10.2f %-15.2f %-15.2f\n",
            name,
            mean_print,
            std_print,
            val_corr_xy,
            val_corr_xspread
        )
    end

    println("\nOther statistics (Attempting to match Bornstein's style):")
    @printf(
        "Mean debt (%% of output): %.2f\n",
        isnan(mean_debt_over_y_all_windows) ? NaN : mean_debt_over_y_all_windows * 100
    )
    # Mean spread is already in the table, can be re-iterated if desired.
    # mean_spread_overall_print = isnan(mean_spread) ? NaN : mean_spread * 100
    # @printf("Mean spread (annualized %%): %.2f\n", mean_spread_overall_print)
    @printf("Default probability (overall %%): %.2f\n", overall_default_prob * 100)
    @printf("Number of valid default events sampled: %d\n", processed_event_count)
    println("Note: 'Std(x)' and some 'Mean' values are scaled by 100.")
end


# Example of how you might call it (assuming ae, y_vec etc. are defined):
# generate_table4_stats_bornstein_method(ae, y_vec, B_vec, q_vec, default_vec)
T = 50000
Random.seed!(1234)
y_vec_λ, B_vec_λ, q_vec_λ, default_vec_λ = AEsimulate(ae_λ, T)
Random.seed!(1234)
y_vec_full_info, B_vec_full_info, q_vec_full_info, default_vec_full_info =
    AEsimulate(ae_full_info, T)
generate_table4_stats_bornstein_method(ae_λ, y_vec_λ, B_vec_λ, q_vec_λ, default_vec_λ)

generate_table4_stats_bornstein_method(
    ae_full_info,
    y_vec_full_info,
    B_vec_full_info,
    q_vec_full_info,
    default_vec_full_info,
)
