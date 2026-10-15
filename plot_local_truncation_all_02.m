function plot_local_truncation_all_02
    % Plotting parameters
    h_list = logspace(-4, -1, 20);
    tspan = [0 10];
    X0 = 1;
    t_ex = 0.492;
    X_ex = solution02(t_ex);

    % Initialize variables for loop
    error_list_euler_f = zeros(size(h_list));
    global_truncation_euler_f = zeros(size(h_list));
    error_list_midpoint_ex = zeros(size(h_list));
    global_truncation_midpoint_ex = zeros(size(h_list));

    error_list_euler_b = zeros(size(h_list));
    global_truncation_euler_b = zeros(size(h_list));
    error_list_midpoint_im = zeros(size(h_list));
    global_truncation_midpoint_im = zeros(size(h_list));

    error_list_diff = zeros(size(h_list));
    h_avg_list = zeros(size(h_list));

    num_evals_midpoint = 0;

    num_evals_euler_list_f = zeros(size(h_list));
    num_evals_midpoint_list_ex = zeros(size(h_list));
    num_evals_euler_list_b = zeros(size(h_list));
    num_evals_midpoint_list_im = zeros(size(h_list));

    for i = 1:length(h_list)
        h = h_list(i);

        %forward euler
        [G, ~] = forward_euler_step(@rate_func02,t_ex,X_ex,h);
        X_euler = solution02(t_ex+h);
        error_list_euler_f(i) = norm(G-X_euler);
        error_list_diff(i) = norm(X_euler-X_ex);


        [~,X_list,h_avg,num_evals] = forward_euler_fixed_step_integration(@rate_func02,tspan,solution02(tspan(1)),h);
        h_avg_list(i) = h_avg;
        num_evals_euler_list_f(i) = num_evals;
        global_truncation_euler_f(i) = norm(X_list(end,:)' - solution02(tspan(end)));

        % backward euler 
        [G, ~] = backward_euler_step(@rate_func02,t_ex,X_ex,h);
        error_list_euler_b(i) = norm(G-X_euler);

        [~,X_list,~,num_evals] = implicit_fixed_step_integration(@rate_func02,@backward_euler_step, tspan,solution02(tspan(1)),h);
        num_evals_euler_list_b(i) = num_evals;
        global_truncation_euler_b(i) = norm(X_list(end,:)' - solution02(tspan(end)));

        % explicit midpt
        [G, ~] = explicit_midpoint_step(@rate_func02,t_ex,X_ex,h);
        error_list_midpoint_ex(i) = norm(G-X_euler);

        [~,X_list,~,num_evals] = explicit_midpoint_fixed_step_integration(@rate_func02,tspan,solution02(tspan(1)),h);
        num_evals_midpoint_list_ex(i) = num_evals;
        global_truncation_midpoint_ex(i) = norm(X_list(end,:)' - solution02(tspan(end)));

        % implicit midpt
        [G, ~] = implicit_midpoint_step(@rate_func02,t_ex,X_ex,h);
        error_list_midpoint_im(i) = norm(G-X_euler);

        [~,X_list,~,num_evals] = implicit_fixed_step_integration(@rate_func02,@implicit_midpoint_step, tspan,solution02(tspan(1)),h);
        num_evals_midpoint_list_im(i) = num_evals;
        global_truncation_midpoint_im(i) = norm(X_list(end,:)' - solution02(tspan(end)));
    end

    % Find local p values
    idx = h_list>= 1e-3 & h_list<=1e-1;
    mdl_local_euler_f = polyfit(log(h_list(idx)), log(error_list_euler_f(idx)),1);
    mdl_local_euler_b = polyfit(log(h_list(idx)), log(error_list_euler_b(idx)),1);
    mdl_local_midpoint_ex = polyfit(log(h_list(idx)), log(error_list_midpoint_ex(idx)),1);
    mdl_local_midpoint_im = polyfit(log(h_list(idx)), log(error_list_midpoint_im(idx)),1);
    mdl_local_diff = polyfit(log(h_list(idx)), log(error_list_diff(idx)),1);

    p_local_euler_f = mdl_local_euler_f(1)
    p_local_euler_b = mdl_local_euler_b(1)
    p_local_midpoint_ex = mdl_local_midpoint_ex(1)
    p_local_midpoint_im = mdl_local_midpoint_im(1)
    p_local_diff = mdl_local_diff(1)

    % Find global p values
    mdl_global_euler_f = polyfit(log(h_list(idx)), log(global_truncation_euler_f(idx)),1);
    mdl_global_euler_b = polyfit(log(h_list(idx)), log(global_truncation_euler_b(idx)),1);
    mdl_global_midpoint_ex = polyfit(log(h_list(idx)), log(global_truncation_midpoint_ex(idx)),1);
    mdl_global_midpoint_im = polyfit(log(h_list(idx)), log(global_truncation_midpoint_im(idx)),1);

    p_global_euler_f = mdl_global_euler_f(1)
    p_global_euler_b = mdl_global_euler_b(1)
    p_global_midpoint_ex = mdl_global_midpoint_ex(1)
    p_global_midpoint_im = mdl_global_midpoint_im(1)

    
    % Plot local truncation error for all methods
    figure();
    loglog(h_list, error_list_euler_f, 'ro')
    hold on;
    loglog(h_list, error_list_euler_b, 'ms')
    loglog(h_list, error_list_midpoint_ex, 'bs')
    loglog(h_list, error_list_midpoint_im, 'gd')

    title('Local Truncation Error', ...
        'Interpreter', 'Latex', 'FontSize', 22)
    xlabel('Step Size: h (-)', ...
        'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Local Truncation Error (-)', ...
        'Interpreter', 'Latex', 'FontSize', 18)

    legend('Forward Euler', 'Backward Euler', 'Explicit Midpoint', 'Implicit Midpoint', ...
        'Interpreter', 'Latex', 'Location', 'northwest')
    hold off

    % Plot global error vs step size
    figure();
    loglog(h_list, global_truncation_euler_f, 'ro')
    hold on;
    loglog(h_list, global_truncation_euler_b, 'ms')
    loglog(h_list, global_truncation_midpoint_ex, 'bs')
    loglog(h_list, global_truncation_midpoint_im, 'gd')

    title('Global Error vs Step Size', ...
        'Interpreter', 'Latex', 'FontSize', 22)
    xlabel('Step Size: h (-)', ...
        'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Global Error (-)', ...
        'Interpreter', 'Latex', 'FontSize', 18)

    legend('Forward Euler', 'Backward Euler', ...
        'Explicit Midpoint', 'Implicit Midpoint', ...
        'Interpreter', 'Latex', 'Location', 'northwest')
    hold off


    % Plot global error vs rate function calls
    figure();
    loglog(num_evals_euler_list_f, global_truncation_euler_f, 'ro')
    hold on;
    loglog(num_evals_euler_list_b, global_truncation_euler_b, 'ms')
    loglog(num_evals_midpoint_list_ex, global_truncation_midpoint_ex, 'bs')
    loglog(num_evals_midpoint_list_im, global_truncation_midpoint_im, 'gd')

    title('Global Error vs Rate Function Calls', ...
        'Interpreter', 'Latex', 'FontSize', 22)
    xlabel('Number of Rate Function Calls (-)', ...
        'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Global Error (-)', ...
        'Interpreter', 'Latex', 'FontSize', 18)

    legend('Forward Euler', 'Backward Euler', ...
        'Explicit Midpoint', 'Implicit Midpoint', ...
        'Interpreter', 'Latex', 'Location', 'northwest')
    hold off


end