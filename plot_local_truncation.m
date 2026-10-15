function plot_local_truncation
    % Plotting parameters
    h_list = logspace(-5, -1, 100);
    tspan = [0 10];
    X0 = 1;
    t_ex = 0.492;
    X_ex = solution02(t_ex);

    % Initialize variables for loop
    error_list_euler = zeros(size(h_list));
    global_truncation_euler = zeros(size(h_list));
    error_list_midpoint = zeros(size(h_list));
    global_truncation_midpoint = zeros(size(h_list));
    error_list_diff = zeros(size(h_list));
    h_avg_list = zeros(size(h_list));
    num_evals_midpoint = 0;
    num_evals_euler_list = zeros(size(h_list));
    num_evals_midpoint_list = zeros(size(h_list));

    % Find truncation for euler step
    % Cycle through h values
    for i = 1:length(h_list)
        [G, ~] = forward_euler_step(@rate_func02,t_ex,X_ex,h_list(i)); % run euler step
        X_euler = solution02(t_ex+h_list(i)); % calculate X(t+h)
        local_error = norm(G-X_euler); % calculate error for local truncation
        error_list_euler(i) = local_error;
        error_list_diff(i) = norm(X_euler-X_ex);

        [~,X_list,h_avg,num_evals_forward_euler] = ...
            forward_euler_fixed_step_integration(@rate_func02,tspan,solution02(tspan(1)),h_list(i));
        h_avg_list(i) = h_avg;
        num_evals_euler_list(i) = num_evals_forward_euler;
        global_truncation_euler(i) = norm(X_list(end,:)' - solution02(tspan(end)));

    end

    idx = h_list>= 1e-3 & h_list<=1e-1;
    mdl_local_euler = polyfit(log(h_list(idx)),log(error_list_euler(idx)),1); % calculate fit line
    mdl_global_euler = polyfit(log(h_list(idx)),log(global_truncation_euler(idx)),1);
    mdl_local_diff = polyfit(log(h_list(idx)), log(error_list_diff(idx)),1);
    p_local_euler = mdl_local_euler(1);

    % Find truncation for explicit midpoint
    % Cycle through h values
    for j = 1:length(h_list)
        [G,~] = explicit_midpoint_step(@rate_func02,t_ex,X_ex,h_list(j)); % run midpoint step
        X_midpoint = solution02(t_ex+h_list(j)); % calculate X(t+h)
        error_list_midpoint(j) = norm(G-X_midpoint);
        
        [~, X_list, ~, num_evals_mid] = explicit_midpoint_fixed_step_integration(@rate_func02,tspan,solution02(tspan(1)),h_list(j));
        global_truncation_midpoint(j) = norm(X_list(end,:)' - solution02(tspan(end)));
        num_evals_midpoint_list(j) = num_evals_mid;

    end

    mdl_local_midpoint = polyfit(log(h_list(idx)),log(error_list_midpoint(idx)),1); % calculate fit line
    mdl_global_midpoint = polyfit(log(h_list(idx)),log(global_truncation_midpoint(idx)),1); % calculate fit line
    p_local_midpoint = mdl_local_midpoint(1);
    p_local_diff = mdl_local_diff(1);

    mdl_evals_euler = polyfit(log(num_evals_euler_list(idx)), log(global_truncation_euler(idx)), 1);
    mdl_evals_midpoint = polyfit(log(num_evals_midpoint_list(idx)), log(global_truncation_midpoint(idx)), 1);
    p_evals_euler = mdl_evals_euler(1);
    p_evals_midpoint = mdl_evals_midpoint(1); 
    
    % Plot Euler method local error for all methods
    figure();
    loglog(h_list, error_list_euler, 'ro')
    hold on;
    loglog(h_list, error_list_midpoint, 'bs')
    loglog(h_list, error_list_diff, 'k^')
    loglog(h_list, exp(polyval(mdl_local_euler, log(h_list))), 'r-')      
    loglog(h_list, exp(polyval(mdl_local_midpoint, log(h_list))), 'b-')
    loglog(h_list, exp(polyval(mdl_local_diff, log(h_list))), 'k-')
    title('Local Truncation Error', 'Interpreter', 'Latex', 'FontSize', 22)
    xlabel('Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Local Truncation Error (-)', 'Interpreter', 'Latex', 'FontSize', 18)
    legend('Forward Euler', 'Explicit Midpoint', '$|X(t+h)-X(t)|$', 'FE fit', 'EM fit', 'Diff fit', ...
        'Interpreter', 'Latex', 'Location', 'northwest')
    hold off

    % Plot Euler method global error
    figure();
    loglog(h_list, global_truncation_euler, 'ro')
    hold on;
    title('Global Truncation Error vs. Step Size', 'Interpreter', 'Latex', 'FontSize', 22)
    xlabel('Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 18)
    ylabel('Global Truncation Error (-)', 'Interpreter', 'Latex', 'FontSize', 18)
    hold on;
    loglog(h_list, exp(polyval(mdl_global_euler, log(h_list))), 'r-')  
    loglog(h_list, global_truncation_midpoint, 'bo')
    loglog(h_list, exp(polyval(mdl_global_midpoint, log(h_list))), 'b-')
    legend('Forward Euler', 'Forward Euler Fit','Explicit Midpoint','Explicit Midpoint Fit', 'Interpreter', 'Latex', 'Location', 'northwest')
    hold off

    % % Plot Midpoint method global error
    % figure();
    % loglog(h_list, global_truncation_midpoint, 'ro')
    % title('Global truncation Error for Explicit Midpoint', 'Interpreter', 'Latex', 'FontSize', 20)
    % xlabel('Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 15)
    % ylabel('Error Size (-)', 'Interpreter', 'Latex', 'FontSize', 15)
    % hold on;
    % loglog(h_list, exp(polyval(mdl_global_midpoint, log(h_list))), 'b-')
    % legend('Sample Points', 'Fit Line', 'Interpreter', 'Latex')
    % hold off

%plot global truncation w/ rate functions
evals_fit_euler = num_evals_euler_list(idx);
evals_fit_midpoint =num_evals_midpoint_list(idx);
  figure();
  loglog(num_evals_euler_list, global_truncation_euler, 'ro', 'MarkerSize', 6)
  hold on;
  loglog(num_evals_midpoint_list, global_truncation_midpoint, 'bo', 'MarkerSize', 6)
  loglog(evals_fit_euler, exp(polyval(mdl_evals_euler, log(evals_fit_euler))), 'r-', 'LineWidth', 2)
  loglog(evals_fit_midpoint, exp(polyval(mdl_evals_midpoint, log(evals_fit_midpoint))), 'b-', 'LineWidth', 2)
  title('Global Truncation Error vs. Rate Function Calls', 'Interpreter', 'Latex', 'FontSize', 20)
  xlabel('Rate Function Calls (-)', 'Interpreter', 'Latex', 'FontSize', 18)
  ylabel('Global Truncation Error (-)', 'Interpreter', 'Latex', 'FontSize', 18)
  legend('Forward Euler', 'Explicit Midpoint', 'Forward Euler Fit', 'Explicit Method Fit','Interpreter', 'Latex', 'FontSize', 14, 'Location', 'northeast')
    hold off
end