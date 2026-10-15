% Create a figure that compares the numerical solutions to
% eqn. 5 with a target time step value of href = .38,
% generated using  forward Euler, backward Euler, 
% explicit midpoint, and implicit midpoint.

function plot_stability
    % initialize variables
    h_list = [0.38, 0.45];
    h_ref = h_list(2);
    tspan = [0 20];
    t_list = linspace(tspan(1), tspan(2), 500);
    x0 = solution01(tspan(1));
    x_exact = solution01(t_list);

    %forward euler
    [t_fe, X_fe,~,~] = forward_euler_fixed_step_integration( ...
    @rate_func01,tspan, x0, h_ref);

    % backward euler 
    [t_be, X_be,~,~] = implicit_fixed_step_integration( ...
    @rate_func01,@backward_euler_step, tspan, x0, h_ref);

    % explicit midpt
    [t_ex ,X_ex,~,~] = explicit_midpoint_fixed_step_integration( ...
    @rate_func01,tspan, x0, h_ref);

    % implicit midpt
    [t_im ,X_im,~,~] = implicit_fixed_step_integration( ...
    @rate_func01,@implicit_midpoint_step, tspan, x0, h_ref);

    figure;

    % Forward euler subplot
    subplot(4, 1, 1);
    plot(t_list, x_exact, 'b-');
    hold on;
    plot(t_fe, X_fe, 'r--');
    hold  off;
    ylabel('x(t) (-)', 'Interpreter', 'Latex', 'FontSize', 10);
    title('Forward Euler', ...
        'Interpreter', 'Latex', 'FontSize', 15);
    legend('Analytical', 'Numerical', ...
        'Interpreter', 'Latex', 'Location', 'East');    
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 10);
    

    % Backward euler subplot
    subplot(4, 1, 2);
    plot(t_list, x_exact, 'b-');
    hold on;
    plot(t_be, X_be, 'r--');
    hold  off;
    ylabel('x(t) (-)', 'Interpreter', 'Latex', 'FontSize', 10);
    title('Backward Euler', ...
        'Interpreter', 'Latex', 'FontSize', 15);
    legend('Analytical', 'Numerical', ...
        'Interpreter', 'Latex', 'Location', 'East');   
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 10);


    % Explicit midpoint subplot
    subplot(4, 1, 3);
    plot(t_list, x_exact, 'b-');
    hold on;
    plot(t_ex, X_ex, 'r--');
    hold  off;
    ylabel('x(t) (-)', 'Interpreter', 'Latex', 'FontSize', 10);
    title('Explicit Midpoint', ...
        'Interpreter', 'Latex', 'FontSize', 15);
    legend('Analytical', 'Numerical', ...
        'Interpreter', 'Latex', 'Location', 'East');  
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 10);


    % Explicit midpoint subplot
    subplot(4, 1, 4);
    plot(t_list, x_exact, 'b-');
    hold on;
    plot(t_im, X_im, 'r--');
    hold  off;
    ylabel('x(t) (-)', 'Interpreter', 'Latex', 'FontSize', 10);
    title('Implicit Midpoint', ...
        'Interpreter', 'Latex', 'FontSize', 15);
    legend('Analytical', 'Numerical', ...
        'Interpreter', 'Latex', 'Location', 'East');   
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 10);
    
    sgtitle('Integration Stability, $h_{ref}$ = 0.45', ...
    'interpreter', 'latex', 'fontsize', 20)
end