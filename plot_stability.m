% Create a figure that compares the numerical solutions to
% eqn. 5 with a target time step value of href = .38,
% generated using  forward Euler, backward Euler, 
% explicit midpoint, and implicit midpoint.

function plot_stability
    % initialize variables
    h_list = [0.38, 0.45];
    h_ref = h_list(1);
    tspan = [0 20];
    t_list = linspace(tspan(1), tspan(2), 200);
    x0 = solution01(tspan(1));
    x_exact = solution01(t_list);
    h_ref_title = num2str(h_ref);


    myfig = gcf(); %set myfig to the current figure
    exportgraphics(myfig , 'stability_38.png', 'Resolution', 600); %saves plot to plot_name at high resolution

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

    % Forward euler subplot
    subplot(4, 1, 1);
    plot(t_list, x_exact, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_fe, X_fe, 'r.-', 'MarkerSize', 10);
    hold off;
    ylabel('$x(t)$ (-)', 'Interpreter', 'Latex', 'FontSize', 15);
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 15);
    set(gca, 'FontSize', 12);
    set(gca,'TickLabelInterpreter','latex');
    title('Forward Euler', ...
        'Interpreter', 'Latex', 'FontSize', 20);

    % Backward euler subplot
    subplot(4, 1, 2);
    plot(t_list, x_exact, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_be, X_be, 'r.-', 'MarkerSize', 10);
    hold  off;
    ylabel('$x(t)$ (-)', 'Interpreter', 'Latex', 'FontSize', 20);
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 20);
    set(gca, 'FontSize', 12);
    set(gca,'TickLabelInterpreter','latex'); 
    title('Backward Euler', ...
        'Interpreter', 'Latex', 'FontSize', 20);


    % Explicit midpoint subplot
    subplot(4, 1, 3);
    plot(t_list, x_exact, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_ex, X_ex, 'r.-', 'MarkerSize', 10);
    hold  off;
    ylabel('$x(t)$ (-)', 'Interpreter', 'Latex', 'FontSize', 20);
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 20);
    set(gca, 'FontSize', 12)
    set(gca,'TickLabelInterpreter', 'latex')    
    title('Explicit Midpoint', ...
        'Interpreter', 'Latex', 'FontSize', 20);


    % Implicit midpoint subplot
    subplot(4, 1, 4);
    plot(t_list, x_exact, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_im, X_im, 'r.-', 'MarkerSize', 10);
    hold  off;
    ylabel('$x(t)$ (-)', 'Interpreter', 'Latex', 'FontSize', 20);
    set(gca, 'FontSize', 12);
    set(gca,'TickLabelInterpreter','latex');
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 20);
    title('Implicit Midpoint', ...
        'Interpreter', 'Latex', 'FontSize', 20);

    % Make legend
    lgnd = legend('Analytical', 'Numerical', 'Fontsize', 15, ...
        'Interpreter', 'Latex', 'Location', 'East');  
    lgnd.Position(1) = 0.15;
    lgnd.Position(2) = 0.85;
    
    sgtitle(['Integration Stability, $h_{ref}$ = ', h_ref_title], ...
    'interpreter', 'latex', 'fontsize', 25)
    
end