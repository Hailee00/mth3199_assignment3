
function plotting_midpoint_step
    tspan = [0 10];
    X0 = 1;
    h_values = [0.38, 0.425, .45];
    
    t_exact = linspace(tspan(1),tspan(2), 1000);
    X_exact = solution01(t_exact);
    
    figure;
    hold on;
    
    plot(t_exact,X_exact,'k-','LineWidth',2);
    
    for i = 1:length(h_values)
    
        [t_list,X_list,h_avg,num_evals] = ...
            explicit_midpoint_fixed_step_integration(...
            @rate_func01,...
            tspan,X0,h_values(i));
    
        plot(t_list,X_list,'.-', 'MarkerSize', 20, 'LineWidth', 1.5);
    
    end
    
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);
    
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 17);
    ylabel('$x(t) (-)$', 'Interpreter', 'Latex', 'FontSize', 17);
    title('Explicit Midpoint Example Integration Plot', ...
        'Interpreter', 'Latex', 'FontSize', 20);
    legend('Exact',...
        '$h = 0.38$', '$h = 0.425$', '$h = 45$', 'Interpreter', 'Latex', 'Location', 'southwest', 'Fontsize', 15);

    myfig = gcf(); %set myfig to the current figure
  % exportgraphics(myfig , 'explicit_mdpt_plot.png', 'Resolution', 600); %saves plot to plot_name at high resolution
end