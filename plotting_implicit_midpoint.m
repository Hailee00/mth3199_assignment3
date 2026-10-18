
function plotting_implicit_midpoint
    tspan = [0 10];
    X0 = 1;
    h_values = [0.3, 1.25, 1.75];


    t_exact = linspace(tspan(1),tspan(2),1000);
    X_exact = solution01(t_exact);

    figure;
    hold on;

    plot(t_exact,X_exact,'k-','LineWidth',2);

    for i = 1:length(h_values)

        [t_list,X_list,h_avg,num_evals] = ...
            implicit_fixed_step_integration(...
            @rate_func01,...
            @implicit_midpoint_step,...
            tspan,X0,h_values(i));

        plot(t_list,X_list,'.-', 'MarkerSize', 20, 'LineWidth', 1.5);

    end

    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 17);
    ylabel('$x(t) (-)$', 'Interpreter', 'Latex', 'FontSize', 17);
    title('Implicit Midpoint Example Integration Plot', ...
        'Interpreter', 'Latex', 'FontSize', 20);

    legend('Exact','$h = 0.38$', '$h = 1.25$', '$h = 1.75$', ...
        'Interpreter', 'Latex', 'Location', 'northeast', 'FontSize', 14);

    hold off;

     myfig = gcf(); %set myfig to the current figure
  % exportgraphics(myfig , 'implicit_mdpt_plot.png', 'Resolution', 600); %saves plot to plot_name at high resolution
end
