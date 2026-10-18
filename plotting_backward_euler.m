
function plotting_backward_euler
    tspan = [0 10];
    X0 = 1;
    h_values = [0.25, 0.85, 2];

    t_exact = linspace(tspan(1),tspan(2),1000);
    X_exact = solution01(t_exact);

    figure;
    hold on;

    plot(t_exact,X_exact,'k-','LineWidth',2);

    for i = 1:length(h_values)

        [t_list,X_list,h_avg,num_evals] = ...
            implicit_fixed_step_integration(...
            @rate_func01,...
            @backward_euler_step,...
            tspan,X0,h_values(i));

         plot(t_list,X_list,'.-', 'MarkerSize', 20, 'LineWidth', 1.5);

    end
    
    set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12);

    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 18);
    ylabel('$x(t) (-)$', 'Interpreter', 'Latex', 'FontSize', 18);
    title('Backward Euler Example Integration Plot', ...
        'Interpreter', 'Latex', 'FontSize', 22);

   legend('Exact',...
        '$h = 0.25$', '$h = 0.85$', '$h = 2.0$', 'Interpreter', 'Latex', 'Location', 'northeast', 'Fontsize', 13);

    hold off;
end
