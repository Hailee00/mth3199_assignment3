
function plotting_midpoint_step
    tspan = [0 10];
    X0 = 1;
    h_values = [0.4, 0.2];
    
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
    
        plot(t_list,X_list,'o-');
    
    end
    
    xlabel('Time (-)', 'Interpreter', 'Latex', 'FontSize', 15);
    ylabel('x(t) (-)', 'Interpreter', 'Latex', 'FontSize', 15);
    title('Explicit Midpoint Example Integration Plot', ...
        'Interpreter', 'Latex', 'FontSize', 20);
    legend('Exact',...
        'h = 0.4',...
        'h = 0.2', 'Interpreter', 'Latex', 'Location', 'Best', 'Fontsize', 15);
end