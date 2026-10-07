function plot_local_trunkation
    h_list = logspace(10e-5, 10, 100);
    tspan = [0 10];
    X0 = 1;
    t_ex = 0.492;
    error_list = zeros(size(h_list));

    for i = 1:length(h_list)
        [G,~] = forward_euler_step(@rate_func01,t_ex,X0,h_list(i));
        X = rate_func01(t_ex+h_list(i), X0);
        error_list(i) = norm(G-X);

        global_trunkation_euler(i) = sum(error_list(i))
    end

    mdl = polyfit(log(h_list),log(error_list),1)

    figure();
    hold on;
    loglog(h_list, error_list, 'ro')
    title('Local Trunkation Error for Forward Euler', 'Interpreter', 'Latex', 'FontSize', 20)
    xlabel('Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 15)
    ylabel('Error Size (-)', 'Interpreter', 'Latex', 'FontSize', 15)
    loglog(h_list, h_list.*mdl(1)+h_list.*mdl(2))
    % set(gca,'xtick',[10^0:6:10^10], 'ytick', [10^0:7:10^12])
    hold off

    figure();
    hold on;
    loglog(h_list, global_trunkation_euler, 'ro')
    title('Global Trunkation Error', 'Interpreter', 'Latex', 'FontSize', 20)
    xlabel('Step Size: h (-)', 'Interpreter', 'Latex', 'FontSize', 15)
    ylabel('Error Size (-)', 'Interpreter', 'Latex', 'FontSize', 15)
    loglog(h_list, h_list.*mdl(1)+h_list.*mdl(2))
    % set(gca,'xtick',[10^0:6:10^10], 'ytick', [10^0:7:10^12])
    hold off



end