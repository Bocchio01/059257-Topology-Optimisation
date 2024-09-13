clc
clear variables
close all

%% Optimization problem definition

f  = @(x1, x2) 3*x1 + x2;
g1 = @(x1, x2) (2*x1+x2)./(2*x1.*(x1+x2)) - 1;
g2 = @(x1, x2) 1./(sqrt(2)*(x1+x2)) - 1;
g3 = @(x1, x2) x2./(2*x1.*(x1+x2)) - 1;


%% Problem solution
% From the plot, we observe that solution must be found on the g1(x1, x2)
% constraint that minimize f(x1, x2)

x = linspace(0, 1, 100);
y = linspace(0, 1, 100);
[X, Y] = meshgrid(x, y);

M = contourc(g1(X, Y), [0,0]);

M(1, 2:end) = M(1, 2:end) / length(x);
M(2, 2:end) = M(2, 2:end) / length(y);

[val, idx] = min(f(M(1, 2:end), M(2, 2:end)));

solution_point.x1 = M(1, idx);
solution_point.x2 = M(2, idx);
solution_point.f = f(solution_point.x1, solution_point.x2);



%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
% set(0, 'defaultaxesfontsize', 15);
% set(0, 'DefaultLineLineWidth', 2);


plot_struct.flags = true;
plot_struct.export_flag = false;
% plot_struct.export_path = 'latex/img/MATLAB';
% plot_struct.data = cell(0);

if (plot_struct.flags)
    run("figures\fig_01_problem_visualization.m");
end

if (plot_struct.export_flag)
    pause(1);
    for plot_idx = 1:numel(plot_struct.data)

        current_plot = plot_struct.data{plot_idx};
        tile = current_plot{1};
        local_path = current_plot{2};

        filename = [plot_struct.export_path local_path '.pdf'];
        exportgraphics(tile, filename, 'ContentType', 'vector');

    end
end


