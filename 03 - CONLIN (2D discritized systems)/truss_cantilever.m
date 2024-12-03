clc
clear variables
close all

%% Structure and Applied loads

[undeformed_structure, structure] = deal(load_structure('models/truss_cantilever.inp'));

structure.nodes(49).forceY = -6000;


%% Optimization problem constraints

% alpha = 0.20; % []
lower_bound = 1e-4; % [mm^2]
upper_bound = 8e+2;  % [mm^2]

max_volume = upper_bound * sum([structure.elements.L], 'all');
max_volume_constrained = max_volume * alpha;


%% Problem solution

epsilon = 1e-5;
N = 100;
lambda = 1;
results = struct( ...
    'C', zeros(N, 1), ...
    'V', zeros(N, 1), ...
    'alpha', zeros(N, 1), ...
    'dCdA', zeros(N, length(structure.elements)), ...
    'A', zeros(N, length(structure.elements)), ...
    'stress', zeros(N, length(structure.elements)));

for k = 1 : N

    if(mod(k, 25) == 0)
        disp(['Iteration #', num2str(k)]);
    end

    % FEM Analysis
    structure = structure.runFEM();
    results.C(k) = structure.F' * structure.U;

    % Sensitivity Analysis
    for ii = 1:length(structure.elements)

        k0 = structure.elements(ii).getStiffnessMatrix / structure.elements(ii).A;
        ue = structure.U(structure.getDOFidxs(structure.elements(ii)));

        results.dCdA(k, ii) = - ue' * k0 * ue;

    end

    % Dual problem solution
    handler = @(lambda) dphi( ...
        lambda, ...
        [structure.elements.A]', ...
        lower_bound, ...
        upper_bound, ...
        results.dCdA(k, :)', ...
        [structure.elements.L]', ...
        max_volume_constrained);

    % Solve sum(L * A(lambda*)) - V_{0} == 0 and compute lambda*
    lambda = fsolve(handler, lambda, optimoptions('fsolve', 'Display', 'off'));
    [~, A] = handler(lambda);

    % Variables updates
    for ii = 1:length(structure.elements)
        structure.elements(ii).A = A(ii);
    end
    results.V(k) = A' * [structure.elements.L]';
    results.alpha(k) = results.V(k) / max_volume;
    results.A(k, :) = A;
    results.stress(k, :) = structure.computeStresses();

    % Check convergence condition
    if(norm(results.dCdA(k, :)) < epsilon)
        break;
    end

end

clear ii
clear epsilon N A handler lambda
clear ue k0 dofIndices node2Index node1Index


%% Results

fprintf([ ...
    'Cross-sectional areas elements [64,77,82]: [%.2f %.2f %.2f] [mm^2]\n' ...
    'Stresses elements [64,77,82]: [%.2f %.2f %.2f] [N/mm^2]\n' ...
    ], results.A(k, [64,77,82]), results.stress(k, [64,77,82]));



%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 15);
set(0, 'DefaultLineLineWidth', 2);

try

    plot_struct.export_path = '../latex/img/MATLAB/truss_cantilever';
    plot_struct.data = cell(0);

    % Undeformed
    % figure_undeformed_structure = figure('Name', 'Undeformed structure');
    %
    % nexttile
    % hold on
    % grid on
    %
    % plot_structure(undeformed_structure.runFEM(), struct('color', '#D95319', 'LineWidth', 2, 'force', true));
    %
    % axis equal
    % axis padded
    % title('Undeformed structure')
    % xlabel('[mm]')
    % ylabel('[mm]')
    % ylim([-500 4500])
    % xlim([-500 8500])


    % Results
    figure_topology_results = figure('Name', 'Topology optimization results');
    tile = tiledlayout(2, 3);

    % Compliance
    nexttile(tile, 1)
    hold on
    grid on

    plot(results.C(1:k), 'ob');

    title(['Optimization compliance @alpha=' num2str(alpha, '%.2f')])
    xlabel('Iteration #')
    ylabel('Compliance [Nmm]')

    % Volume
    nexttile(tile, 4)
    hold on
    grid on

    plot(results.alpha(1:k) * 100, 'ob');

    ylim([0 30])

    title(['Volume fraction @alpha=' num2str(alpha, '%.2f')])
    xlabel('Iteration #')
    ylabel('Volume fraction []')

    ytickformat('percentage')


    % Deformed
    nexttile(tile, 2, [2, 2])
    hold on
    grid on

    plot_structure(structure);

    title(['Deformed structure @alpha=' num2str(alpha, '%.2f')])
    xlabel('[mm]')
    ylabel('[mm]')


    % plot_struct.data{end+1} = {figure_undeformed_structure, '/undeformed_structure'};
    plot_struct.data{end+1} = {figure_topology_results, ['/results_alpha' num2str(alpha*100)]};

    export_pdf_figure(plot_struct);
    clear plot_struct plot_idx current_plot local_path filename tile

catch
    disp('Could not export figure')
end



%% Functions

function [y, x] = dphi(lambda, A, lower_bound, upper_bound, dCdA, l_elements, max_volume_constrained)

lambda(lambda < 0) = 0;

q0 = -A.^2 .* dCdA;
x = sqrt(q0 ./ (lambda * l_elements));

x(x > upper_bound) = upper_bound;
x(x < lower_bound) = lower_bound;

y = l_elements' * x - max_volume_constrained;

end
