clc
clear variables
close all

for w1 = [0.0 0.2 0.5 0.8 1.0]
%% Structure and Applied loads

[undeformed_structure, structure_brk, structure_acc] = deal(load_structure('models/motorcycle_frame.inp'));

% Braking
structure_brk.nodes(8).forceX = -13000;
structure_brk.nodes(8).forceY = +2000;
structure_brk.nodes(9).forceX = +9000;
structure_brk.nodes(9).forceY = +0;

% Acceleration
structure_acc.nodes(1).forceX = +2000;
structure_acc.nodes(1).forceY = -8000;


%% Optimization problem constraints

alpha = 0.1; % []
lower_bound = 1e-3; % [mm^2]
upper_bound = 700;  % [mm^2]
% w1 = 0.5;
w2 = 1 - w1;

max_volume = upper_bound * sum([structure_brk.elements.L], 'all');
max_volume_constrained = max_volume * alpha;


%% Problem solution

epsilon = 1e-5;
N = 200;
lambda = 1;
results = struct( ...
    'C', zeros(N, 1), ...
    'V', zeros(N, 1), ...
    'alpha', zeros(N, 1), ...
    'dCdA', zeros(N, length(structure_brk.elements)), ...
    'A', zeros(N, length(structure_brk.elements)), ...
    'stress', zeros(N, length(structure_brk.elements)));

for k = 1 : N

    if(mod(k, 25) == 0)
        disp(['Iteration #', num2str(k)]);
    end

    % FEM Analysis
    structure_brk = structure_brk.runFEM();
    structure_acc = structure_acc.runFEM();
    results.C(k) = w1 * structure_brk.F' * structure_brk.U + ...
                   w2 * structure_acc.F' * structure_acc.U;

    % Sensitivity Analysis
    for ii = 1:length(structure_brk.elements)

        k01 = structure_brk.elements(ii).getStiffnessMatrix / structure_brk.elements(ii).A;
        ue1 = structure_brk.U(structure_brk.getDOFidxs(structure_brk.elements(ii)));

        k02 = structure_acc.elements(ii).getStiffnessMatrix / structure_acc.elements(ii).A;
        ue2 = structure_acc.U(structure_acc.getDOFidxs(structure_acc.elements(ii)));

        results.dCdA(k, ii) = w1 * (- ue1' * k01 * ue1) + ...
                              w2 * (- ue2' * k02 * ue2);

    end

    % Dual problem solution
    handler = @(lambda) dphi( ...
        lambda, ...
        [structure_brk.elements.A]', ...
        lower_bound, ...
        upper_bound, ...
        results.dCdA(k, :)', ...
        [structure_brk.elements.L]', ...
        max_volume_constrained);
    
    % Solve sum(L * A(lambda*)) - V_{0} == 0 and compute A(lambda*)
    lambda = fsolve(handler, lambda, optimoptions('fsolve', 'Display', 'off'));
    [~, A] = handler(lambda);

    % Variables updates
    for ii = 1:length(structure_brk.elements)
        structure_brk.elements(ii).A = A(ii);
        structure_acc.elements(ii).A = A(ii);
    end
    results.V(k) = A' * [structure_brk.elements.L]';
    results.alpha(k) = results.V(k) / max_volume;
    results.A(k, :) = A;
    % results.stress(k, :) = structure_brk.computeStresses();

    % Check convergence condition
    if(norm(results.dCdA(k, :)) < epsilon)
        break;
    end

end

clear ii
clear epsilon N A handler lambda
clear ue1 ue2 k01 k02 dofIndices node2Index node1Index


%% Plots

reset(0)
set(0, 'DefaultFigureNumberTitle', 'off');
set(0, 'DefaultFigureWindowStyle', 'docked');
set(0, 'defaultaxesfontsize', 15);
set(0, 'DefaultLineLineWidth', 2);

try

    plot_struct.export_path = 'latex/img/MATLAB/motorcycle_frame';
    plot_struct.data = cell(0);

    % Undeformed
    figure_undeformed_structure = figure('Name', 'Undeformed structure');

    nexttile
    hold on
    grid on

    plot_structure(undeformed_structure.runFEM(), struct('color', '#D95319', 'LineWidth', 2, 'force', true));

    axis equal
    axis padded
    title('Undeformed structure')
    xlabel('[mm]')
    ylabel('[mm]')
    ylim([-100 550])
    xlim([-250 850])
    


    % Compliance
    figure_comliance = figure('Name', 'Compliance');
    tile = tiledlayout(2, 1);

    % Compliance
    nexttile;
    hold on
    grid on

    plot(results.C(1:k), 'ob');

    title(['Optimization compliance @alpha=' num2str(alpha, '%.2f')])
    xlabel('Iteration #')
    ylabel('Compliance [Nmm]')


    % Results
    figure_topology_results = figure('Name', 'Topology optimization results');
    tile = tiledlayout(1, 2);

    % Deformed braking
    nexttile(tile, 1)
    hold on
    grid on

    plot_structure(structure_brk);

    title(['Braking load case @w1=' num2str(w1, '%.1f')])
    xlabel('[mm]')
    ylabel('[mm]')

    % Deformed accelrating
    nexttile(tile, 2)
    hold on
    grid on

    plot_structure(structure_acc);

    title(['Acceleration load case @w2=' num2str(w2, '%.1f')])
    xlabel('[mm]')
    ylabel('[mm]')


    % plot_struct.data{end+1} = {figure_undeformed_structure, '/undeformed_structure'};
    % plot_struct.data{end+1} = {figure_comliance, ['/compliance_alpha' num2str(alpha*100)]};
    plot_struct.data{end+1} = {figure_topology_results, ['/results_w' num2str(w1*100)]};

    export_pdf_figure(plot_struct);
    clear plot_struct plot_idx current_plot local_path filename tile

catch
    disp('Could not export figure')
end

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
