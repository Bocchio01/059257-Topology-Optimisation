clc
clear variables
close all

%% Structure and Applied loads

structure = load_structure('models/simple_structure.inp');
structure.nodes(3).forceX = 5000;
structure.nodes(4).forceX = 5000;
structure.nodes(4).forceY = 5000;


%% Optimization problem constraints

alpha = 0.5; % []
lower_bound = 1e-4; % [mm^2]
upper_bound = 800;  % [mm^2]

max_volume = upper_bound * sum([structure.elements.L], 'all');
max_volume_constrained = max_volume * alpha;


%% Problem solution

epsilon = 1e-5;
N = 10;
lambda = 1;
results = struct( ...
    'C', zeros(N, 1), ...
    'V', zeros(N, 1), ...
    'alpha', zeros(N, 1), ...
    'dCdA', zeros(N, length(structure.elements)), ...
    'A', zeros(N, length(structure.elements)), ...
    'stress', zeros(N, length(structure.elements)));

for k = 1 : N

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
    
    % Solve sum(L * A(lambda*)) - V_{0} == 0 and compute A(lambda*)
    lambda = fsolve(handler, lambda);
    [~, A] = handler(lambda);

    % Variables updates
    for ii = 1:length(structure.elements)
        structure.elements(ii).A = A(ii);
    end
    results.V(k) = A' * [structure1.elements.L]';
    results.alpha(k) = results.V(k) / max_volume;
    results.A(k, :) = A;
    results.stress(k, :) = structure.computeStresses();

    % Check convergence condition
    if(norm(results.dCdA(k, :)) < epsilon)
        break;
    end

end

clear ii k
clear epsilon N A handler lambda
clear ue k0 dofIndices node2Index node1Index


%% Plots

run("figures\routine.m")


%% Functions

function [y, x] = dphi(lambda, A, lower_bound, upper_bound, dCdA, l_elements, max_volume_constrained)

lambda(lambda < 0) = 0;

q0 = -A.^2 .* dCdA;
x = sqrt(q0 ./ (lambda * l_elements));

x(x > upper_bound) = upper_bound;
x(x < lower_bound) = lower_bound;

y = l_elements' * x - max_volume_constrained;

end
