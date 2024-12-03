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


%% Functions

function [y, x] = dphi(lambda, A, lower_bound, upper_bound, dCdA, l_elements, max_volume_constrained)

lambda(lambda < 0) = 0;

q0 = -A.^2 .* dCdA;
x = sqrt(q0 ./ (lambda * l_elements));

x(x > upper_bound) = upper_bound;
x(x < lower_bound) = lower_bound;

y = l_elements' * x - max_volume_constrained;

end