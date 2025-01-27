function x = TO2_distributed(nx, ny, alpha, p, rmin)

x(1:ny,1:nx) = alpha;
change = 1.;

while change > 0.01
    xold = x;

    % FEM Analysis
    U = compute_U(nx, ny, x, p);

    % Sensitivity Analysis
    [KE] = compute_k_local();
    c = 0.;
    for ely = 1:ny
        for elx = 1:nx
            n1 = (ny+1)*(elx-1)+ely;
            n2 = (ny+1)* elx   +ely;
            Ue = U([ ...
                2*n1-1; ...
                2*n1; ...
                2*n2-1; ...
                2*n2; ...
                2*n2+1; ...
                2*n2+2; ...
                2*n1+1; ...
                2*n1+2 ...
                ], 1);
            c = c + x(ely,elx)^p * Ue' * KE * Ue;
            dc(ely,elx) = -p*x(ely,elx)^(p-1) * Ue' * KE * Ue;
        end
    end

    % Filtering Sensitivities
    dc = apply_filter_to_dc(nx, ny, rmin, x, dc);

    % Design variable update
    x = optimality_criteria(nx, ny, x, alpha, dc);

    % colormap(gray)
    % imagesc(-x)
    % axis equal
    % axis tight
    % pause(1e-6);

    change = max(max(abs(x - xold)));
end

end

%%%%%%%%%% OPTIMALITY CRITERIA UPDATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function xnew = optimality_criteria(nx,ny,x,alpha,dc)

l1 = 0; l2 = 100000; move = 0.2;
while (l2-l1 > 1e-4)
    lmid = 0.5*(l2+l1);
    xnew = max(0.001,max(x-move,min(1.,min(x+move,x.*sqrt(-dc./lmid)))));
    if sum(sum(xnew)) - alpha*nx*ny > 0;
        l1 = lmid;
    else
        l2 = lmid;
    end
end
end


function dcn = apply_filter_to_dc(nx, ny, rmin, x, dc)

dcn = zeros(ny, nx);

for x_idx = 1:nx
    for y_idx = 1:ny
        
        sum = 0.0;
        
        for k = max(x_idx-floor(rmin),1):min(x_idx+floor(rmin),nx)
            for l = max(y_idx-floor(rmin),1):min(y_idx+floor(rmin),ny)
        
                fac = rmin-sqrt((x_idx-k)^2+(y_idx-l)^2);
                sum = sum + max(0, fac);
                dcn(y_idx,x_idx) = dcn(y_idx,x_idx) + max(0,fac)*x(l,k)*dc(l,k);

            end
        end
        
        dcn(y_idx,x_idx) = dcn(y_idx,x_idx)/(x(y_idx,x_idx)*sum);

    end
end

end



function U = compute_U(nx, ny, x, p)

k_local = compute_k_local();

K = sparse(2*(nx+1)*(ny+1), 2*(nx+1)*(ny+1));
F = sparse(2*(ny+1)*(nx+1), 1);
U = zeros(2*(ny+1)*(nx+1), 1);

for x_idx = 1:nx
    for y_idx = 1:ny
        n1 = (ny+1)*(x_idx-1) + y_idx;
        n2 = (ny+1)* x_idx    + y_idx;
        edof = [2*n1-1; 2*n1; 2*n2-1; 2*n2; 2*n2+1; 2*n2+2; 2*n1+1; 2*n1+2];
        K(edof, edof) = K(edof, edof) + x(y_idx,x_idx)^p*k_local;
    end
end

% Define loads & boundary conditions
F(2*(nx+1)*(ny+1),1) = -1;
fixeddofs   = 1:2*(ny+1);
alldofs     = 1:2*(ny+1)*(nx+1);
freedofs    = setdiff(alldofs,fixeddofs);

% Solve for displacements
U(freedofs, :) = K(freedofs,freedofs) \ F(freedofs,:);
U(fixeddofs, :) = 0;

end


function k_local = compute_k_local()

E = 1;
nu = 0.3;

k = [1/2-nu/6 1/8+nu/8 -1/4-nu/12 -1/8+3*nu/8 -1/4+nu/12 -1/8-nu/8  nu/6 1/8-3*nu/8];

k_local = E/(1-nu^2) * [ ...
    k(1) k(2) k(3) k(4) k(5) k(6) k(7) k(8)
    k(2) k(1) k(8) k(7) k(6) k(5) k(4) k(3)
    k(3) k(8) k(1) k(6) k(7) k(4) k(5) k(2)
    k(4) k(7) k(6) k(1) k(8) k(3) k(2) k(5)
    k(5) k(6) k(7) k(8) k(1) k(2) k(3) k(4)
    k(6) k(5) k(4) k(3) k(2) k(1) k(8) k(7)
    k(7) k(4) k(5) k(2) k(3) k(8) k(1) k(6)
    k(8) k(3) k(2) k(5) k(4) k(7) k(6) k(1)];

end