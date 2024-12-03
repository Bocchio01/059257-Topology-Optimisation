clc
clear variables
close all

%% Problem definition

E = 210000;
P = 5000;
L = 1;
C0 = 50;

A1 = sym('A1', 'positive');
A2 = sym('A2', 'positive');
A3 = sym('A3', 'positive');
A4 = sym('A4', 'positive');
A5 = sym('A5', 'positive');
lambda = sym('lambda', 'positive');

% Define g0 and g1 as symbolic expressions
g0 = A1 + A2 + A3 + A4 + sqrt(2)*A5;
g1 = 1/A1 + 1/A2 + 4/A3 + 4/A4 + 8*sqrt(2)/A5 - (E*C0)/(P^2*L);

% Define the Lagrangian as a symbolic expression
Lagrangian = g0 + lambda * g1;

% Calculate the gradient of the Lagrangian with respect to A1, A2, A3, A4, A5
gradientL = gradient(Lagrangian, [A1, A2, A3, A4, A5]);

% Solve for lambda using the gradient conditions
A_of_lambda = solve(gradientL == 0, [A1, A2, A3, A4, A5], 'ReturnConditions', true);

% Extract the solutions for A1, A2, A3, A4, A5
A1_lambda = A_of_lambda.A1;
A2_lambda = A_of_lambda.A2;
A3_lambda = A_of_lambda.A3;
A4_lambda = A_of_lambda.A4;
A5_lambda = A_of_lambda.A5;

% Substitute the solutions into the Lagrangian to get it in terms of lambda
phi = subs(Lagrangian, ...
    [A1, A2, A3, A4, A5], ...
    [A1_lambda, A2_lambda, A3_lambda, A4_lambda, A5_lambda]);

% Solve for lambda_star
lambda_star = solve(phi == 0, lambda, 'ReturnConditions', true);

% Extract the value of lambda_star
lambda_value = lambda_star.lambda;

% Substitute lambda_star into A1_lambda, A2_lambda, A3_lambda, A4_lambda, A5_lambda
A1_star = subs(A1_lambda, lambda, lambda_value);
A2_star = subs(A2_lambda, lambda, lambda_value);
A3_star = subs(A3_lambda, lambda, lambda_value);
A4_star = subs(A4_lambda, lambda, lambda_value);
A5_star = subs(A5_lambda, lambda, lambda_value);

% Display the final values of A1, A2, A3, A4, and A5
disp('A1_star = ');
disp(A1_star);
disp('A2_star = ');
disp(A2_star);
disp('A3_star = ');
disp(A3_star);
disp('A4_star = ');
disp(A4_star);
disp('A5_star = ');
disp(A5_star);