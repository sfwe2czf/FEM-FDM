clc;
clear all;
close all;

% Define the domain and boundary conditions
a = 1;
Left = -a/2;
Right = a/2;


%Discretization parameters
N = 50;                        % Numer of intervals
dx = (Right - Left) / N;      % Interval
x = Left : dx : Right;        % Discretized spatial domain
dy = (Right - Left) / N;      % Interval
y = Left : dy : Right;        % Discretized spatial domain
h = a/N;
[X, Y] = meshgrid(x, y);

Radius = sqrt(X.^2 + Y.^2);
Rho2D = -10 * (Radius <= a/6);
figure();
mesh(X, Y, Rho2D);
xlabel('x');
ylabel('y');
title('Source Term (\rho)');

[nodeX, nodeY] = meshgrid(x, y);
nodes = [nodeX(:), nodeY(:)];
elements = delaunay(nodes(:,1), nodes(:,2));

figure();
triplot(elements, nodes(:,1), nodes(:,2), 'k');
xlabel('x');
ylabel('y');
title('Discretized Domain with Triangular Elemnts');
axis equal;

boundaryNodes = find(nodes(:,1) == min(x) | nodes(:,1) == max(x) | nodes(:,2) == min(y) | nodes(:,2) == max(y));
interiorNodes = setdiff(1:size(nodes, 1), boundaryNodes);

figure();
hold on;
plot(nodes(boundaryNodes, 1), nodes(boundaryNodes, 2), 'ro', 'MarkerSize', 8);
plot(nodes(interiorNodes, 1), nodes(interiorNodes, 2), 'ro', 'MarkerSize', 8);
xlabel('x');
ylabel('y');
title('Boundary and Interior Nodes');
legend('Boundary Nodes', 'Interior Nodes');
axis equal;
hold off;

numNodes = size(nodes, 1);
K = sparse(numNodes, numNodes);
F = zeros(numNodes, 1);

function [K_e, F_e] = elementStiffnessLoad(elementCoords, Rho2D, X, Y, h)
    x = elementCoords(:, 1);
    y = elementCoords(:, 2);

    Area = polyarea(x, y);

    b = [y(2) - y(3); y(3) - y(1); y(1) - y(2)];
    c = [x(3) - x(2); x(1) - x(3); x(2) - x(1)];
    B = [b'; c'] / (2 * Area);
    K_e = (B' * B) * Area;

    elementCenter = mean(elementCoords, 1);
    rhoValue = interp2(X, Y, Rho2D, elementCenter(1), elementCenter(2), 'linear', 0);
    F_e = rhoValue * Area / 3 * ones(3, 1);
end

for e1 = 1:size(elements, 1)
    nodeIndicies = elements(e1, :);
    elementCoords = nodes(nodeIndicies, :);

    [K_e, F_e] = elementStiffnessLoad(elementCoords, Rho2D, X, Y, h);

    K(nodeIndicies, nodeIndicies) = K(nodeIndicies, nodeIndicies) + K_e;
    F(nodeIndicies) = F(nodeIndicies) + F_e;
end

K(boundaryNodes, :) = 0;
K(sub2ind(size(K), boundaryNodes, boundaryNodes)) = 1;
F(boundaryNodes) = 0;

V = zeros(numNodes, 1);
tolerance = 1e-12;
maxIterations = 5000;
iterantion = 0;
MSD = 1;
errors = [];

D = diag(diag(K));
L = tril(K, -1);
U = triu(K, 1);

while MSD > tolerance && iterantion < maxIterations
    V_new = D \ (F - (L + U) * V);
    MSD = mean((V_new - V).^2);
    errors = [errors; MSD];
    V = V_new;
    iterantion = iterantion + 1;
end

fprintf('Jacobi method converged in %d iterations.\n', iterantion);

figure();
semilogy(1:iterantion, errors, 'b-', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Mean Squared Distance (Error)');
title('Convergence of Jacobi Iterative Method');
grid on;

V2D = reshape(V, [N+1, N+1]);
figure();
mesh(X, Y, V2D);
xlabel('x');
ylabel('y');
title('Solution (V)');