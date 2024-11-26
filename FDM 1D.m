clc;
clear all;
close all;

% Define the domain and boundary conditions
a = 1;
xLeft = -a;
xRight = a;
VLeft = 0;
VRight = a^2/4;

%Discretization parameters
N = 100;                        % Numer of intervals
dx = (xRight - xLeft) / N;      % Interval
x = xLeft : dx : xRight;        % Discretized spatial domain

% Construct Jacobi Matrix for finite difference method
D = -2 * eye(N-1);              % Diagonal matrix (second derivative)
R = zeros(N-1);                 % Initialize off-diagonal matrix

% Fill the off-diagonal elements
for i = 1 : N-2
    R(i, i+1) = 1;
    R(i+1, i) = 1;
end

A = D + R;                      % Jacobi matrix
H = -inv(D) * R;                % Iteration matrix

% Initialize the right-hand size vector
b = zeros(N-1, 1);
for i = 2 : N - 2
    if x(i) < -a/2
        b(i) = 2;
    elseif x(i) < 0
        b(i) = 1;
    elseif x(i) < a / 2
        b(i) = -1;
    else
        b(i) = -2;
    end    
end
    
b = b * dx^2
b(1) = 2*dx^2 - VLeft;                 % Boundary condition at the left
b(end) = -2*dx^2 - VRight;              % Boundary condition at the right

d = inv(D) * b;

% Iteration and sequential plot
figure();
V_initial = linspace(0, a*a/4, N-1);      % Initial voltage distribution
V_evolution = transpose(V_initial);        % Current voltage estimate

% Iterate to solve for voltage distribution
for iteration = 1 : 10000

    % Update voltafe estimate
    V_evolution = H * V_evolution + d;

    % Plot the voltage distribution at specified iterations
    if ismember(iteration, [1, 1000, 5000, 10000])
        plot(x, [VLeft; V_evolution; VRight], 'Color', [1/10000 * iteration, 0, 0, 1])
        hold on; % Retain current plot for overlaying
    end
end

% Labeling the plot
xlabel('x')
ylabel('V(x)')