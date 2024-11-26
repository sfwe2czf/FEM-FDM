clc;
clear all;
close all;

% Define the domain and boundary conditions
a = 1;
Left = -1/2;
Right = 1/2;
V=0;

%Discretization parameters
N = 50;                        % Numer of intervals
dx = (Right - Left) / N;      % Interval
x = Left : dx : Right;        % Discretized spatial domain
dy = (Right - Left) / N;      % Interval
y = Left : dy : Right;        % Discretized spatial domain

%Create the diagonal and off-diagonal matrices
A_Diag = -4 * eye(N-1) + diag(ones(N-2, 1), 1) + diag(ones(N-2, 1), -1);
A_offDiag = eye(N-1);

%Initialize the final matrix A
A = zeros((N-1)*(N-1), (N-1)*(N-1));

%Fill the diagonal blocks with A_diag
for i = 1:N-1
    row_start = 1 + (i-1) * (N-1);
    row_end = row_start + (N-1) - 1;
    col_start = row_start;
    col_end = row_end;
    A(row_start:row_end, col_start:col_end) = A_Diag;
end

%Fill the off-diagonal blocks
for i = 1:N-2
    row_start = 1 + (i-1) * (N-1);
    row_end = row_start + (N-1) - 1;
    col_start = row_start + (N-1);
    col_end = col_start + (N-1) - 1;

    % Upper off-diagonal block
    A(row_start:row_end, col_start:col_end) = A_offDiag;

    % Lower off-diagonal block
    A(col_start:col_end, row_start:row_end) = A_offDiag;
end


% Construct Jacobi Matrix for finite difference method
D = diag(diag(A));              % Diagonal matrix (second derivative)
R = A - D;                 % Initialize off-diagonal matrix
H = -inv(D) * R;                % Iteration matrix

% Initialize the right-hand size vector
b = zeros((N-1)^2, 1);
for i = 1 : (N - 1)
    for j = 1 : (N - 1)
        if (x(j+1)^2+y(i+1)^2) <= a^2/36
            b((N-1)*(i-1)+j) = -10;
        end 
    end
end
    
b1 = b * -dx^2;
d = inv(D) * b1;

% Iteration and sequential plot
figure();
V_initial = zeros((N-1)^2, 1);      % Initial voltage distribution
V_evolution = V_initial;        % Current voltage estimate

[xx, yy] = meshgrid(x, y);

% Iterate to solve for voltage distribution
for iteration = 1 : 10000

    % Update voltafe estimate
    V_evolution = H * V_evolution + d;

    % Plot the voltage distribution at specified iterations
    if ismember(iteration, [10000])
        V_mat = zeros(N+1, N+1);
        for i = 1 : N+1
            for j = 1 : N+1
                if  i == 1 | j == 1 | i == N+1 | j == N+1
                    V_mat(i, j)=0;
                else
                    V_mat(i, j)=V_evolution((j-2)*(N-1)+i-1);
                end 
            end
        end
        surf(xx, yy, V_mat)
        hold on; % Retain current plot for overlaying
    end

end
% Labeling the plot
xlabel('x')
ylabel('y')
zlabel('V')

figure();
b_mat = zeros(N+1, N+1);
        for i = 1 : N+1
            for j = 1 : N+1
                if(x(i)^2+y(j)^2<a^2/36)
                    b_mat(i, j) = -10;
                else
                    b_mat(i, j) = 0;
                end
            end
        end
surf(xx, yy, b_mat);
xlabel('x')
ylabel('y')
zlabel('\rho')
hold off;