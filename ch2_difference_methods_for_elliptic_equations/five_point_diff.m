%% 考虑 Dirichlet 边界条件的 Poisson 方程
clc, clear

% 定义 f
f = @(x, y) - 2 * pi^2 * sin(pi * x) .* sin(pi * y);

% 边界范围
a = 1; % 0<x<a
b = 1; % 0<y<b

% 划分（本案例中，下面公式按照 h=k 计算）
I = 19; % x 方向 I+1 等分
J = 19; % y 方向 J+1 等分
h = a / (I+1);
k = b / (J+1);

% 初始化
x = linspace(0, a, I + 2);
y = linspace(0, b, J + 2);
F = f(x', y);
U = zeros(I + 2, J + 2);
U(1, :) = 0;
U(end, :) = 0;
U(:, 1) = 0;
U(:, end) = 0;

% 比较两种方法的误差
u = sin(pi * x') .* sin(pi * y); % 精确解
U_jab = jacobi_5p(U, F, h);
U_gs = gauss_seidel_5p(U, F, h);

disp(['Jacobi 迭代法误差：', num2str(rerror_fro(U_jab, u))]);
disp(['Gauss-Seidel 迭代法误差：', num2str(rerror_fro(U_gs, u))]);


% 相对误差
function e = rerror_fro(U_hat, U)
    e = norm(U_hat - U, 'fro') / norm(U, 'fro');
end

% Jacobi 迭代法
function U_new = jacobi_5p(U, F, h, tol, max_iter)
    if nargin < 4
        tol = 1e-6;
    end
    if nargin < 5
        max_iter = 100;
    end
    [I, J] = size(U);
    I = I - 2;
    J = J - 2;
    U_new = zeros(size(U));
    U_new(1, :) = U(1, :);
    U_new(end, :) = U(end, :);
    U_new(:, 1) = U(:, 1);
    U_new(:, end) = U(:, end);

    for k = 1:max_iter
        for i = 2:(I+1)
            for j = 2:(J+1)
                U_new(i, j) = 0.25 * (U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1)) - 0.25 * h^2 * F(i, j);
            end
        end

        % 终止条件
        if norm(U_new - U, inf) < tol
            break;
        end
        U = U_new;
    end
end

% Gauss-Seidel 迭代法
function U_new = gauss_seidel_5p(U, F, h, tol, max_iter)
    if nargin < 4
        tol = 1e-6;
    end
    if nargin < 5
        max_iter = 100;
    end
    [I, J] = size(U);
    I = I - 2;
    J = J - 2;
    U_new = zeros(size(U));
    U_new(1, :) = U(1, :);
    U_new(end, :) = U(end, :);
    U_new(:, 1) = U(:, 1);
    U_new(:, end) = U(:, end);

    for k = 1:max_iter
        for i = 2:(I+1)
            for j = 2:(J+1)
                U_new(i, j) = 0.25 * (U(i+1, j) + U_new(i-1, j) + U(i, j+1) + U_new(i, j-1)) - 0.25 * h^2 * F(i, j);
            end
        end

        % 终止条件
        if norm(U_new - U, inf) < tol
            break;
        end
        U = U_new;
    end
end
