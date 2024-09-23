%% Jacobi 迭代法求线性方程组
function x = jacobi(A, b, tol, max_iter)
    if nargin < 3
        tol = 1e-6; % 默认容忍度
    end
    if nargin < 4
        max_iter = 100; % 默认最大迭代次数
    end
    n = length(b);
    x = zeros(n, 1); % 初始解
    M = diag(diag(A));
    N = M - A;
    S = M \ N; % inv(M) * N

    for k = 1:max_iter
        x_new = S * x + M \ b;

        % 终止条件
        if norm(x_new - x, inf) < tol
            break;
        end
        x = x_new;
    end
end
