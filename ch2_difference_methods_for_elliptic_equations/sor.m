function x = sor(A, b, omega, tol, max_iter)
    % SOR 求解线性方程组 Ax=b
    if nargin < 3 
        omega = 1; % 默认值为 1
    end
    if nargin < 4
        tol = 1e-6; % 默认容忍度
    end
    if nargin < 5
        max_iter = 100; % 默认最大迭代次数
    end
    n = length(b);
    x = zeros(n, 1);
    M = tril(A);
    N = M - A;
    S = omega * M \ N + (1 - omega) * eye(n); % 迭代矩阵
    b1 = omega * M \ b;

    for k=1:max_iter
        x_new = S * x + b1;

        % 终止条件
        if norm(x_new - x, inf) < tol
            break;
        end
        x = x_new;
    end 
end