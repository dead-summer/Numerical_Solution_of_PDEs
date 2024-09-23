%% 求解五点差分格式的快速 DST 算法
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
lambda = 4 / (h ^2) * sin(pi * (1:I) / (2*I +2)).^2;
mu = 4 / (k ^2) * sin(pi * (1:J) / (2*J +2)).^2;
F = - f(x(2:I+1)', y(2:J+1));
V = dst(dst(F)')';
W = 4 * V ./ ((I + 1) * (J + 1) * (lambda' + mu));

% 快速 DST 计算方程的解
U_dst = dst(dst(W)')';

% 计算相对误差
u = sin(pi * x(2:I+1)') .* sin(pi * y(2:J+1)); % 精确解
disp(['DST 快速算法误差：', num2str(rerror_fro(U_dst, u))]);

% 相对误差
function e = rerror_fro(U_hat, U)
    e = norm(U_hat - U, 'fro') / norm(U, 'fro');
end