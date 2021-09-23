clear;
clc;

x0 = [0.1 0.05]; % 起始點
A = []; % 線性不等式拘束條件的係數矩陣
b = []; % 線性不等式拘束條件的係數向量 AX <= b
Aeq = []; % 線性不等式拘束條件的係數向量
beq = []; % 線性等式拘束條件的係數向量 AeqX = beq
ub = [0.5 0.5]; % 設計空間的 upper bounds
lb = [0.001 0.001]; % 設計空間的 lower bounds
options = optimset ('display', 'off', 'Algorithm', 'sqp');%c演算法的參數設定
[x, fval, exitflag] = fmincon(@(x)obj(x), x0, A, b, Aeq, beq, lb, ub, @(x)nonlcon(x), options);
% 執行fmincon，輸出最佳解,x; 最佳目標函數值,fval; 收斂情形,exitflag
% obj為目標函數，nonlcon 為 (非線性) 拘束條件

function f = obj(x)
    ro = 7860;
    L1 = 9.14; L2 = 9.14*sqrt(2);
    A1 = pi()*x(1)^2; A2 = pi()*x(2)^2;
    m1 = A1*L1*ro; m2 = A2*L2*ro;
    f = 6*m1 + 4*m2;
end

function [g, geq] = nonlcon(x)
    [sigma, Q] = sol_TenBarTruss(x(1), x(2));
    sigma_y = 250e06;
    g(1) = sqrt(Q(3)^2+Q(4)^2) - 0.02;
    for i = 1:10
        g(i+1) = abs(sigma(i)) - sigma_y;
    end
    geq = [];
    % 顯示迭代過程及相關結果
    fprintf('r1 = %f m , r2 = %f m\n', x(1), x(2));
    fprintf('displacement at node 2 = %f m\n', sqrt(Q(3)^2+Q(4)^2));
    fprintf('sigma_y(GPa) =')
    disp(sigma/1e06);
    ro = 7860;
    L1 = 9.14; L2 = 9.14*sqrt(2);
    A1 = pi()*x(1)^2; A2 = pi()*x(2)^2;
    m1 = A1*L1*ro; m2 = A2*L2*ro;
    f = 6*m1 + 4*m2;
    fprintf('mass = %f kg\n', f)
end