function [sigma, Q] = sol_TenBarTruss(r1, r2)
    E = 200e09; ro = 7860; sigma_y = 250e06;
    L1 = 9.14; L2 = 9.14*sqrt(2);
    F = 1e07;
    A1 = pi()*r1^2; A2 = pi()*r2^2;
    m1 = A1*L1*ro; m2 = A2*L2*ro;
    nodes = 6;

    K = zeros(nodes*2, nodes*2);   % stiffness matrix 

    K = add_element(K, A1, E, L1, 0, 5, 3); % link 1
    K = add_element(K, A1, E, L1, 0, 3, 1); % link 2
    K = add_element(K, A1, E, L1, 0, 6, 4); % link 3
    K = add_element(K, A1, E, L1, 0, 4, 2); % link 4

    K = add_element(K, A1, E, L1, 90, 4, 3); % link 5
    K = add_element(K, A1, E, L1, 90, 2, 1); % link 6

    K = add_element(K, A2, E, L2, -45, 5, 4); % link 7
    K = add_element(K, A2, E, L2, -45, 3, 2); % link 9

    K = add_element(K, A2, E, L2, 45, 6, 3); % link 8
    K = add_element(K, A2, E, L2, 45, 4, 1); % link 10

    F_matrix = [0; 0; 0; -F; 0; 0; 0; -F; 0; 0; 0; 0];

    Q = zeros(12, 1);
    Q(1:8, :) = K(1:8, 1:8)\F_matrix(1:8, :);

    sigma(1) = compute_stress(Q, E, L1, 0, 5, 3);
    sigma(2) = compute_stress(Q, E, L1, 0, 3, 1);
    sigma(3) = compute_stress(Q, E, L1, 0, 6, 4);
    sigma(4) = compute_stress(Q, E, L1, 0, 4, 2);

    sigma(5) = compute_stress(Q, E, L1, 90, 4, 3);
    sigma(6) = compute_stress(Q, E, L1, 90, 2, 1);

    sigma(7) = compute_stress(Q, E, L2, -45, 5, 4);
    sigma(9) = compute_stress(Q, E, L2, -45, 3, 2);

    sigma(8) = compute_stress(Q, E, L2, 45, 6, 3);
    sigma(10) = compute_stress(Q, E, L2, 45, 4, 1);

    R = K(9:12, :)*Q - F_matrix(9:12, 1);   % reaction force at node5 and node6

    function K = add_element(K, A, E, L, theta, node1, node2)
        c = cosd(theta); s = sind(theta);
        temp = A*E/L*[c^2 c*s; c*s s^2];
        K((2*node1-1):(2*node1), (2*node1-1):(2*node1))...
        = K((2*node1-1):(2*node1), (2*node1-1):(2*node1)) + temp;
        K((2*node2-1):(2*node2), (2*node2-1):(2*node2))...
        = K((2*node2-1):(2*node2), (2*node2-1):(2*node2)) + temp;
        K((2*node1-1):(2*node1), (2*node2-1):(2*node2))...
        = K((2*node1-1):(2*node1), (2*node2-1):(2*node2)) - temp;
        K((2*node2-1):(2*node2), (2*node1-1):(2*node1))...
        = K((2*node2-1):(2*node2), (2*node1-1):(2*node1)) - temp;
    end

    function sigma = compute_stress(Q, E, L, theta, node1, node2)
        c = cosd(theta); s = sind(theta);
        sigma = E/L*[-c -s c s]*[Q(2*node1-1); Q(2*node1); Q(2*node2-1); Q(2*node2)];
    end
end