function[V, U,count] = AEF(N, M, Neuron_Str, Iapp)
    % Neuronal parameters for RS, IB and CH 
    RS = [200, 10, -70.0, -50.0, 2.0, 2.0, 30.0, 0.0, -58.0];
    IB = [130, 18, -58, -50.0, 2, 4, 150, 120, -50];
    CH = [200, 10, -58, -50, 2, 2, 120, 100, -46];
    
    Neuron = zeros(N, 9);
    Vs = zeros(N, 1);
    Us = zeros(N, 1);

    for i = 1 : N
        if(Neuron_Str(i) == "RS")
            Neuron(i, :) = RS;
        end
        if(Neuron_Str(i) == "IB")
            Neuron(i, :) = IB;
        end
        if(Neuron_Str(i) == "CH")
            Neuron(i, :) = CH;
        end
        [Vs(i), Us(i)] = solve_ss(Neuron(i, :));
    end
    
    V = ones(N, M);
    V(:, 1) = Vs;
    U = ones(N, M);
    U(:, 1) = Us;
    
    dt = 1e-4;
    count=0;
    for ii = 2 : M
        [V(:, ii), U(:, ii)] = RK1(V(:, ii - 1), U(:, ii - 1), dt, Neuron, Iapp(1, ii));
        for jj = 1 : N
            if(V(jj, ii) >= 0)
                V(jj, ii - 1) = 0;
                V(jj, ii) = Neuron(jj, 9) * 1e-3; % c
                count=count+1;
                U(jj, ii) = U(jj, ii) + Neuron(jj, 8) * 1e-12; % d 
            end
        end
    end

end

function [V_, U_] = RK1(V, U, h, Neuron, Iapp)
    % Euler Method Implemented
    V_ = V + fV(V, U, Iapp, Neuron) * h;
    U_ = U + fU(V, U, Neuron) * h;
end

function der = fV(V, U, Iapp, Neuron)
    % Neuronal parameters
    C = Neuron(:, 1) * 1e-12;
    Gl = Neuron(:, 2) * 1e-9;
    El = Neuron(:, 3) * 1e-3;
    Vt = Neuron(:, 4) * 1e-3;
    del = Neuron(:, 5) * 1e-3;

    % Derivative 
    der = (-1 * Gl .* (V - El) + Gl .* del .* exp((V - Vt) ./ del) - U + Iapp) ./ C;
end


function der = fU(V, U, Neuron)
    % Neuronal parameters
    El = Neuron(:, 3) * 1e-3;
    a = Neuron(:, 6) * 1e-9;
    tao = Neuron(:, 7) * 1e-3;
    
    % Derivative 
    der = (a .* (V - El) - U) ./ tao;
end

function [Vs, Us] = solve_ss(Neuron)
    % Neuronal parameters
    El = Neuron(3) * 1e-3;
    Vt = Neuron(4) * 1e-3;
    gl = Neuron(2) * 1e-9;
    a = Neuron(6) * 1e-9;
    DT = Neuron(5) * 1e-3;

    % Solving a system of linear equations
    syms x y
    [Vss, Uss] = solve(x - El == y / a, gl * (x - El) + y == gl * DT * exp((x - Vt) / DT));
    Vs = double(Vss);
    Us = double(Uss);
end