function[V,U] = AEF(N,M,Neuron_Str,Iapp)
RS = [200, 10, -70.0, -50.0, 2.0, 2.0, 30.0, 0.0, -58.0];
IB = [130, 18, -58, -50.0, 2, 4, 150, 120, -50];
CH = [200, 10, -58, -50, 2, 2, 120, 100, -46];

Neuron =zeros(N,9);
Vs = zeros(N,1);
Us = zeros(N,1);
for i=1:N
    if(Neuron_Str(i)=="RS")
        Neuron(i,:) = RS;
    end
    if(Neuron_Str(i)=="IB")
        Neuron(i,:) = IB;
    end
    if(Neuron_Str(i)=="CH")
        Neuron(i,:) = CH;
    end
    [Vs(i), Us(i)] = solve_ss(Neuron(i,:));
end
Vs
Us
V = ones(N,M);
V(:, 1) = Vs;
U = ones(N,M);
U(:, 1) = Us;

dt = 1e-4;
%Iapp = [250; 350; 450]*1e-12;

for i = 2:M

    [V(:,i), U(:,i)] = RK1(V(:,i-1), U(:,i-1), dt, Neuron, Iapp);
 
    for j = 1:N
        if(V(j,i) >= 0)
            V(j,i-1)= 0;
            V(j,i) = Neuron(j,9)*1e-3; % c
           
            U(j,i) = U(j,i) + Neuron(j,8)*1e-12; % d 
            
        end

    end

        
   
end
end
function [V_, U_] = RK1(V, U, h, Neuron, Iapp)
    V_ = V + fV(V,U,Iapp,Neuron)*h;
    U_ = U + fU(V,U,Neuron)*h;
end

function der = fV(V, U, Iapp, Neuron)
C = Neuron(:,1)*1e-12;
Gl = Neuron(:,2)*1e-9;
El = Neuron(:,3)*1e-3;
Vt = Neuron(:,4)*1e-3;
del= Neuron(:,5)*1e-3;
der = (-1*Gl.*(V-El)+Gl.*del.*exp((V-Vt)./del)-U+Iapp)./C;

end


function der = fU(V, U, Neuron)
El = Neuron(:,3)*1e-3;
a = Neuron(:,6)*1e-9;
tao= Neuron(:,7)*1e-3;

der = (a.*(V-El)-U)./tao;
end

function [Vs,Us] = solve_ss(Neuron)
syms x y
El = Neuron(3)*1e-3;
Vt = Neuron(4)*1e-3;
gl = Neuron(2)*1e-9;
a = Neuron(6)*1e-9;
DT = Neuron(5)*1e-3;
[Vss, Uss] = solve(x - El == y/a, gl*(x - El) + y == gl*DT*exp((x - Vt)/DT) );
Vs=double(Vss);
Us=double(Uss);
end