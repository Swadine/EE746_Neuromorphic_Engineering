function [V_, U_] = Izhikevich(N,M,Neuron_Str,Iapp)
    RS = [100.0, 0.7, -60.0, -40.0, 0.03, -2.0, -50.0, 100.0, 35.0];
    IB = [150, 1.2, -75, -45, 0.01, 5, -56, 130, 50];
    CH = [50, 1.5, -60, -40, 0.03, 1, -40, 150, 25];
    Neuron =zeros(N,9);
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
    end
    
    V_ss=(Neuron(:,4)+Neuron(:,6)./Neuron(:,2))*1e-3;
    U_ss=Neuron(:,6) .* 1e-12 .* (Neuron(:,4)+Neuron(:,6)./Neuron(:,2)-Neuron(:,3));
    V = ones(N,M);
    U = ones(N,M);
   
    
    V(:,1)=V_ss;
    U(:,1)=U_ss;
    dt = 1e-4;
    %Iapp = [400; 500; 600]*1e-12;
    
    %[V(:,1), U(:,1)] = RK4_VU(V(:,1), U(:,1), dt, Neuron, Iapp);
    for i = 2:M
        [V(:,i), U(:,i)] = RK4_VU(V(:,i-1), U(:,i-1), dt, Neuron, Iapp);
        for j = 1:N
           if(V(j,i) >= Neuron(j,9)*1e-3)
                V(j,i-1)= Neuron(j,9)*1e-3;
                V(j,i) = Neuron(j,7)*1e-3; % c
                U(j,i) = U(j,i) + Neuron(j,8)*1e-12; % d   
           end
        end
    end
    V_=V;
    U_=U;
end

function [V_, U_] = RK4_VU(V, U, h, Neuron, Iapp)
    V1 = fV(V, U, Iapp, Neuron);
    V2 = fV(V+V1*(h/2), U, Iapp, Neuron);
    V3 = fV(V+V2*(h/2), U, Iapp, Neuron);
    V4 = fV(V+V3*h, U, Iapp, Neuron);
    V_ = V + h*((1/6)*V1 + (1/3)*V2 + (1/3)*V3 + (1/6)*V4);
    
    U1 = fU(V_, U, Neuron);
    U2 = fU(V_, U+U1*(h/2), Neuron);
    U3 = fU(V_, U+U2*(h/2), Neuron);
    U4 = fU(V_, U+U*h, Neuron);
 
    U_ = U + h*((1/6)*U1 + (1/3)*U2 + (1/3)*U3 + (1/6)*U4);
end

function der = fV(V, U, Iapp, Neuron)
C_ = (Neuron(:,1)*1e-12);
kz = Neuron(:,2)*1e-6;
Er = Neuron(:,3)*1e-3;
Et = Neuron(:,4)*1e-3;

%(V-Er).*(V-Et).*kz - U + Iapp

%x = kz .* (V-Er).*(V-Et);
der = 1./C_ .* ( kz.*((V-Er).*(V-Et)) - U + Iapp);
%disp("hi")
%der = V-Et
end


function der = fU(V, U, Neuron)
Er = Neuron(:,3)*1e-3;
a = Neuron(:,5)*1e3;
b = Neuron(:,6)*1e-9;
%c = Neuron(7)*1e-3;
%d = Neuron(8)*1e-12;
%v_peak = Neuron(9)*1e-3;

der = a.* ( b.*(V-Er) - U );
end