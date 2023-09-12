ms = 1e-3;
T = 500 * ms;
dt = 0.1 * ms;
lambda = 1; % spikes / second
P_spike_dt = lambda * dt; % 10 * 0.1ms = 0.001
tau_m = 15 * ms;
tau_s = tau_m / 4;
I_0 = 1e-12; % 1 pA
%w_0 = 250;
%sigma_w = 50;
w_0 = 50;
sigma_w = 5;
Ns = 100;
N = 5; % Random number generation till N decimal digits
lr = 2;

synapses = randn(Ns,1) * sigma_w + w_0;
spikes = zeros(Ns, floor(T / dt));
%spike_times = ones(Ns,t);
Iapp = zeros(1, floor(T / dt));

time_scale = dt:dt:T;

for  i = 1:Ns
    for ii = time_scale
        r = round(rand, N);
        if(r <= P_spike_dt)
            spikes(i, floor(ii / dt)) = 1;
            %spike_times(i,j) = i1;
        end 
           
        x = 0;
        for jj = dt:dt:ii
            if(spikes(i, floor(jj / dt)) == 1)
                x = x + exp((jj - ii) / tau_m) - exp((jj - ii) / tau_s);
            end
        end

        Iapp(1, floor(ii / dt)) = Iapp(1, floor(ii / dt)) + I_0 * synapses(i) * x;
    end
end

Neuron = "RS";
[V, U] = AEF(1, floor(T / dt), Neuron, Iapp);

[V_max, t_max] = max(V);
t_max = t_max * dt

tk = t_max * ones(Ns, 1);

for i = 1:Ns
    for ii = t_max:-dt:dt
        if(spikes(i, round(ii / dt, 0)) == 1)
            tk(i, 1) = ii;
            break;
        end
    end
end

del_tk = t_max - tk;
del_wk = synapses .* (exp(-del_tk / tau_m) - exp(-del_tk / tau_s)) * lr;

for i = 1:Ns
    for ii = time_scale       
        x = 0;
        for jj = dt:dt:ii
            if(spikes(i, floor(jj / dt)) == 1)
                x = x + exp((jj - ii) / tau_m) - exp((jj - ii) / tau_s);
            end
        end
    
        Iapp(1, floor(ii / dt)) = Iapp(1, floor(ii / dt)) + I_0 * del_wk(i, 1) * x;
    end
end

Neuron = "RS";
[V_, U_] = AEF(1, floor(T / dt), Neuron, Iapp);

tiledlayout(2,2)
nexttile
plot(time_scale, V)
title("Membrane potential old")
nexttile
plot(time_scale, V_)
title("Membrane potential")
nexttile
plot(time_scale, Iapp)
title("Iapp")
nexttile
plot(time_scale, spikes)
title("Spikes[Taylor's Version]")




