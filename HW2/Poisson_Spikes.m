ms = 1e-3;
T = 500 * ms;
dt = 0.1 * ms;
lambda = 10; % spikes / second
P_spike_dt = lambda * dt; % 10 * 0.1ms = 0.001
tau_m = 15 * ms;
tau_s = tau_m / 4;
I_0 = 1e-12; % 1 pA
w_e = 500; 
N = 5; % Random number generation till N decimal digits

spikes = zeros(1, floor(T / dt));
spike_times = [];
Iapp = zeros(1, floor(T / dt));

time_scale = dt:dt:T;

for ii = time_scale
    r = round(rand, N);
    if(r <= P_spike_dt)
        spikes(1, floor(ii / dt)) = 1;
        spike_times = [spike_times, ii];
    end 
    
    x = 0;
    for jj = spike_times
        x = x + exp((jj - ii) / tau_m) - exp((jj - ii) / tau_s);
    end

    Iapp(1, floor(ii / dt)) = I_0 * w_e * x;
end

Neuron = "RS";
[V, U] = AEF(1, floor(T / dt), Neuron, Iapp);

tiledlayout(2,2)
nexttile
plot(time_scale, V)
title("Membrane potential")
nexttile
plot(time_scale, Iapp)
title("Iapp")
nexttile
plot(time_scale, spikes)
title("Spikes")





