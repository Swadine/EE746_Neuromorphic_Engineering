% Parameters
T = 500; %ms
delt = 0.1;%ms
lambda = 10 / 1000;%per ms

tau = 15;
tau_s = tau/4;
Io = 1e-12; %pA
we=500;

num_steps=5000;

%time instants
time_instants =[];
Iapp = zeros(1,num_steps);
stimuli = zeros(1,num_steps);

% Simulate
for step = 1:num_steps
    random_number = rand();
    event_probability = lambda * delt;

    if random_number < event_probability
        time_instants = [time_instants, step * delt];
        stimuli(step)=1;
    end
    %Iapp calc
    for i =time_instants
        Iapp(step) =Iapp(step)+ Io*we*(exp(-(step*delt-(i))/tau) - exp(-(step*delt-(i))/tau_s));
    end
end

%disp((Iapp));
fprintf('number of stimuli=:');
disp(length(time_instants));
fprintf('tk(ms):\n');
disp(time_instants);

%feed in AEF
Neuron_Str="RS";
N=1;
M=5000;

t=1:1:num_steps;

[V,U,num]=AEF(1,M,Neuron_Str,Iapp);
fprintf('Number of spikes=');
disp(num);

tiledlayout(3,1)
nexttile
plot(t/10,stimuli)
title('stimuli');

nexttile
plot(t/10,Iapp(1,:))
title('Current');

nexttile
plot(t,V(1,:))
title('Voltage');