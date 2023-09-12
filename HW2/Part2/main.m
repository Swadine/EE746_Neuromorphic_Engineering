% Parameters
T = 500; %ms
delt = 0.1;%ms
lambda = 1/ 1000;%per ms

tau = 15;
tau_s = tau/4;
Io = 1e-12; %pA
we=500;

num_steps=5000;

Iapp = zeros(100,num_steps);

Ns=100;
mu = 50;          % Mean
sigma = 5;     % Standard Deviation

weights = mu + sigma * randn(1, Ns);

spike_train=zeros(Ns,num_steps);
%disp((weights));
% Simulate
count=1;
for w=weights
time_instants =[];
for step = 1:num_steps
    random_number = rand();
    event_probability = lambda * delt;
    if random_number < event_probability
        time_instants = [time_instants, step * delt];
        spike_train(count,step)=1;
    end
    %Iapp calc

    for i =time_instants
        Iapp(count,step) =Iapp(count,step)+ Io*w*(exp(-(step*delt-(i))/tau) - exp(-(step*delt-(i))/tau_s));
    end
end
count=count+1;
end

Inew = zeros(1,num_steps);
for i=1:100
    Inew=Inew+Iapp(i,:);
end    


Iapp2 = zeros(100,num_steps);
mu2 = 250;          % Mean
sigma = 50;     % Standard Deviation

weights2 = mu2 + sigma2 * randn(1, Ns);
count=1;
for w2=weights2
time_instants =[];
for step = 1:num_steps
   if spike_train(count,step)==1
        time_instants = [time_instants, step * delt];
    end
    %Iapp calc

    for i =time_instants
        Iapp2(count,step) =Iapp2(count,step)+ Io*w2*(exp(-(step*delt-(i))/tau) - exp(-(step*delt-(i))/tau_s));
    end
end
count=count+1;
end

Inew2 = zeros(1,num_steps);
for i=1:100
    Inew2=Inew2+Iapp2(i,:);
end  

%feed in AEF
Neuron_Str="RS";
N=1;
M=5000;

[V,U,num_spikes]=AEF(1,M,Neuron_Str,Inew);
[V2,U2,num_spikes2]=AEF(1,M,Neuron_Str,Inew2);

fprintf('Number of spikes=');
disp(num_spikes);

fprintf('Number of spikes (b)=');
disp(num_spikes2);

t=1:1:num_steps;

tiledlayout(4,1)
nexttile
figure(1)
plot(t/10,V)
title('Membrane Potential')

nexttile
plot(t/10,Inew)
title('Iapp')

nexttile
plot(t/10,Inew2)
title('Iapp2')

nexttile
plot(t/10,V2)
title('Membrane Potential')
%nexttile
%for i=1:length(spike_train(:,1))
%plot(spike_train(i,:))
%title('spikes')
%hold on
%end