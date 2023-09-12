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

lr =1;


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

Neuron = "RS";
[V, U,Num_spikes] = AEF(1,num_steps, Neuron, Inew);

%initialise array of arrays
dwk={};
dtk={};
Num_spikes2=0;
iterations=0;
while Num_spikes2==0
[V_max, t_max] = max(V);
%t_max = t_max * delt
%t_max*delt
tk = t_max * ones(Ns, 1);

for i = 1:Ns
    for ii = t_max:-1:1
        if(spike_train(i, ii) == 1)
            tk(i, 1) = ii;
            break;
        end
    end
end

del_tk = t_max - tk;
dtk=[dtk,{del_tk}];

del_wk = weights .* (exp(-del_tk*delt / tau) - exp(-del_tk*delt / tau_s)) * lr;
dwk=[dwk,{del_wk}];

for i=1:Ns
    weights(i)=weights(i)+del_wk(i);
end

Iapp2=zeros(Ns,num_steps);
count=1;
for w=weights
time_instants =[];
for step = 1:num_steps
   if spike_train(count,step)==1
        time_instants = [time_instants, step * delt];
    end
    %Iapp calc

    for i =time_instants
        Iapp2(count,step) =Iapp2(count,step)+ Io*w*(exp(-(step*delt-(i))/tau) - exp(-(step*delt-(i))/tau_s));
    end
end
count=count+1;
end

Inew2 = zeros(1,num_steps);
for i=1:100
    Inew2=Inew2+Iapp2(i,:);
end


Neuron = "RS";
[V_, U_,Num_spikes2] = AEF(1, num_steps, Neuron, Inew2);

iterations=iterations+1;
end

fprintf('Tmax=');
disp(t_max*delt);
fprintf('Old number of Spikes=');
disp(Num_spikes);
fprintf('New number of Spikes=');
disp(Num_spikes2);
fprintf('number of iterations=');
disp(iterations);

t=1:1:num_steps;

tiledlayout(3,1)
nexttile
plot(t, V)
title("Membrane potential old")
nexttile
plot(t, V_)
title("Membrane potential")
nexttile
plot(t, Inew)
title("Iapp")

figure
tiledlayout(1,1)
nexttile
for i=1:iterations
scatter(dwk{i},dtk{i}, 'Marker', '.')
hold on
legendLabels{i} = ['iter ', num2str(i)];
end
title("dwk vs dtk")
legend(legendLabels);


