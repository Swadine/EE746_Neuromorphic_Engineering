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
mu = 200;          % Mean
sigma = 20;     % Standard Deviation

weights = mu + sigma * randn(1, Ns);

spike_train1=zeros(Ns,num_steps);

lr =1;


count=1;
for w=weights
time_instants =[];
for step = 1:num_steps
    random_number = rand();
    event_probability = lambda * delt;
    if random_number < event_probability
        time_instants = [time_instants, step * delt];
        spike_train1(count,step)=1;
    end
    %Iapp calc

    for i =time_instants
        Iapp(count,step) =Iapp(count,step)+ Io*w*(exp(-(step*delt-(i))/tau) - exp(-(step*delt-(i))/tau_s));
    end
end
count=count+1;
end

Inew1 = zeros(1,num_steps);
for i=1:100
    Inew1=Inew1+Iapp(i,:);
end 

spike_train2=zeros(Ns,num_steps);
count=1;
Iapp2 = zeros(100,num_steps);
count=1;
for w=weights
time_instants =[];
for step = 1:num_steps
    random_number = rand();
    event_probability = lambda * delt;
    if random_number < event_probability
        time_instants = [time_instants, step * delt];
        spike_train2(count,step)=1;
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
    Inew2=Inew1+Iapp2(i,:);
end 

spike_time1=[];
spike_time2=[];
Neuron = "RS";
[V1, U1,N1,spike_time1] = AEF(1,num_steps, Neuron, Inew1);
[V2, U2,N2,spike_time2] = AEF(1,num_steps, Neuron, Inew2);

%plot a part
t=1:1:5000;
figure
tiledlayout(2,1)
nexttile
plot(t/10, V1)
title("Membrane potential S1")
nexttile
plot(t/10, V2)
title("Membrane potential S2")


% b part remove S2 spikes
%initialise array of arrays
dwk={};
dtk={};
Num_spikes2=N2;
iterations=0;
weights2=weights;

while Num_spikes2>0
%for loops=1:10

for time=spike_time1

t_max=time;
%t_max
tk = t_max * ones(Ns, 1);

for i = 1:Ns
    for ii = t_max:-1:1
        if(spike_train2(i, ii) == 1)
            tk(i, 1) = ii;
            break;
        end
    end
end

del_tk = t_max - tk;
dtk=[dtk,{del_tk}];

del_wk = weights2 .* (exp(-del_tk*delt / tau) - exp(-del_tk*delt / tau_s)) * lr;
dwk=[dwk,{del_wk}];

for i=1:Ns
    weights2(i)=weights2(i) - del_wk(i);
end

for i=1:length(weights2)
    if weights2(i)<10
        weights2(i)=10;
    end
end

Iapp2=zeros(Ns,num_steps);
count=1;
for w=weights2
time_instants =[];
for step = 1:num_steps
   if spike_train2(count,step)==1
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

spt=[];
Neuron = "RS";
[VS2, U_,Num_spikes2,spt] = AEF(1, num_steps, Neuron, Inew2);

iterations=iterations+1;

if Num_spikes2==0
    break;
end

end

end



w1=weights2;


%c part
[VS1, U_,NumS1,spt] = AEF(1, num_steps, Neuron, Inew2);
if NumS1>0
    w2=w1;
else
Num_spikes2=0;
iterations=0;
w2=w1;
while Num_spikes2==0
[V_max, t_max] = max(V1);
%t_max = t_max * delt
%t_max*delt
tk = t_max * ones(Ns, 1);

for i = 1:Ns
    for ii = t_max:-1:1
        if(spike_train1(i, ii) == 1)
            tk(i, 1) = ii;
            break;
        end
    end
end

del_tk = t_max - tk;


del_wk = w2 .* (exp(-del_tk*delt / tau) - exp(-del_tk*delt / tau_s)) * lr;


for i=1:Ns
    w2(i)=w2(i)+del_wk(i);
end
Iapp2=zeros(Ns,num_steps);
count=1;
for w=w2
time_instants =[];
for step = 1:num_steps
   if spike_train1(count,step)==1
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

end

figure
tiledlayout(2,1)
nexttile
plot(t/10, VS2)
title("Vs2")
hold on 
plot(t/10, V2)

nexttile
plot(t/10, V_)
title("Vs1")
hold on 
plot(t/10, V1)


figure
tiledlayout(1,1)
nexttile
plot(weights)
hold on
plot(w2)
legend('old weights','new weights','Location','northeast')


w2

%dpart

% remove S1 spikes
%initialise array of arrays

Num_spikes2=N1;
iterations=0;
weights2=weights;

while Num_spikes2>0
%for loops=1:10

for time=spike_time1

t_max=time;
%t_max
tk = t_max * ones(Ns, 1);

for i = 1:Ns
    for ii = t_max:-1:1
        if(spike_train1(i, ii) == 1)
            tk(i, 1) = ii;
            break;
        end
    end
end

del_tk = t_max - tk;


del_wk = weights2 .* (exp(-del_tk*delt / tau) - exp(-del_tk*delt / tau_s)) * lr;


for i=1:Ns
    weights2(i)=weights2(i) - del_wk(i);
end

for i=1:length(weights2)
    if weights2(i)<10
        weights2(i)=10;
    end
end

Iapp2=zeros(Ns,num_steps);
count=1;
for w=weights2
time_instants =[];
for step = 1:num_steps
   if spike_train1(count,step)==1
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

spt=[];
Neuron = "RS";
[VS1, U_,Num_spikes2,spt] = AEF(1, num_steps, Neuron, Inew2);

iterations=iterations+1;

if Num_spikes2==0
    break;
end

end

end



w1=weights2;


%c part
[VS1, U_,NumS1,spt] = AEF(1, num_steps, Neuron, Inew2);
if NumS1>0
    w2=w1;
else
Num_spikes2=0;
iterations=0;
w2=w1;
while Num_spikes2==0
[V_max, t_max] = max(V1);
%t_max = t_max * delt
%t_max*delt
tk = t_max * ones(Ns, 1);

for i = 1:Ns
    for ii = t_max:-1:1
        if(spike_train2(i, ii) == 1)
            tk(i, 1) = ii;
            break;
        end
    end
end

del_tk = t_max - tk;


del_wk = w2 .* (exp(-del_tk*delt / tau) - exp(-del_tk*delt / tau_s)) * lr;


for i=1:Ns
    w2(i)=w2(i)+del_wk(i);
end
Iapp2=zeros(Ns,num_steps);
count=1;
for w=w2
time_instants =[];
for step = 1:num_steps
   if spike_train2(count,step)==1
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

end

w2

figure
tiledlayout(2,1)
nexttile
plot(t/10, VS1)
title("Vs1")
hold on 
plot(t/10, V1)

nexttile
plot(t/10, V_)
title("Vs2")
hold on 
plot(t/10, V2)


figure
tiledlayout(1,1)
nexttile
plot(weights)
hold on
plot(w2)
legend('old weights','new weights','Location','northeast')
