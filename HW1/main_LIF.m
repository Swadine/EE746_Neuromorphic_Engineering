N=10;
M=5000;
[V,I,avg]=LIF(N,M);
t=1:1:M;
tiledlayout(5,1)
nexttile
figure(1)

plot(t,V(2,:))
title("Neuron 2")
nexttile
plot(t,V(4,:))
title("Neuron 4")
nexttile
plot(t,V(6,:))
title("Neuron 6")
nexttile
plot(t,V(8,:))
title("Neuron 8")

figure(2)
plot(I(:,1),avg,'-o')
title("Average time interval between spikes vs Iapp,k")