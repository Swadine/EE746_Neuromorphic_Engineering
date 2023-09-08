N=3;
M=5000;
Iapp = [400; 500; 600]*1e-12;
Neuron_Str=["CH","CH","CH"];
[V,U]=Izhikevich(N,M,Neuron_Str,Iapp);
t=1:1:M;
tiledlayout(2,1)
nexttile
figure(1)
plot(t,V(3,:))
title(Neuron_Str(1)+","+Iapp(3)*1e12+"pA")
ylabel("Membrane Potential")
nexttile
plot(t,U(3,:))
ylabel("Membrane Potential Energy")
