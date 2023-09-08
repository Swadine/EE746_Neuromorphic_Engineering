% Constants
Cm = 1.0;       % Membrane capacitance (uF/cm^2)
gNa = 120;      % Sodium conductance (mS/cm^2)
ENa = 50;      % Sodium reversal potential (mV)
gK = 36;        % Potassium conductance (mS/cm^2)
EK = -77;       % Potassium reversal potential (mV)
gL = 0.3;       % Leak conductance (mS/cm^2)
EL = -55;      % Leak reversal potential (mV)
P = 0;
% Initial conditions
V0 = -65;       % Initial membrane potential (mV)
m0 = 0.05;      % Initial m gate value
h0 = 0.6;       % Initial h gate value
n0 = 0.32;      % Initial n gate value

% Time parameters
T = 30;         % Total simulation time (ms)
dt = 0.01;      % Time step (ms)
t = 0:dt:4*T;     % Time vector

% External current
I0 = 15;        % Amplitude of the current step (uA/cm^2)
t_on = 2*T;
t_off = 3*T;
Iext = I0 * (t >= t_on & t < t_off);

% Initialize vectors to store variables over time
V = zeros(size(t));
m = zeros(size(t));
h = zeros(size(t));
n = zeros(size(t));
INa = zeros(size(t));
IK = zeros(size(t));
IL = zeros(size(t));
PNa = zeros(size(t));
PK = zeros(size(t));
PL = zeros(size(t));

% Set initial values
V(1) = V0;
m(1) = m0;
h(1) = h0;
n(1) = n0;
 
% Simulate the Hodgkin-Huxley model
for i = 1:length(t) - 1
    % Compute gating variable derivatives
    alpha_m = (0.1 * (V(i) + 40)) / (1 - exp(-(V(i) + 40) / 10));
    beta_m = 4 * exp(-(V(i) + 65) / 18);
    alpha_h = 0.07 * exp(-(V(i) + 65) / 20);
    beta_h = 1 / (1 + exp(-(V(i) + 35) / 10));
    alpha_n = (0.01 * (V(i) + 55)) / (1 - exp(-(V(i) + 55) / 10));
    beta_n = 0.125 * exp(-(V(i) + 65) / 80);

    % Update gating variables
    m(i + 1) = m(i) + dt * (alpha_m * (1 - m(i)) - beta_m * m(i));
    h(i + 1) = h(i) + dt * (alpha_h * (1 - h(i)) - beta_h * h(i));
    n(i + 1) = n(i) + dt * (alpha_n * (1 - n(i)) - beta_n * n(i));

    % Compute membrane current
    INa(i+1) = gNa * m(i)^3 * h(i) * (V(i) - ENa);
    IK(i+1) = gK * n(i)^4 * (V(i) - EK);
    IL(i+1) = gL * (V(i) - EL);

    % Update membrane potential
    V(i + 1) = V(i) + dt * ((Iext(i) - INa(i) - IK(i) - IL(i)) / Cm);
end

%instantaneous power
PNa=INa.*(V-ENa);
PK=IK.*(V-EK);
PL=IL.*(V-EL);
P_membrane=((Iext - INa - IK - IL)).*V;

%c part- power for membrane
range=t>60 & t<70;
for i=6000:7000
    P=P+P_membrane(i);
end
-P*1e-15

% Plot the membrane potential with spikes
figure;
tiledlayout(5,1)
nexttile
plot(t, Iext)
title('Hodgkin-Huxley Neuron Model');
nexttile
plot(t, V);
legend('Membrane Potential (mV)','Location','northwest');
nexttile
plot(t, INa)
legend('INa','Location','northwest');
nexttile
plot(t, IK)
legend('IK','Location','northwest');
nexttile
plot(t, IL)
legend('IL','Location','northwest');

figure;
tiledlayout(1,1)
%nexttile
%plot(t(range), Iext(range))
%title('Hodgkin-Huxley Neuron Model');
nexttile
plot(t(range), PNa(range));
legend('PNa','PL','Location','northeast');
hold on
%nexttile
plot(t(range), PK(range))
legend('PK','Location','northeast');
%nexttile
plot(t(range), PL(range))
%nexttile
plot(t(range), P_membrane(range))
legend('PNa','Pk','PL','Pmembrane','Location','northeast');
title('Instantaneous Powers of Channels and Membrane');
ylabel('Power(nW/cm^2)');
hold off
%legend('Pmembrane','Location','northeast');