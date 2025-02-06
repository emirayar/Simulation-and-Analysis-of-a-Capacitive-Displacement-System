% Simulation and Analysis of a Capacitive Displacement System 
% with Spring-Mass Dynamics and Electrical Coupling

% Ahmet Arif Tezgel
%200027037

% Emir Ayar
%200030473

clear;
clc;

%% Parameters
epsilon = 8.854 * 10^-12; % Dielectric constant (F/m)
A = 0.0001;          % Capacitor plate area (m^2)
d = 0.18;          % Initial distance between plates (m)
R1 = 0.002;          % Resistance R1 (Ohm)
R2 = 0.001;          % Resistance R2 (Ohm)
Vin = 12;           % Applied voltage (V)
m = 2.228 * 10^-5;  % Mass of the moving plate (kg)
k = 19.9;           % Spring constant (N/m)

%% Initial Conditions
x0 = 0.001;        % Initial displacement (m)
v0 = 0;            % Initial speed (m/s)
Vcap0 = 0;         % Initial capacitor voltage (V)

%% Time Settings
t_end = 0.01;      % Simulation time (s)
dt = 1e-5;         % Time step (s)
t = 0:dt:t_end;    % Time vector

%% Preparation of Variables
N = length(t);     % Number of time steps
x = zeros(N,1);    % Displacement
v = zeros(N,1);    % Velocity
Vcap = zeros(N,1); % Capacitor voltage
Ic = zeros(N,1);   % Current
C = zeros(N,1);    % Capacitance
F_spring = zeros(N,1); % Spring force
F_cap = zeros(N,1); 
Q = zeros(N,1);      % Charge accumulated in the capacitor
E_mech = zeros(N,1);   % Mechanical energy
E_elec = zeros(N,1);   % Electrical energy
dVcap_dt = zeros(N,1);
a = zeros(N,1);
%% Initial Conditions
x(1) = x0;
v(1) = v0;
Vcap(1) = Vcap0;

%% Main Simulation Cycle
for i = 1:N-1
    % Capacitance calculation
    C(i) = epsilon * A / (d + x(i));
    
    % Since the capacitor voltage derivative is a current divider circuit, 
    % VR1 = Vc and I1= Ic In the first step the current is assumed to be zero
    if i > 1
        dVcap_dt(i) = (Vin - R1 * Ic(i)) / (R1 + R2); % Current is 
        % calculated according to the value of the previous step
    else
        dVcap_dt(i) = 0;
    end
    
    % Current calculation (Ic = dQ/dt, Q = C*V)
    Ic(i) = C(i) * dVcap_dt(i); % Capacitor current
    % Total load (Q = integral(Ic) * dt)
    Q(i+1) = Q(i) + Ic(i) * dt;
    
    % Capacitor voltage (Vcap = Q/C)
    %Vcap(i+1) = Vcap(i) + dVcap_dt(i) * dt;
    Vcap(i+1) = Q(i+1) / C(i);

   
    if i > 1
        F_cap(i) = -((Vcap(i)^2) / (2 * d^2)) * ((C(i) - C(i-1)) / (x(i) - x(i-1)));
    else
        F_cap(i) = 0;  % In the first step, F_cap can be assumed zero.
    end


    % Spring force (F = -k*x)
    F_spring(i) = -k * x(i);


    a(i) = (((Vin * (R1 / (R1 + R2)))^2 / 2) * (epsilon * A) ./ (d + x(i)).^4 + -k * x(i)) / m; % Acceleration update
    v(i+1) = v(i) + a(i) * dt; % Speed update
    x(i+1) = x(i) + v(i+1) * dt; % Displacement update
    
    % Energy calculations
    E_mech(i) = 0.5 * m * v(i)^2 + 0.5 * k * x(i)^2; % Kinetic + Potential Energy
    E_elec(i) = 0.5 * C(i) * Vcap(i)^2; % Electrical energy of the capacitor
end

%% Last Step Values
C(end) = epsilon * A / (d + x(end));
F_spring(end) = -k * x(end);
E_mech(end) = 0.5 * m * v(end)^2 + 0.5 * k * x(end)^2;
E_elec(end) = 0.5 * C(end) * Vcap(end)^2;
Vcap(end) = Q(end) / C(end);
F_cap(i) = -((Vcap(end)^2) / (2 * d^2)) * ((C(end) - C(end-1)) / (x(end) - x(end-1)));

%% Graphics
figure;

subplot(4,3,1);
plot(t, x);
title('Displacement (x)');
xlabel('Time (s)');
ylabel('x (mm)');

subplot(4,3,2);
plot(t, v);
title('Velocity (v)');
xlabel('Time (s)');
ylabel('V (mm/s)');

subplot(4,3,3);
plot(t, Vcap);
title('Capacitor Voltage (Vcap)');
xlabel('Time (s)');
ylabel('Vcap (V)');

subplot(4,3,4);
plot(t, Ic);
title('Capacitor Current (Ic)');
xlabel('Time (s)');
ylabel('Ic (A)');

subplot(4,3,5);
plot(t, F_spring);
title('Spring Force (F_spring)');
xlabel('Time (s)');
ylabel('F (N) 10^-3');

subplot(4,3,6);
plot(t, E_mech);
title('Mechanical Energy (E_mech)');
xlabel('Time (s)');
ylabel('Enerji (J)');

subplot(4,3,7);
plot(t, F_cap);
title('Electrical Force (F_cap)');
xlabel('Time (s)');
ylabel('F(N)');

subplot(4,3,8);
plot(t, Q);
title('Total Load on Capacitor (Q)');
xlabel('Time (s)');
ylabel('Q (C)');

subplot(4,3,9);
plot(t, a);
title('Equation of Motion (a)');
xlabel('Time (s)');
ylabel('a (mm/s^2)');

subplot(4,3,10);
plot(t, E_elec);
title('Electrical Energy (E_elec)');
xlabel('Time (s)');
ylabel('Enerji (J)');

