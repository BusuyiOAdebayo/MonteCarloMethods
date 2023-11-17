%% Monte Carlo simulation of 1D random walk
% Jixin Chen @ Ohio University started 10/31/2019 driving to Youngstown
% University seminar talk
% Jixin Chen 2020-05-25
% Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free
% Documentation License, Version 1.3 or any later version published by the Free Software Foundation.
clear;
%% calculate diffusion constant
T = 300; % temperature unit K.
R = 8.314; % gas constant unit J/K/mol.
Nav = 6.02e23; % Avogadro constant.
kB = 1.38e-23; % Boltzmann constant.
yita = 8.9e-4; % water viscosity unit Pa*s.
rou = 0.8; % solute density unit g/mL.
Mw = 1271; % molecular weight g/mol.
%----convert to SI unit----
rou = rou*1000; % density unit kg/m^3.
Mw = Mw/1000; % molecular weight kg/mol.
%---
D = kB*T/(6*pi*yita*(3*Mw/4/pi/Nav/rou)^(1/3)); % unit m^2/s.
%% 1D diffusion
% ------- L
% o
% o
% o
% -------- 0
conc = 1e-6; % concentration unit mole/L.
L = 100000; % space between two surfaces unit nm.
area = 1000*2000; % surface area of the surface unit nm^2.
dt = 1; % simulated time step unit ms.
lt = 100; % total simulated time unit s.
%----convert to SI unit -----
conc = conc*1000*Nav; % conc unit m^-3.
distance = conc^(-1/3); % unit m
distime = distance^2/2/D;
L = L*1e-9; % unit m.
area = area*1e-18; % unit m^2.
dt = dt*1e-3; %unit s.
%----
% numm = round(conc*L*area); % number of molecule in the volume.
S17
numm = 10000;
numt = round(lt/dt); % total time step index.
meandis = L/numm; % average space between two molecules in z.
sigma = sqrt(2*D*dt);
% meandisv = conc^(-1/3); % average space between two molecules in volume.
meantime = meandis^2/2/D; % average time to travel the average space.
% check meantime and dt. dt should be much smaller (>10 times) than meantime.
display(['meantime/dt = ', num2str(meantime/dt)]);
display(['simu length = ', num2str(lt/dt)]);
display(['num molecules = ', num2str(numm)]);
%------- initialize -----
mcst = randn(numm, numt); %Monte Carlo space and time. first column not used.
%mcst = mcst*1.57*sigma; % ??????????? why 1.57 ???????????? golden number 1.618, Great pyrimid side
height ratio 1.57
% mcst = mcst2*2*sigma^2; % ??????????? why 1.57 ???????????? golden number 1.618, Great pyrimid side
height ratio 1.57
% mcst = sqrt(abs(mcst)).*mcst2./abs(mcst2);
mcst = mcst*sigma; %
traj = zeros(numm, numt); % trajectories of all molecules.
% traj(:,1) = L*rand(numm, 1); % time 0 all molecules randomly located
traj(:,1) = zeros(numm, 1); % middel space for all initial locations.
%------- random walk at the step time resolution ------
rng('shuffle');
for i = 2:numt
 traj(:, i) = traj(:, i-1) + mcst(:, i); % random walk
% indl = find(traj(:, i)< 0); % hit lower wall
% indu = find(traj(:, i)> L); % hit upper wall
% hitl(indl, i) = 1; traj(indl, i) = -traj(indl, i); % reflect and record hit.
% hitu(indu, i) = 1; traj(indu, i) = 2*L - traj(indu, i); % reflect and record hit.
end
figure; hist(traj(:,2), 100); title('step 1');
hold on;
t = dt;
x = -L*1E9:L*1E9;
x = x*1E-9;
pdf = numm/100*3*exp(-(x).^2/4/D/t);
plot(x, pdf);
%figure; hist(traj(:,3), 100); title('step 2');
figure; hist(traj(:,11), 100); title('step 10');
hold on;
t = dt*10;
x = -L*1E9:L*1E9;;
x = x*1E-9;
pdf = numm/100*3*exp(-(x).^2/4/D/t);
plot(x, pdf);
S18
%figure; hist(traj(:,3), 100); title('step 2');
figure; hist(traj(:,101), 100); title('step 100');
hold on;
t = dt*100;
x = -L*1E9:L*1E9;;
x = x*1E-9;
pdf = numm/100*3*exp(-(x).^2/4/D/t);
plot(x, pdf);
figure; hist(traj(:,1001), 100); title('step 1000');
hold on;
t = dt*1000;
x = -L*1E9:L*1E9;;
x = x*1E-9;
pdf = numm/100*3*exp(-(x).^2/4/D/t);
plot(x, pdf);
figure; hist(traj(:,10001), 100); title('step 10,000');
hold on;
t = dt*10000;
x = -L*1E9:L*1E9;
x = x*1E-8;
pdf = numm/100*3*exp(-(x).^2/4/D/t);
plot(x, pdf);
figure; hist(traj(:,100000), 100); title('step 99999');
hold on;
t = dt*99999;
x = -L*1E9:L*1E9;
x = x*1E-7;
pdf = numm/100*3*exp(-(x).^2/4/D/t);
plot(x, pdf);
