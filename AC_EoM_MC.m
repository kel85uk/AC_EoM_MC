%
clear all; clc; close all;

matlabpool(4);

Nsamples = 2500; %Number of MC trials
Nsigma = 2; %Sigma confidence levels
Clmax = 1.8; %Max lift coefficient
C_tsfc = 0.00001565; %TSFC
Th_max = 830e3; %Max thrust (N)
coord0 = [6.8 103.52]*pi/180;
X0 = [0 0 35000]; %Expectation of initial position [m m ft]
dX0 = [1e3 1e3 500]; %Uncertainty in initial position [m m ft]
V0 = 474; % Expectation of IAS (knots)
dV0 = 50; % Uncertainty in IAS (knots)
psy0 = 25*pi/180; %Expectation of heading (deg)
dpsy0 = 45*pi/180; %Uncertainty in heading (deg)
hdot0 = 0; %Expectation in climb rate (m/s)
dhdot0 = 5; %Uncertainty in climb rate (m/s)
W0 = 229520; %Expectation of mass of aircraft (kg)
dW0 = 1000; %Uncertainty in mass of aircraft (kg)
S = 427.8; %Expectation of wing area (m^2)
dS = 10; %Uncertainty in wing area (m^2)
Cd0 = 0.026; %Expectation of parasitic drag
dCd0 = 0.005; %Uncertainty in parasitic drag
K = 0.0431; %Expectation of drag polar coefficient
dK = 0.005; %Uncertainty in drag polar coefficient
Cl = 0.5; %Expectation of cruise lift coefficient
dCl = 0.05; %Uncertainty in cruise lift coefficient
%Cd = Cd0 + K*Cl^2
Th = 0; %Expectation of thrust (N)
dTh = 1e3; %Uncertainty in thrust (N)
C = C_tsfc; % Expectation of TSFC
dC = C_tsfc/10; %Uncertainty in TSFC
epsi = 0; %Expectation of thrust line angle (deg)
depsi= 1*pi/180; %Uncertainty in thrust line angle (deg)
mu = 0; %Expectation of roll angle (deg)
dmu = 5*pi/180; %Uncertainty in roll angle (deg)
unknowns_mean = [X0,V0,psy0,hdot0,W0,S,Cd0,K,Cl,Th,C,epsi,mu];
unknowns_var = [dX0,dV0,dpsy0,dhdot0,dW0,dS,dCd0,dK,dCl,dTh,dC,depsi,dmu];
unknowns_var = unknowns_var.^2;
unknowns_max = [2e3,2e3,45000,500,1e7,10,233600,450,0.1,0.1,1.8,10e3,0.2,15*pi/180,40*pi/180];
unknowns_min = [-2e3,-2e3,25000,250,-1e7,-10,100000,400,0.01,0.01,0,0,0,-15*pi/180,-40*pi/180];
% WGS84 parameters (Flat Earth assumption)
a = 6378.137e3; % in meters
f=1/298.257223563; % Shape factor
e2 = f*(2-f);
lat0 = coord0(1);
R1=a*(1-e2)/(1-e2*(sin(lat0))^2)^(3/2);
R2=a/sqrt(1-e2*(sin(lat0))^2);
X_samples = zeros(Nsamples,2);
time_samples = zeros(Nsamples,1);
coord_samples = zeros(Nsamples,3);

parfor NN = 1:Nsamples
	uncertains = rtnorm(unknowns_min,unknowns_max,unknowns_mean,unknowns_var);
	[TE,XE,t,X] = simulator(uncertains);
	%distance_North=R1*dlat;
	dlat = XE(1)/R1;
	dlon = XE(2)/(R2*cos(lat0));
	X_samples(NN,:) = [XE(1),XE(2)];
	%distance_East=R2*cos(lat0)*dlon;
	time_samples(NN) = TE;
	coord = (coord0 + [dlat dlon])*180/pi;
	coord_samples(NN,:) = [coord,0];
	NN
end

matlabpool close;

expect_time = mean(time_samples);
expect_coord = mean(coord_samples,1);
expect_X = mean(X_samples,1);
var_coord = var(coord_samples,0,1);
var_X = var(X_samples,0,1);
sigma_coord2 = Nsigma*sqrt(var_coord);
sigma_X2 = Nsigma*sqrt(var_X);
distance_uncertainty = norm(sigma_X2,2);
coord_write = coord_samples;
coord_write(:,[1 2]) = coord_write(:,[2 1]);
dlmwrite('coord_data.dat',coord_write);
fprintf('Expected crash site = %f N +- %f N, %f E +- %f E with total crash time = % s and radius = %f km \n',expect_coord(1),sigma_coord2(1),expect_coord(2),sigma_coord2(2),expect_time,distance_uncertainty/1000);


