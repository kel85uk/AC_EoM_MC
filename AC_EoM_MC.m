%
clear all; clc; close all;

matlabpool(4);

Nsamples = 10000;

Clmax = 1.8;
C_tsfc = 0.00001565;
Th_max = 830e3;
coord0 = [6.8 103.52]*pi/180;
X0 = [0 0 35000];
dX0 = [1e3 1e3 500];
V0 = 474;
dV0 = 50;
psy0 = 25*pi/180;
dpsy0 = 5*pi/180;
hdot0 = 0;
dhdot0 = 5;
W0 = 229520;
dW0 = 1000;
S = 427.8; %427.8; % m^2
dS = 10;
Cd0 = 0.026; %0.026; %0.5;
dCd0 = 0.005;
K = 0.0431; %0.0431;
dK = 0.005;
Cl = 0.5; %0.4;
dCl = 0.05;
Cd = Cd0 + K*Cl^2
Th = 0;
dTh = 1e3;
C = C_tsfc;%uncertains(13); %0.2;
dC = C_tsfc/10;
epsi = 0; %uncertains(14); %0;
depsi= 1*pi/180;
mu = 0; %uncertains(15); %0;
dmu = 5*pi/180;
unknowns_mean = [X0,V0,psy0,hdot0,W0,S,Cd0,K,Cl,Th,C,epsi,mu];
unknowns_var = [dX0,dV0,dpsy0,dhdot0,dW0,dS,dCd0,dK,dCl,dTh,dC,depsi,dmu];
unknowns_var = unknowns_var.^2;
unknowns_max = [2e3,2e3,45000,500,89*pi/180,10,233600,450,0.1,0.1,1.8,100e3,0.2,10*pi/180,10*pi/180];
unknowns_min = [-2e3,-2e3,25000,250,0,-10,100000,400,0.01,0.01,0,0,0,-10*pi/180,-10*pi/180];

parfor NN = 1:Nsamples
	uncertains = rtnorm(unknowns_min,unknowns_max,unknowns_mean,unknowns_var);
	[TE,XE,t,X] = simulator(uncertains);

	a = 6378.137e3;
	f=1/298.257223563;
	e2 = f*(2-f);
	lat0 = coord0(1);
	R1=a*(1-e2)/(1-e2*(sin(lat0))^2)^(3/2);
	R2=a/sqrt(1-e2*(sin(lat0))^2);
	%distance_North=R1*dlat;
	dlat = XE(1)/R1;
	dlon = XE(2)/(R2*cos(lat0));
	%distance_East=R2*cos(lat0)*dlon;

	coord = (coord0 + [dlat dlon])*180/pi;
	coord_samples(:,:,NN) = coord;
	NN
end

matlabpool close;

expect_coord = mean(coord_samples,3);
var_coord = var(coord_samples,0,3);
dist_var = sqrt(var_coord(1)^2 + var_coord(2)^2);
%R1 and R2 are called the meridional radius of curvature and the radius of curvature in the prime vertical, respectively.

%a is the equatorial radius of the earth (=6378.137000km for WGS84), and e^2=f*(2-f) with the flattening f=1/298.257223563 for WGS84.

%In the spherical model used elsewhere in the Formulary, R1=R2=R, the earth's radius. (using R=1 we get distances in radians, using R=60*180/pi distances are in nm.)


