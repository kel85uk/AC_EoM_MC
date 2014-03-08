function [TE,XE,t,X] = simulator(uncertains)
%	clear all; clc; close all;

	% Initial conditions
	hft         =   uncertains(3); %35000;                       % Altitude above Sea Level, ft
	hm          =   hft * 0.3048;                % Altitude above Sea Level, m

	VKIAS       =   uncertains(4); %474;                         % Indicated Airspeed, kt
	VmsIAS      =   VKIAS * 0.5154;              % Indicated Airspeed, m/s

	[airDens,airPres,temp,soundSpeed] = Atmos(hm);

	qBarFixed   =   0.5*1.225*VmsIAS^2;        % Dynamic Pressure at sea level, N/m^2

	V           =   sqrt(2*qBarFixed/airDens);	% True Airspeed, TAS, m/s	(relative to air mass)
	TASms       =   V;
	hdot    =	uncertains(6); %0;                  % Altitude rate, m/s [Inertial vertical flight path angle = 0 if hdot = 0]
	
%	Initial Conditions Depending on Prior Initial Conditions

	gama	=	57.29578 * atan(hdot / sqrt(V^2 - hdot^2));
						% Inertial Vertical Flight Path Angle, deg
	qbar	= 	0.5 * airDens * V^2;	
						% Dynamic Pressure, N/m^2
	IAS		=	sqrt(2 * qbar / 1.225);
						% Indicated Air Speed, m/s
	Mach	= 	V / soundSpeed;	
						% Mach Number
						
% Aircraft parameters
	S = uncertains(8); %427.8; % m^2
	Cd0 = uncertains(9); %0.5;
	K = uncertains(10);
	Cl = uncertains(11); %0.4;
	Th = uncertains(12); %0;
	C = uncertains(13); %0.2;
	epsi = uncertains(14); %0;
	mu = uncertains(15); %0;
% Simulation parameters
	g = 9.81;
	X0 = [uncertains(1) uncertains(2) hm V uncertains(5) gama uncertains(7)*g];% 233600*g];
	tspan = [0 6000];
% create an options variable
	options = odeset('Events',@event_function);
% note the extra outputs
	[t,X,TE,XE] = ode15s(@EOM,tspan,X0,options,Cd0,K,Cl,S,Th,epsi,mu,C,g);
%	fprintf('Aircraft crashed at x=%f m, y=%f m, z=%f m, t=%f s \n',XE(1),XE(2),XE(3),TE);
end


function [xdot] = EOM(t,X,Cd0,K,Cl,S,Th,epsi,mu,C,g)
	V = X(4);
	x = X(1);
	y = X(2);
	z = X(3);
	psy = X(5);
	gam = X(6);
	W = X(7);

	[rho,P,T,a] = Atmos(z);
	
	Cd = Cd0 + K*Cl^2;

	D = 0.5*rho*Cd*V^2*S;
	L = 0.5*rho*Cl*V^2*S;

	xdot(1) = V*cos(gam)*cos(psy);
	xdot(2) = V*cos(gam)*sin(psy);
	xdot(3) = V*sin(gam);
	xdot(4) = g/W*(Th*cos(epsi) - D - W*sin(gam));
	xdot(5) = g/(W*V*cos(gam))*(Th*sin(epsi) + L)*sin(mu);
	xdot(6) = g/(W*V)*((Th*sin(epsi) + L)*cos(mu) - W*cos(gam));
	xdot(7) = -C*Th;
	xdot = xdot';
end

function [value,isterminal,direction] = event_function(t,X,Cd0,K,Cl,S,Th,epsi,mu,C,g)
% when value is equal to zero, an event is triggered.
% set isterminal to 1 to stop the solver at the first event, or 0 to
% get all the events.
%  direction=0 if all zeros are to be computed (the default), +1 if
%  only zeros where the event function is increasing, and -1 if only
%  zeros where the event function is decreasing.
	value = X(3);  % when value = 0, an event is triggered
	isterminal = 1; % terminate after the first event
	direction = 0;  % get all the zeros
end
