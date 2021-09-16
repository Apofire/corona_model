
N = 1.3e9; % total population
%y0 = [N-1 1 0 0 0];
load('y4.mat')
y0 = [y4(11,1) y4(11,2) y4(11,3) y4(11,4) y4(11,5)];
%t = linspace(106,145,39);
t = linspace(107,199,93);
%alpha = 0.035; % death rate in India (probability in percentage)
%alpha = -1.8594e-8*t^3 + 1.7153e-6*t^2 + 1.7693e-05*t + 0.034
T_r = 10; % number of days to recover after infection
gamma = 1/T_r;%1/T_r; % recovery rate
%r_exp = (5.1968.*exp(-0.011.*t))
T_inc = 5; % incubation period (average is 5-6 days but can extend to 14 days)
delta = 1/T_inc; % 1/incubation period (number of infections per day)

T_d = 13; % number of days to die after infection
rho = 1/T_d; % 1/number of days to die after infection (number of deaths per day)

[t,y] = ode45(@(t,y) odefun(t,y,alpha,gamma,delta,rho,N,R0),t,y0);
y_end = round(y,0);
save y_end;
%load('covidcases.mat')
% S = y(1) Susceptibile
% E = y(2) Exposed
% I = y(3) Infected
% R = y(4) Recovered
% D = y(5) Dead

function dydt = odefun(t,y,alpha,gamma,delta,rho,N,R0)
   % R0 = -0.0001.*(t.^2) - 0.0176.*t  + 4.6766;
   if t<145
    dydt = zeros(5,1);
    dydt(1) = - (5.1968.*exp(-0.011.*t))*gamma * y(1) * y(3)/N ;
    dydt(2) = (5.1968.*exp(-0.011.*t))*gamma * y(1) * y(3)/N  - (delta * y(2));
    dydt(3) = delta * y(2) - (1-(-1.8594e-8*t^3 + 1.7153e-6*t^2 + 1.7693e-05*t + 0.034)) * (gamma) * y(3) - (-1.8594e-8*t^3 + 1.7153e-6*t^2 + 1.7693e-05*t + 0.034) * rho * y(3);
    dydt(4) = (1--1.8594e-8*t^3 + 1.7153e-6*t^2 + 1.7693e-05*t + 0.034) * (gamma) * y(3);
    dydt(5) = (-1.8594e-8*t^3 + 1.7153e-6*t^2 + 1.7693e-05*t + 0.034) * rho * y(3);
   else
     dydt = zeros(5,1);
     dydt(1) = - (5.1968.*exp(-0.011.*t))*gamma * y(1) * y(3)/N ;
     dydt(2) = (5.1968.*exp(-0.011.*t))*gamma * y(1) * y(3)/N  - (delta * y(2));
     dydt(3) = delta * y(2) - (1-(0.0159)) * (gamma) * y(3) - (0.0159) * rho * y(3);
     dydt(4) = (1-0.0159) * (gamma) * y(3);
     dydt(5) = (0.0159) * rho * y(3);
   end
end
