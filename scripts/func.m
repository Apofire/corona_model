%% Haar Wavelet method to solve system of linear differential equations

                % KAUSHIK IYER
                
%           y1'(t) = y3(t) - cos(t), y1(0) = 1
%           y2'(t) = y3(t) - exp(t), y2(0) = 0
%           y3'(t) = y1(t) - y2(t),  y3(0) = 2
             
%           the exact solution is 
%        y1(t) = exp(t)     y2(t) = sin(t)      y3(t) = exp(t) + cos(t)

addpath '/Users/kaushikiyer/Desktop/haar wavelets/codes/WHM'
% step 1
% collocation points
J    = input('Enter the level of decomposition:');  % level of decomposition
N    = 2^(J+1); % N = 2M                            % number of basis functions
j    = 1:N;                                         % index of grid points
t    = (2.*j - 1) ./ (2.*N);                                % grid points

q0 = 1;
lambda = 0.05;
alpha = 0.003;
lambda1 = 0.0001;
s = 0.01;
L = 1000;
theta = 0.00001;
u = 0.2;
M = 2000;
pie = 0.01;
phi = 0.0002;
pi1 = 0.01;
X0 = 500;
N0 = 270;
F0 = 230;

%generating the Haar matrix and integral of Haar matrix of size (2M x 2M)
H = zeros(N,N);                                     % initialising the matrix with zeros
P = zeros(N,N);
for i = 1:N
    H(i,:) = haar_matrix(t,i,J); 
    P(i,:) = integral_of_H(t,i,J);
end
%disp(H);




% finding the filter coefficients
a_x = zeros(1,length(t));                            % coefficient function matrix at the collocation points
a_n = zeros(1,length(t));
a_f = zeros(1,length(t));
E  = ones(1,length(t)); 

%a3 = fsolve(@(a3) fun(F_t,G_t,E,H,P,a3),zeros(1,length(t)));
% using the initial conditions to find the first estimate of the function
% and using it to find the filter coefficients.
a_x = fsolve(@(a_x) myfun1(E,H,P,a_x,a_n,a_f,q0,lambda,phi,u,M,L,alpha,lambda1,pie,pi1,N0,X0,F0),zeros(1,length(t)));


% constructing the approximated solution 
approx_sol_yx = a_x*P + X0.*E;
approx_sol_yn = a_n*P + N0.*E;
approx_sol_yf = a_f*P + F0.*E;

% true solution
f = @(t,x) [ones(1,length(t)) + 0.05*x(2) - 0.003*x(1) - 0.0001*(x(1).*(x(3)));
    0.01*x(2).*((ones(1,length(t)) - x(2)./1000)) - 0.00001*x(1).*x(2) + 0.01*0.0007*(x(2).*x(3));
    0.2*x(3).*(ones(1,length(t)) - x(3)./2000) - 0.0007*(x(2).*x(3)) + 0.01*0.0001*(x(1).*x(3))];

[t,x] = ode45(f,t,[500 270 230]);


% error
error = zeros(length(t),3);
error(:,1) = abs(x(:,1)' - approx_sol_yx);
error(:,2) = abs(x(:,2)' - approx_sol_yn);
error(:,3) = abs(x(:,3)' - approx_sol_yf);



 %% Plot graphics 
% fig:01
figure('color','w')
plot(t,approx_sol_yx,'g',t,x(:,1),'rs')
hold on
plot(t,approx_sol_yn,'g',t,x(:,2),'rs')
hold on
plot(t,approx_sol_yf,'g',t,x(:,3),'rs')
xlabel('$t$'); ylabel('$u(t)$');
title(['J = ' num2str(J) ', ' '2M = ' num2str(N)])
legend('Wavelet','Exact')


% fig:02
figure('color','w')
plot(t,error,'r.-')
xlabel('$t$'); ylabel('Absolute Error');
title('Absolute Error: $\max|y_{numeric} - y_{analytic}|$')





%% Functions

function y = haar_matrix(t,i,J)
% Function to generate the haar function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the haar function

if i == 1
    m = 0;
    k = 0;
    
else
    IMat = zeros([J+1 2^J]);
    IMask = IMat;

    ind_s = 1;
    for ind_j = 0:J    
        for ind_i = 0:2^ind_j-1
            ind_s = ind_s + 1;
            IMask(ind_j+1,ind_i+1) = ind_s;     
            IMat(ind_j+1,ind_i+1) = ind_i+ind_j+1;
        end % for i
    end % for j

    [ind_j, ind_i] = find(IMask == i);
    m = 2^(ind_j - 1);
    k = ind_i - 1;
    j = ind_j -1;
end

y = zeros([length(t) 1]);

if i == 1    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = 1;
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(t)
        if (alpha <= t(i) && (t(i) < beta))
            y(i) = 1;
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = -1;
        else
            y(i) = 0;
        end
    end
end


end


function y = integral_of_H(t,i,J)
% Function to generate the integral of the haar function for the given interval for a
% given index 'i'

% inputs
% x - the length of the vector containing the collocation points
% i - index of the haar function
% J - Maximum resolution
% outputs
% y = column vector containing the integral of the haar function

if i == 1
    m = 0;
    k = 0;
else
    IMat = zeros([J+1 2^J]);
    IMask = IMat;

    ind_s = 1;
    for ind_j = 0:J    
        for ind_i = 0:2^ind_j-1
            ind_s = ind_s + 1;
            IMask(ind_j+1,ind_i+1) = ind_s;     
            IMat(ind_j+1,ind_i+1) = ind_i+ind_j+1;
        end % for i
    end % for j

    [m, k] = find(IMask == i);
    m = 2^(m - 1);
    k = k - 1;
end


y = zeros([1 length(t)]);

if i == 1    
    for i = 1:length(t)
        if (0 <= t(i) && (t(i) < 1))
            y(i) = t(i);
        else
            y(i) = 0;
        end
    end
else
    alpha = k / m;
    beta = (k + 0.5) / m;
    gamma = (k + 1) / m;
    for i = 1:length(t)
        if (alpha <= t(i) && (t(i) < beta))
            y(i) = t(i) - alpha;
        elseif (beta <= t(i) && (t(i) < gamma))
            y(i) = gamma - t(i);
        else
            y(i) = 0;
        end
    end
end



end

function F = myfun1(E,H,P,a_x,a_n,a_f,q0,lambda,phi,u,M,L,alpha,lambda1,pie,pi1,N0,X0,F0)
F = a_x*H - q0.*E - lambda.*((a_n)*P + N0.*E) - alpha.*(a_x*P + X0.*E) + lambda.*(a_x*P + X0.*E).*(a_f*P + F0.*E);
end


