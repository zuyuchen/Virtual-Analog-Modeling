% -------------------------------
% PBMMI - State-space VA Modeling Assignment BB2
% Look-up table with Lagrange cubic inerpolation
%
% Zuyu Chen, 15 April, 2025
% -------------------------------
clc
clear all
close all

% Simulation parameters
SR = 88.2e3;                            % sample rate [Hz]
dur = 2;                                % duration of sim [s]
Nf = round(SR*dur);                     % number of samples in sim

% Settings
plotting = true;                        % plotting on/off
audio = true;                           % play output sounds on/off
method = 'capped';                      % Newton's Method version
asymmetric = true;                      % Assymetric Clipping flag
Lp = -10;                                 % lower boundary of p in the look-up table
Up = 10;                                  % upper boundary of p in the look-up table
Np = 1000;                               % length of the look-up table (number of p-v pairs)
mode = 'cubic';                         % Interpolation mode
Q = 100;                                % Interpolation samples (resolution)

% Input - create sine wave input
f0 = 500;                               % frequency [Hz]
amp = 2;                                % amplitude of input
tvec = dur*(0:Nf-1)'/Nf;                % time vector [s]
u = amp*sin(2*pi*f0*tvec);              % input vector

% Physical parameters
r = [10; 51; 4.7; 1000]*1e3;            % resistance of R [Ohms]
r_dist = 500e3;                        % Overdrive pot 
c = [1; 51e-6; 47e-3]*1e-6;           % capacitance of C [F]
Is = 2.52e-9;                           % diode saturation current (A)
Vt = 25.83e-3;                          % diode thermal voltage (V)
Ni = 1.752;                             % diode ideality factor

% Newton-raphson paratemers
tol = 1e-9;                             % convergence tolerance
max_iter = 50;                          % maximum allowed iterations per time-step

% Derived parameteres
k = 1/SR;                               % time interval

% Scheme parameters
A =  -[1/(r(1)*c(1)) 0 0; 
   1/(r(3)*c(2)) 1/((r(2)+r_dist)*c(2)) 1/(r(3)*c(2)); 
   1/(r(3)*c(3)) 0 1/(r(3)*c(3))];
B = [1/(r(1)*c(1)); 1/(r(3)*c(2)); 1/(r(3)*c(3))];
C = [0; -1/c(2); 0];
D = [0 1 0];
E = 0;
F = 0;
L = [-1 1 0];
M = 1;
N = 0;

H_m = (2/k)*eye(3) + A;
H_p = (2/k)*eye(3) - A;
K = D*(H_m\C) + F;

% Nonlinear function and its derivative
if asymmetric
    g = @(v, p) Is*K*(exp(v/(2*Ni*Vt)) - exp(-v/(Ni*Vt))) - v + p;
    dg_dv = @(v) Is*K*((1/(2*Ni*Vt))*exp(v/(2*Ni*Vt)) + (1/(Ni*Vt))*exp(-v/(Ni*Vt))) - 1;
else
    g = @(v, p) 2*Is*K*sinh(v/(Vt*Ni)) - v + p;         
     dg_dv = @(v) (2 * Is * K / (Vt * Ni)) * cosh(v / (Vt * Ni)) - 1;% Q1) TODO: complete this function to calculate dg/dv for a given v
end

%% precompute the look-up table to store p-v pairs
%
pax = linspace(Lp, Up, Np);        % p-axis (array of p values)
vax = 0*pax;                       % v-axis (array to store solutions)
for n = 1:length(pax)
[v, ~, ~, ~] = newtonRaphsonMethod(pax(n), method, g, dg_dv);
vax(n) = v;
end

% plot p vs v
figure
plot(pax, vax)
xlabel('p')
ylabel('v')
title('p - v relationship')

%% Newton Raphson method 
%
function [v, guesses, residual, iter] = newtonRaphsonMethod(p, method, g, dg_dv)

% Newton's method settings
tol = 1e-10;                    % convergence tolerance
max_iter = 100;                 % max number of iterations
guesses = zeros(max_iter,1);    % array to store each individual guess
residual = guesses;             % array to store residual at each guess
v_cap = .5e-1;                   % step size cap

% Newton's method
v = 0;                      % initial guess
step = 1;                   % initial step value
iter = 0;                   % iteration counter
while iter < max_iter && abs(step) > tol %%%% Q2) TODO: define two conditions which must BOTH be true for the while loop to excute
    res =  g(v, p); % Q3) TODO: residual calculation (g(v) = ), use previously defined functions
    J =  dg_dv(v);  % Q4) TODO: jacobian calculation (derivative)
    step = J \ res; % Q5) TODO: Newton step
    switch method
        case 'basic'
            v = v - step;               % update current estimate
            iter = iter + 1;
            guesses(iter) = v;          % we can keep track of our guesses to observe convergence
            residual(iter) = res;
        case 'damped'

            iter = iter + 1;
            residual(iter) = res;
            rate = 0.5;                 % rate for reducing step size
            next_res = res;             % initialize the next residual value
            
            % subloop for finding the next decreasing residual value
            while abs(next_res) >= abs(res)

                step = rate * step;      % reduce the step size by half
                v_temp = v - step; 
                next_res = g(v_temp, p); % update current residual   

            end

            v = v_temp;
            guesses(iter) = v;          % keep track of the guesses to observe convergence
            res = next_res;             % shift the residual value
         case 'capped'
             if abs(step) > v_cap
                 step = sign(step) * v_cap;
             end
                v = v - step;               % update current estimate
                iter = iter + 1;
                guesses(iter) = v;          % we can keep track of our guesses to observe convergence
                residual(iter) = res;
    end  
end
end 

%% Cubic Lagrange Interpolation
%
function T = lagrangeInterp(P, Q)

m = (0:P-1)';            % Indices of N samples
a_m_vec = -(P-1)/2 + m;  % vector for the N sample locations
T = zeros(Q, P);         % Initialize the table     

% evaluates the polynomials for Q samples in the range of [-1/2, 1/2)  

% create the sample location vector for evaluating the polynomials
a = (-Q/2 + (0:Q-1))/Q ; 
% loop along the columns
for m = 0:P-1
    % the m-th sample location 
    a_m = -(P-1)/2 + m; 
    % sample location vector eliminating the current sample location
    a_q = a_m_vec(a_m_vec ~= a_m);
    % calculates the Lagrange polynomial and assign it to a column
    % in the table
    T(:,m+1) = prod(a - a_q)./prod(a_m - a_q);
end

end

%% Sinusoid Test
% Initialization 
x = zeros(3,1);             % initial state
v = x(2);                   % initial guess
                            
if asymmetric               % compute f(v)
    f1 = Is*(exp(v/(2*Ni*Vt)) - exp(-v/(Ni*Vt))); % Nonliear Function for Asymmetric Clipping
else
    f1 = 2*Is*sinh(v/(Vt*Ni)); % symmetric clipping (odd harmonics only)
end

y = zeros(Nf, 1);           % ouput vector
p = zeros(Nf, 1);
% Main loop: compute the output
for n = 1:Nf-1

    p(n) = D*(H_m\(H_p*x + B*u(n) + C*f1)) + (D*(H_m\B) + E)*u(n+1);
    [v, ~, ~, ~] = newtonRaphsonMethod(p(n), method, g, dg_dv);

    % switch mode
    % 
    %     case 'round'
    %         % just round the p to the nearest p stored in the p-v table and
    %         % use the corresponding v as the approximate v
    %         [~, ind] = min(abs(p - pax));
    %         v = vax(ind); 
    % 
    %     case 'cubic'
    %         P = 4; 
    %         % in case p happens to be a precomputed value
    %         [val, ind] = min(abs(p(n) - pax)); % position of the nearest p stored in the p-v table
    %         if val == 0
    %             v = vax(ind);
    %             break
    %         end
    % 
    %         % if not, apply Lagrange Cubic Interpolation
    %         % =======================================
    % 
    %         T = lagrangeInterp(P, Q);    % Lagrange inteprolation table (Q * 4)
    %         a = (-Q/2 + (0:Q-1))/Q ;     % sample locations
    % 
    %         delta_pax = pax(2) - pax(1);          % interval of pax vector
    % 
    %         p_ind = ((p(n) - Lp)/(Up - Lp)) * Np;    % index of the p in "p-v look-up table" (in the range of 1:Np)
    % 
    %         v_interp = [vax(floor(p_ind)-1); vax(floor(p_ind));   % the P-point interpolation values
    %                     vax(floor(p_ind)+1); vax(floor(p_ind)+2)];
    % 
    %         a_p = p_ind - floor(p_ind) - 0.5;     % the zero-centered fraction ((0.5, 0.5)) of non-integer index 
    % 
    %         [~, row_idx] = min(abs(a_p - a));     % the nearest position in the Lagrange iterpolation samples
    % 
    %         % compute v as the weighted average of the P-point interpolation values
    %         v = T(row_idx, :) * v_interp;
    % 
    % end
    
    % update f(v)
    if asymmetric
        f = Is*(exp(v/(2*Ni*Vt)) - exp(-v/(Ni*Vt))); % Nonliear Function for Asymmetric Clipping
    else
        f = 2*Is*sinh(v/(Vt*Ni));  % symmetric clipping (odd harmonics only)
    end
    
    % state update
    x = H_m \ (H_p*x + B*(u(n) + u(n+1)) + C*(f1 + f));
    
    % write to output     
    y(n+1) = u(n) - x(1) + x(2);    
    
    % shift state
    f1 = f;

end


if plotting

% Plot input vs output over time to illustrate clipping effect
subplot(2,1,1)
plot(tvec, u, 'b')  % input vs time 
hold on
plot(tvec, y, 'r')
legend('u (input)', 'y (output)')
xlabel('$\\Time(s)$',  'Interpreter', 'latex')
ylabel('$\\Voltage(v)$', 'Interpreter', 'latex')
if asymmetric
    title(sprintf('Assymetric Clipping f0 = %d amp = %d', f0, amp))
else
    title(sprintf('Symmetric Clipping f0 = %d amp = %d', f0, amp))
end
xlim([0, 3/f0])

% Freuqnecy domain plot showing the harmonic distortion
subplot(2,1,2)
yfft = 10*log10(abs(fft(y)));
plot([0:Nf-1]'*(SR/Nf), yfft - max(yfft))
xlim([0, SR/2]);
xlabel('$frequency(Hz)$', 'Interpreter', 'latex')
ylabel('$magnitude$', 'Interpreter', 'latex')
legend('Normalized Output')
if asymmetric
    title(sprintf('Aliased All-Integer Harmonic Distortion f0 = %d amp = %d', f0, amp))
else 
    title(sprintf('Aliased Odd Harmonic Distortion f0 = %d amp = %d', f0, amp))
end

% Plot the residual and the guesses for computing the last sample output
figure
subplot(2,1,1)
plot(guesses)
xlabel('Iteration no.')
ylabel('Estimate for v')
xlim([1 iter])
title(" Residual 'g' vs Guesses 'v'")

subplot(2,1,2)
plot(residual)
xlabel('Iteration no.')
ylabel('Residual')
xlim([1 iter])

end

% audio out
if audio
soundsc(y, SR)
pause(3)
end

%% Audio File Test

[u, Fs] = audioread('A6_PG_Chen_s2751685_Input.wav');
u = mean(u, 2); % averaging the stereo to mono
u = 10*(u/max(abs(u))); % peak normalized and amplified by 10 times
Nf = length(u);

% Initialization 
v = 0;                             % initial guess
y = zeros(Nf, 1);           % ouput vector

% Main loop: compute the output
for n = 1:Nf-1
if asymmetric
    f = Is*(exp(v/(2*Ni*Vt)) - exp(-v/(Ni*Vt))); % Nonliear Function for Asymmetric Clipping
else
    f = 2*Is*sinh(v/(Vt*Ni));  % symmetric clipping (odd harmonics only)
end
p = (H_p*v + B*(u(n) + u(n+1)) + C*f)/H_m;
[v, guesses, residual, iter] = newtonRaphsonMethod(p, K, method);
y(n+1) = v; % write to output
end

if plotting 
% Plot input vs output over time to illustrate clipping effect
figure
tvec = linspace(0, Nf/Fs, Nf);
plot(tvec, u, 'b')  % input vs time 
hold on
plot(tvec, y, 'r')
legend('u (input)', 'y (output)')
xlabel('$\\Time(s)$',  'Interpreter', 'latex')
ylabel('$\\Voltage(v)$', 'Interpreter', 'latex')
if asymmetric
    title(sprintf('Assymetric Clipping of Normalized Input f0 = %d amp = %d', f0, amp))
else
    title(sprintf('Symmetric Clipping of Normalized Input f0 = %d amp = %d', f0, amp))
end
end

if audio
player = audioplayer(y, Fs);
play(player);
% stop(player);
% audiowrite('A6_PG_Chen_s2751685_Output_Asymmetric.wav', y, Fs)
% audiowrite('A6_PG_Chen_s2751685_Input.wav', u, Fs)
end
