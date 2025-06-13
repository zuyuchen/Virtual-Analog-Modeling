  
% -------------------------------
% PBMMI - State-space VA Modeling Assignment BB1
% Assymmetric Clipping with damped and capped Newton's methods
%
% Zuyu Chen, 14 April, 2025
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

% Input - create sine wave input
f0 = 500;                               % frequency [Hz]
amp = 2;                                % amplitude of input
tvec = dur*(0:Nf-1)'/Nf;                % time vector [s]
u = amp*sin(2*pi*f0*tvec);              % input vector

% Physical parameters
r = 1e3;                                % resistance of R [Ohms]
c = 33e-9;                              % capacitance of C [F]
Is = 2.52e-9;                           % diode saturation current (A)
Vt = 25.83e-3;                          % diode thermal voltage (V)
Ni = 1.752;                             % diode ideality factor

% Newton-raphson paratemers
tol = 1e-9;                             % convergence tolerance
max_iter = 50;                          % maximum allowed iterations per time-step

% Derived parameteres
k = 1/SR;                               % time interval

% Scheme parameters
A =  -1/(r*c);
B = 1/(r*c);
C = -1/c;
D = 1;
% E = 0;
F = 0;

H_m = (2/k)*1 + A;
H_p = (2/k)*1 - A;
K = D*C/H_m + F;

% Nonlinear function and its derivative
if asymmetric
    g = @(v, p) Is*K*(exp(v/(2*Ni*Vt)) - exp(-v/(Ni*Vt))) - v + p;
    dg_dv = @(v) Is*K*( (1/(2*Ni*Vt))*exp(v/(2*Ni*Vt)) + (1/(Ni*Vt))*exp(-v/(Ni*Vt)) ) - 1;
else
    g = @(v, p) 2*Is*K*sinh(v/(Vt*Ni)) - v + p;         
    dg_dv = @(v) (2 * Is * K / (Vt * Ni)) * cosh(v / (Vt * Ni)) - 1;% Q1) TODO: complete this function to calculate dg/dv for a given v
end
%% function of Newton Raphson method 
function [v, guesses, residual, iter] = newtonRaphsonMethod(p, method, g, dg_dv)
 
    % Newton's method settings
    tol = 1e-10;                    % convergence tolerance
    max_iter = 100;                 % max number of iterations
    guesses = zeros(max_iter,1);    % array to store each individual guess
    residual = guesses;             % array to store residual at each guess
    v_cap = .5e-1;                  % step size cap

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


%% Sinusoid Test
% Initialization 
v = 0;                      % initial guess
y = zeros(Nf, 1);           % ouput vector

% Main loop: compute the output
% Dampped Newton's method
iter_damped = zeros(Nf-1, 1);
disp("damped")
tic
for n = 1:Nf-1

    if asymmetric
        f = Is*(exp(v/(2*Ni*Vt)) - exp(-v/(Ni*Vt))); % Nonliear Function for Asymmetric Clipping
    else
        f = 2*Is*sinh(v/(Vt*Ni));  % symmetric clipping (odd harmonics only)
    end
    p = (H_p*v + B*(u(n) + u(n+1)) + C*f)/H_m;
    [v, guesses, residual, iter] = newtonRaphsonMethod(p, 'damped', g, dg_dv);
    iter_damped(n) = iter;   % accummulate iterations for each sample
    y(n+1) = v; % write to output

end
toc

avg_iter_damped = mean(iter_damped);   % average number of iterations per sample
max_iter_damped = max(iter_damped);    % global maximum of iterations

% Plot the residual and the guesses for computing the last sample output
figure
subplot(2,1,1)
plot(guesses)
xlabel('Iteration no.')
ylabel('Estimate for v')
xlim([1 iter])
title("damped: Residual 'g' vs Guesses 'v'")

subplot(2,1,2)
plot(residual)
xlabel('Iteration no.')
ylabel('Residual')
xlim([1 iter])
% Capped Newton's method
iter_capped = zeros(Nf-1, 1);
disp("capped")
tic
for n = 1:Nf-1

    if asymmetric
        f = Is*(exp(v/(2*Ni*Vt)) - exp(-v/(Ni*Vt))); % Nonliear Function for Asymmetric Clipping
    else
        f = 2*Is*sinh(v/(Vt*Ni));  % symmetric clipping (odd harmonics only)
    end
    p = (H_p*v + B*(u(n) + u(n+1)) + C*f)/H_m;
    [v, guesses, residual, iter] = newtonRaphsonMethod(p, 'capped', g, dg_dv);
    iter_capped(n) = iter;
    y(n+1) = v; % write to output

end
toc

avg_iter_capped = mean(iter_capped);   % average number of iterations per sample
max_iter_capped = max(iter_capped);    % global maximum 

% Plot the residual and the guesses for computing the last sample output
figure
subplot(2,1,1)
plot(guesses)
xlabel('Iteration no.')
ylabel('Estimate for v')
xlim([1 iter])
title("capped: Residual 'g' vs Guesses 'v'")

subplot(2,1,2)
plot(residual)
xlabel('Iteration no.')
ylabel('Residual')
xlim([1 iter])

fprintf("Avegerage number of iterations per sample: damped: %.2f, capped: %.2f", avg_iter_damped, avg_iter_capped);
fprintf("\nGlobal maximum iterations: damped: %.2f, capped: %.2f", max_iter_damped, max_iter_capped);

if plotting
    figure
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

end

% audio out
if audio
    soundsc(y, SR)
    pause(3)
end
%% Audio File Test

[u, Fs] = audioread('A6_PG_Chen_s2751685_Input.wav');
u = mean(u, 2); % averaging the stereo to mono
u = 10*(u/max(abs(u))); % peak normalized and amplified by 15 times
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
    [v, guesses, residual, iter] = newtonRaphsonMethod(p, method, g, dg_dv);
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
        title(sprintf('Assymetric Clipping of Normalized Input with amp = %d', 10))
    else
        title(sprintf('Symmetric Clipping of Normalized Input with amp = %d', 10))
    end
end

if audio
    player = audioplayer(y, Fs);
    play(player);
    % stop(player);
    % audiowrite('A6_PG_Chen_s2751685_Output_Asymmetric.wav', y, Fs)
    % audiowrite('A6_PG_Chen_s2751685_Input.wav', u, Fs)
end
