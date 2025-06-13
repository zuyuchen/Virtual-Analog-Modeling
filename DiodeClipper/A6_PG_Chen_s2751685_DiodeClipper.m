% -------------------------------
% PBMMI - State-space VA Modeling Assignment 
% state-space simulation of the Diode Clipper circuit
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

%% function of Newton Raphson method 
% output: v, g, residual, iter
% input: p
function [v, guesses, residual, iter] = newtonRaphsonMethod(p, K)
    % Diode Parameters
    Is = 2.52e-9;                   % diode saturation current (A)
    Vt = 25.83e-3;                  % diode thermal voltage (V)
    Ni = 1.752;                     % ideality factor
    
    % Newton's method settings
    tol = 1e-10;                    % convergence tolerance
    max_iter = 100;                 % max number of iterations
    guesses = zeros(max_iter,1);    % array to store each individual guess
    residual = guesses;             % array to store residual at each guess
    
    
    % Nonlinear function and its derivative
    g = @(v, p) 2*Is*K*sinh(v/(Vt*Ni)) - v + p;         
    dg_dv = @(v) (2 * Is * K / (Vt * Ni)) * cosh(v / (Vt * Ni)) - 1;% Q1) TODO: complete this function to calculate dg/dv for a given v

    % Newton's method
    v = 0;                      % initial guess
    step = 1;                   % initial step value
    iter = 0;                   % iteration counter
    while iter < max_iter && abs(step) > tol %%%% Q2) TODO: define two conditions which must BOTH be true for the while loop to excute
        res =  g(v, p); % Q3) TODO: residual calculation (g(v) = ), use previously defined functions
        J =  dg_dv(v);  % Q4) TODO: jacobian calculation (derivative)
        step = J \ res; % Q5) TODO: Newton step
        v = v - step;               % update current estimate
        iter = iter + 1;
        guesses(iter) = v;          % we can keep track of our guesses to observe convergence
        residual(iter) = res;
    end
end 


%%
% Initialization 
v = 0;                      % initial guess
y = zeros(Nf, 1);           % ouput vector

% Main loop: compute the output
for n = 1:Nf-1
    f = 2*Is*sinh(v/(Vt*Ni));
    p = (H_p*v + B*(u(n) + u(n+1)) + C*f)/H_m;
    [v, guesses, residual, iter] = newtonRaphsonMethod(p, K);
    y(n+1) = v; % write to output
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
    title(sprintf('Symmetric Clipping f0 = %d amp = %d', f0, amp))
    xlim([0, 3/f0])
    
    % Freuqnecy domain plot showing the harmonic distortion
    subplot(2,1,2)
    yfft = 10*log10(abs(fft(y)));
    plot([0:Nf-1]'*(SR/Nf), yfft - max(yfft))
    xlim([0, SR/2]);
    xlabel('$frequency(Hz)$', 'Interpreter', 'latex')
    ylabel('$magnitude$', 'Interpreter', 'latex')
    legend('Normalized Output')
    title(sprintf('Aliased Odd Harmonic Distortion f0 = %d amp = %d', f0, amp))

    % Plot residual and guesses for the last sample
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

%% test wav file

[u, Fs] = audioread('A6_PG_Chen_s2751685_Input.wav');
u = mean(u, 2); % averaging stereo to mono
u = 10*u/max(abs(u));
Nf = length(u);

% Initialization 
v = 0;                             % initial guess
y = zeros(Nf, 1);           % ouput vector

% Main loop: compute the output
for n = 1:Nf-1
    f = 2*Is*sinh(v/(Vt*Ni));
    p = (H_p*v + B*(u(n) + u(n+1)) + C*f)/H_m;
    [v, guesses, residual, iter] = newtonRaphsonMethod(p, K);
    y(n+1) = v; % write to output
end

% Plot input vs output over time to illustrate clipping effect
if plotting
    figure
    tvec = linspace(0, Nf/Fs, Nf);
    plot(tvec, u, 'b')  % input vs time 
    hold on
    plot(tvec, y, 'r')
    legend('u (input)', 'y (output)')
    xlabel('$\\Time(s)$',  'Interpreter', 'latex')
    ylabel('$\\Voltage(v)$', 'Interpreter', 'latex')
    title(sprintf('diode clipping of normalized input with amp = %d', 10));
end

if audio
    player = audioplayer(y, Fs);
    play(player);
    audiowrite('A6_PG_Chen_s2751685_Output.wav', y, Fs)
    % audiowrite('A6_PG_Chen_s2751685_Input.wav', u, Fs)
    % stop(player);
end
