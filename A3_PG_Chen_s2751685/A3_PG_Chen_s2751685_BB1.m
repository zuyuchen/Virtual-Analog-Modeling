%++++++++++++++++++++++++++++++++++++++++
% Moog VCF
% BB1 - Calculating and plotting the transfer functions 
% for the numerical schemes
% 
% Chen s2751685
% 18 Feb 2025
%++++++++++++++++++++++++++++++++++++++++


clc
clear all
close all

%++++++++++++++++++++++++++++++++++++++++
% input parameters

SR = 44100;% sample rate
f0 = 120 ; % resonant filter frequency [Hz]
r = 0.7 ;  % feedback coeff [choose 0 \leq r \leq 1]

% derived parameters
om0 = 2*pi*f0;          % angular resonant frequency (rad/s)
k = 1/SR;               % time step

% stability check for **Forward Euler Method**
k_max = 2/(om0*(1 + 2*sqrt(r) - 2*(r^.25)));
assert(k <= k_max, "Error, unstable for the Forward Euler Method")

% initialize
I = eye(4) ;
A = om0*[-1 0 0 -4*r; 1 -1 0 0; 0 1 -1 0; 0 0 1 -1];
b = om0*[1 0 0 0]';
c = [0 0 0 1]';

% prepare for plotting
N = 10000;                              % number of frequency values
fmax = SR/2;                            % maximum frequency to compute

f = (0:N-1)'*fmax/N;                    % frequency vector
omega = (0:N-1)'*2*pi*fmax/N;           % angular frequency vector

Hf = zeros(N,1);                        % Forward Euler transfer function
Hb = zeros(N,1);                        % Backwar Euler transfer function
Hc = zeros(N,1);                        % Continuous-time transfer function
%++++++++++++++++++++++++++++++++++++++++
% main loop
tic;
for n = 1 : N
    z = exp(1i*omega(n)*k);
    Hf(n) = c'*(((z - 1)*eye(4) - k*A)\(k*b));
    Hb(n) = c'*(((z - 1)*eye(4) - k*z*A)\(k*z*b));
    Hc(n) = c'*((1i*omega(n)*eye(4) - A)\b);
end
toc ;

% plot magnitude response of the transfer functions

figure
semilogx(f, 20*log10(abs(Hf)))
title(sprintf('log plot of H(w), r = %.2f, f0 = %d', r, f0))
xlabel('log frequency(Hz)')
ylabel('magnitude dB')

hold on
semilogx(f, 20*log10(abs(Hb)))

hold on 
semilogx(f, 20*log10(abs(Hc)))

legend('Hf', 'Hb', 'Hc')
ylim([-200, 10])
grid on