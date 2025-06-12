%++++++++++++++++++++++++++++++++++++++++
% Moog VCF
%
% Chen s2751685
% 11 Feb 2025
%++++++++++++++++++++++++++++++++++++++++

clc
clear all
close all

%++++++++++++++++++++++++++++++++++++++++
% input parameters
SR = 44100 ; % sample rate [Hz]
Tf = 2 ; % total simulation time [s]
f0 = 120 ; % resonant filter frequency [Hz]
r = 0.7 ;  % feedback coeff [choose 0 \leq r \leq 1]

% derived parameters
om0 = 2*pi*f0;          % angular resonant frequency (rad/s)
Nf = floor(Tf*SR);      % total number of samples
k = 1/SR;               % time step

% stability check for **Forward Euler Method**
k_max = 2/(om0*(1 + 2*sqrt(r) - 2*(r^.25)));
assert(k <= k_max, "Error, unstable for the Forward Euler Method")

% initialize
I = eye(4) ;
A = om0*[-1 0 0 -4*r; 1 -1 0 0; 0 1 -1 0; 0 0 1 -1];
b = om0*[1 0 0 0]';
c = [0 0 0 1]';

Bf = I + k*A;   % 4 by 4 matrix used for matrix multiplication in FE
Bb = I - k*A;   % 4 by 4 matrix used for matrix inversion or linear system solution in BE

xf = zeros(4,1);    % initialize 4 by 1 vector for FE state
xb = zeros(4,1);    % initialize 4 by 1 vector for BE state

yf = zeros(Nf, 1);   % initialize the output vector for FE 
yb = zeros(Nf, 1);   % initialize the output vector for BE

u = [1; zeros(Nf-1, 1)];% initialize the input vector as an impulse vector

% tvec = (0:Nf-1)'*k ;     %-- time vector for plots
fvec = (0:Nf-1)'*SR/Nf;  %-- frequency vector for plots

Hc = zeros(Nf,1);                       % continuous-time transfer function 
omega = (0:Nf-1)'*2*pi*SR/Nf;           % angular frequency vector
%++++++++++++++++++++++++++++++++++++++++
% main loop
tic;
for n = 1 : Nf

    % update state xf from n to n+1, with current sample of u as input (FE)
    xf = Bf*xf + k*b*u(n);
    % update state xb from n to n+1, with current sample of u as input (BE)
    xb = Bb \ (xb + k*b*u(n));
    % write sample to output vector yf (FE)
    yf(n) = c' * xf;
    % write sample to output vector yb (BE)
    yb(n) = c' * xb;
    % write sampel to continuous-time transfer funciton (CT)
    Hc(n) = c'*((1i*omega(n)*eye(4) - A)\b);

end
toc ;
%++++++++++++++++++++++++++++++++++++++++

% discrete transfer functions as computed from IR
Hf = fft(yf) ;
Hb = fft(yb) ;

% plot magnitude response of the transfer functions

figure
semilogx(fvec, 20*log10(abs(Hf)))
% findpeaks(20*log10(abs(Hf)), fvec)
title(sprintf('log plot of H(w), r = %.2f, f0 = %d', r, f0))
xlabel('log frequency(Hz)')
ylabel('magnitude dB')

hold on
semilogx(fvec, 20*log10(abs(Hb)))
% findpeaks(20*log10(abs(Hb)), fvec)
xlabel('log frequency(Hz)')
ylabel('magnitude dB')

hold on 
semilogx(fvec, 20*log10(abs(Hc)))
% findpeaks(20*log10(abs(Hc)), fvec)
xlabel('log frequency(Hz)')
ylabel('magnitude dB')

legend('Hf', 'Hb', 'Hc')
xlim([0, SR/2])
ylim([-200, 10])
grid on


