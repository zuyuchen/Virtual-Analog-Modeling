%++++++++++++++++++++++++++++++++++++++++
% Moog VCF
% BB2 - calculating and plotting the transfer function for both the linear 
% and the non-linear Moog VCF system using trapezoidal integration method 
% for the discrete simulation and Newton-Raphson algorithm for solviing the
% difference equation in the non-linear case. It also includes the plot for 
% a sinusoidal input 
%
% Chen s2751685
% 18 Feb 2025
%++++++++++++++++++++++++++++++++++++++++

clc
clear all
close all

%++++++++++++++++++++++++++++++++++++++++
% input parameters
SR = 44100 ; % sample rate [Hz]
Tf = 2 ; % total simulation time [s]
f0 = 50 ; % resonant filter frequency_NL[Hz]
r = 0.9 ;  % feedback coeff [choose 0 \leq r \leq 1]

% derived parameters
om0 = 2*pi*f0;          % angular resonant frequency_NL(rad/s)
Nf = floor(Tf*SR);      % total number of samples
k = 1/SR;               % time step

% stability_NLcheck for **Forward Euler Method**
% k_max = 2/(om0*(1 + 2*sqrt(r) - 2*(r^.25)));
% assert(k <= k_max, "Error, unstable for the Forward Euler Method")

% initialize
I = eye(4) ;
A = om0*[-1 0 0 -4*r; 1 -1 0 0; 0 1 -1 0; 0 0 1 -1];
b = om0*[1 0 0 0]';
c = [0 0 0 1]';

Bt_1 = I + k*A/2; % 4 by_NL4 matrix used for matrix multiplication in TI
Bt_2 = I - k*A/2; % 4 by_NL4 matrix used for matrix inversion or linear system solution in TI

x_linear_imp = zeros(4, 1);    % initialize 4 by_NL1 vector for linear state of impulse input
x_linear_sine = zeros(4, 1);   % initialize 4 by_NL1 vector for linear state of sinusoidal input
x_NL_Imp = zeros(4, 1);        % initialize 4 by_NL1 vector for non linear state of impulse input
x_NL_Sine = zeros(4, 1);       % initialize 4 by_NL1 vector for non linear state of sinusoidal input

y_linear_imp = zeros(Nf, 1);      % initialize the output vector for linear reponse of impulse input
y_linear_sine = zeros(Nf, 1);     % initialize the output vector for linear reponse of sinusoidal input
y_NL_Imp= zeros(Nf, 1);           % initialize the output vector for non linear response of impulse input
y_NL_Sine= zeros(Nf, 1);          % initialize the output vector for non linear response of sinusoidal input

gain = 5;                               % input gain
uImp = gain*[1; zeros(Nf-1, 1)];        % initialize the impulse input vector
uSine = gain*sin(2*pi*f0*(1:Nf)/SR);    % initialize the sinusoidal input vectpr

u0_Imp = 0;                             % intial sample of the impulse input
u0_Sine = 0;                            % initial sample of the sinusoidal input
fvec = (0:Nf-1)'*SR/Nf;                 % frequency vector for plots

Hc = zeros(Nf,1);                       % continuous-time transfer function 
omega = (0:Nf-1)'*2*pi*SR/Nf;           % angular frequency_NLvector
%++++++++++++++++++++++++++++++++++++++++
% main loop
tic;
for n = 1 : Nf

    % Linear response
    %
    % update state x_linearfrom n to n+1, with current sample of u as input (TI)
    x_linear_imp= Bt_2 \ (Bt_1*x_linear_imp + k*b*(u0_Imp + uImp(n))/2);
    x_linear_sine = Bt_2 \ (Bt_1*x_linear_sine+ k*b*(u0_Sine + uSine(n))/2);

    %--------------------------------------------------------------------
    % Non-linear response
    %
    % Since the non linear function f(x, u) now involves the current state 
    % x, so we have to solve for the root of the trapezoid difference equation.
    % One common method is to use the root finding algorithm Newton-Raphson
    % method.
    x_NL_Imp_next= NewtonRaphsonMOOG(x_NL_Imp, uImp(n), u0_Imp, k, om0, r, 100, 1e-8);
    x_NL_Sine_next = NewtonRaphsonMOOG(x_NL_Sine, uSine(n), u0_Sine, k, om0, r, 100, 1e-8);

    % update the input sample
    u0_Imp = uImp(n);
    u0_Sine = uSine(n);
    % update the state for the non-linear response
    x_NL_Imp = x_NL_Imp_next;
    x_NL_Sine = x_NL_Sine_next;
    
    % write sample to output vector y_linear_imp
    y_linear_imp(n) = c' * x_linear_imp;
    % write sample to output vector y_linear_sine
    y_linear_sine(n) = c' * x_linear_sine;
    % write ssample to output vector y_NL_Imp
    y_NL_Imp(n) = c' * x_NL_Imp_next;
     % write ssample to output vector y_NL_Sine
    y_NL_Sine(n) = c' * x_NL_Sine_next;
    % write sampel to continuous-time transfer funciton (CT)
    Hc(n) = c'*((1i*omega(n)*eye(4) - A)\b);
    
end
toc ;
%++++++++++++++++++++++++++++++++++++++++

% discrete transfer functions as computed from IR
H_linear = fft(y_linear_imp) ; % linear FR
H_NL = fft(y_NL_Imp)         ; % non-linear FR

% sinusoidal response
Y_linear_sine = fft(y_linear_sine);
Y_NL_Sine = fft(y_NL_Sine); 
%% 

% plot magnitude response 

figure

% plot the magnitude response of an impulse input

subplot(2,1,1)
semilogx(fvec, 20*log10(abs(H_linear)))   % linear response
hold on
semilogx(fvec, 20*log10(abs(H_NL)))       % non-linear response
hold on 
semilogx(fvec, 20*log10(abs(Hc)))         % CT response (linear)

title(sprintf('log plot of impulse input, r = %.2f, f0 = %d', r, f0))
xlabel('log frequency(Hz)')
ylabel('magnitude dB')
legend('Linear', 'NonLinear', 'ContinuousTime')
xlim([0, SR/2])
ylim([-200, 30])
grid on


% plot the magnitude response of a sinusoidal input

subplot(2,1,2)
semilogx(fvec, 20*log10(abs(Y_linear_sine)))  % linear response
hold on
semilogx(fvec, 20*log10(abs(Y_NL_Sine)))      % non-linear resposne


title(sprintf('log plot of sinusoidal input, r = %.2f, f0 = %d', r, f0))
xlabel('log frequency(Hz)')
ylabel('magnitude dB')
legend('Linear', 'NonLinear')
xlim([0, SR/2])
ylim([-50, 120])
grid on

%% 

function x_next = NewtonRaphsonMOOG(x, u, u0, k, om0, r, max_iter, eps)
    % initialize guess for x
    x_next= x;

    for n = 1:max_iter
        % Compute f_next(x_next, u)
        f_next= om0*[-tanh(x_next(1)) + tanh(u - 4*r*x_next(4));
                 -tanh(x_next(2)) + tanh(x_next(1));
                 -tanh(x_next(3)) + tanh(x_next(2));
                 -tanh(x_next(4)) + tanh(x_next(3))];
        % Compute f(x, u0)
        f = om0*[-tanh(x(1)) + tanh(u0 - 4*r*x(4));
                 -tanh(x(2)) + tanh(x(1));
                 -tanh(x(3)) + tanh(x(2));
                 -tanh(x(4)) + tanh(x(3))];

        % Compute G_x_next
        G_x = x_next- x - (k/2)*(f_next + f);

        % Compute the jacobian matrix of f, J_f(f_next, u)
        J_f = om0 * diag([
            -sech(x_next(1))^2, ...
            -sech(x_next(2))^2, ...
            -sech(x_next(3))^2, ...
            -sech(x_next(4))^2
        ]); 
        
        % Adjust off-diagonal terms
        J_f(2,1) = om0 * sech(x_next(1))^2;
        J_f(3,2) = om0 * sech(x_next(2))^2;
        J_f(4,3) = om0 * sech(x_next(3))^2;
        J_f(1,4) = om0 *  -4*r * sech(u - 4*r*x_next(4))^2;

        % Compute Jacobian of G(x)
        J_G = eye(4) - (k/2) * J_f;

        % updata x_next
        delta_x = J_G \ G_x;
        x_next= x_next - delta_x;

        % check for convergence
        % if norm(delta_x) < eps
        if norm(x_next - x) < eps
            break;
        end
    end

end