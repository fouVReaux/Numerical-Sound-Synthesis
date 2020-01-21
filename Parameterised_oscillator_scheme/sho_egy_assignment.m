% oscillator assignment
% S. Bilbao, 14 December 2019

% In this assignment, a code is provided for the parameterised oscillator scheme,
% as given in Eq. (58) in the notes. The scheme depends on an angular frequency w0, as well as
% a free paramater alpha.
% You should derive a scheme of the form of:
% u^{n+1} = a*u^{n}-u^{n-1}, for some constant a, which you will determine
% below.

% There are three questions Q1 through Q3 below. Please edit this Matlab
% code and return it with the questions answered!

clear all
close all

% parameters

Fs = 44100;
f0 = 1842;
Tf = 1;
u0 = 1;
v0 = 0;
alpha = 0.8;

% derived quantities

k = 1/Fs;
w0 = 2*pi*f0;
Nf = floor(Tf*Fs);

% Q1: Insert proper value of the parameter a to be used in the update in the
% main loop (see comments at top of assignment).

a = (2/(k^2) - (w0^2)*(alpha + 1)/2)/(1/(k^2) + (w0^2)*(1-alpha)/4)

% Q2: check stability condition, which gives bound on k, given alpha and w0 (see the notes). Here, perform a check to make sure that the stability condition is
% satisfied, for any values of the parameters. If it is not, the code
% should exit with an error message.

if k <= 2/w0*sqrt(alpha)
    fprintf('INFO : Stability condition : satisfied!\n');
else
    error('WARNING! : Stability condition : NOT satisfied!\n');
    return;
end

% initialize

u2 = u0;
u1 = u0+k*v0;
H = zeros(Nf,1);
out = zeros(Nf,1);

% main loop

tic

for n=1:Nf
    u = a*u1-u2;
    out(n) = u2;

    % Q3: derive an expression for conserved energy and insert here below:
    H(n) = 0.5*(u1-u2)^2/k^2 + 0.5*w0^2*(alpha*u1*u2 + (1 - alpha)*(0.5*(u1+u2)^2));

    u2 = u1;
    u1 = u;
end

toc

% plot

tax = [0:Nf-1]'*k;

Herr = (H-H(1))/H(1);
plot(tax,Herr,'k.');
xlabel('t');
ylabel('H');
title('Energy variation');

% play sound

soundsc(out,Fs)
