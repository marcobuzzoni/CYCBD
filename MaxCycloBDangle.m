function [h,s,kappa,W,count,err] = MaxCycloBDangle(x,N,alpha,fs,tPulse,NPulse,param,p)
% [h,s,kappa,count,err] = MaxCycloBDangle(x,N,alpha,fs,param,p)
% SISO BLIND DECONVOLUTION VIA MAXIMUM CYCLOSTATIONARITY (ANGLE/TIME DOMAIN
% VERSION)
%  Blind deconvolution of signal x by finding the optimal inverse filter (FIR form) 
%  that maximizes the cyclostationarity of the output through an weighting matrix
%  defined in the time/angle domain.
%
%--------
% Inputs
%--------
%
% x.........observed signal
% N.........FIR filter length
% alpha.....cyclic order set (in the form of a vector)
% fs........sampling frequency of x
% tPulse....time occurrences of the angular reference
% Npulse....number of occurrences per revolution
% param.....structure of setting parameters organized as follows:
%                param.ER......minimal relative error on result (default value = 1e-3)
%                param.iter....maximum number of iterations (default value  = 50)
% p.........cyclostationarity order to maximize (default = 2)
%
%---------
% Outputs
%---------
% 
% s.........blindly deconvolved signal
% h_final...optimal inverse FIR filter at convergence
% kappa.....value of criterion at convergence
% W.........weights used in the criterion at convergence
% count.....number of iteration to convergence
% err.......relative error on result as a function of iterations
%
%-----------
% Reference
%-----------
%  M. Buzzoni, J. Antoni, G. D'Elia, "Blind deconvolution based on
%  cyclostationarity maximization and its application to fault
%  identification", Journal of Sound and Vibration, 2018, Accepted 
%
%-------------------------------------------------
% Code by J. Antoni and M. Buzzoni, Dicember 2017
%-------------------------------------------------

if nargin < 2
    disp('Matlab function aborted: FIR filter length is missing.')
    return
end

if nargin < 3
    disp('Matlab function aborted: cyclic frequency set is missing.')
    return
end

if nargin < 4
    disp('Matlab function aborted: sampling frequency is missing.')
    return
end

if nargin < 5
    disp('Matlab function aborted: tacho signal is missing.')
    return
end

if nargin < 6
    disp('Matlab function aborted: tacho pulses per revolution is missing')
    return
end

L = length(x);
x = x - mean(x);

if nargin < 7
    param.RE = 1e-3;
    param.iter = 50;
end

% inverse FIR filters initialization
h = lpc(x,N-1);
h = h(:);

XX = CorrMatrix(x,[],N);

test = 0;
count = 1;
kappa_old = 0;
err = zeros(param.iter,1);

while test == 0
    s = filter(h,1,x);
    W = abs(s(N:L)).^p;
    W = periodicAngle(W,tPulse,fs,alpha,NPulse);
    W = W/(mean(W).^(p/2));
    XWX = CorrMatrix(x,W,N);
    [h,kappa] = eigs(XWX,XX,1); 
    err(count) = abs(kappa-kappa_old)/abs(kappa_old);
    if ((err(count) < param.RE) || count >= param.iter)
        test = 1;
    end
    count = count + 1;
    kappa_old = kappa;
end
count = count - 1;
s = s(N+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = CorrMatrix(x,W,N)
% Compute correlation matrix (NxN) of signal x with weights W

L = length(x);
R = zeros(N);

if isempty(W)
    W = ones(L-N+1,1);
end

W = W(:);
x = x(:);
for i = 1:N
    R(i,i) = mean(abs(x(N+1-i:L+1-i)).^2.*W);
    for j = i+1:N
        R(i,j) = mean(x(N+1-i:L+1-i).*conj(x(N+1-j:L+1-j)).*W);
        R(j,i) = conj(R(i,j));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extracton of time components periodic in the angle domain
%
% M. Buzzoni
% 09/03/2017
%
%   y = periodicAngle(y,tPulse,Fs,O,Npulse)
%
% y = input signal
% tPulse = time of the tacho impulses
% Fs = sampling frequency
% O = orders to extract
% Npulse = number of tacho impulses per revolution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = periodicAngle(x,tPulse,fs,O,NPulse)

L = length(x);
t = (0:L-1)./fs; % time
T = t(end);
tPulse = tPulse(tPulse < T(end));

fr = 1./diff(tPulse)/NPulse; % instantaneous frequency
fr = interp1(tPulse(2:end),fr,t,'spline');
thetaDot = fr.*2.*pi;
theta = cumtrapz(thetaDot)./fs; % angles
Theta = theta(end);
% figure
if length(theta) == length(x) && length(thetaDot) == length(x)
    y = mean(x);
    for k = 1:length(O)
        ck = 2/(Theta*fs).*sum(x'.*exp(-1i.*O(k).*theta).*thetaDot);
        yk = ck.*exp(1i.*k.*theta);
        yk = real(yk);
%         subplot(211), plot(yk)
        y = y + yk;
%         subplot(212), plot(y)
%         pause(.25)
    end
else
    disp('error')
end

% hard thresholding for improving results,
th = quantile(y,3); th = th(3);
y(y < th) = 0;

