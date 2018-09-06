function [h_final,s,kappa,W,count,err] = MaxCycloBD(x,N,alpha,fs,param,p)
% [h,s,kappa,count,err] = MaxCycloBD(x,N,alpha,fs,param,p)
% SISO BLIND DECONVOLUTION VIA MAXIMUM CYCLOSTATIONARITY
%  Blind deconvolution of signal x by finding the optimal inverse filter (FIR form) 
%  that maximizes the cyclostationarity of the output.
%
%--------
% Inputs
%--------
%
% x.........observed signal
% N.........FIR filter length
% alpha.....cyclic frequency set (in the form of a vector)
% fs........sampling frequency of x
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
    disp('Matlab function aborted:FIR filter length is missing.')
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

L = length(x);
x = x - mean(x);

if nargin < 5
    param.RE = 1e-3;
    param.iter = 50;
end

% inverse FIR filters initialization
h = lpc(x,N-1);
h = h(:);

if nargin < 6 
    p = 2;
end

XX = CorrMatrix(x,[],N);

test = 0;
count = 1;
kappa_old = 0;
err = zeros(param.iter,1);

while test == 0
 s = filter(h,1,x);
    W = abs(s(N:L)).^p;
    W = Periodic(W,alpha,fs);
    W = W/(mean(W).^(p/2));
    XWX = CorrMatrix(x,W,N);
    [h,kappa] = eigs(XWX,XX,1);
    err(count) = abs(kappa-kappa_old)/abs(kappa_old);
    if (err(count) < param.RE) || count >= param.iter 
        test = 1;
    end
    count = count + 1;
    kappa_old = kappa;    
%     s = filter(h,1,x);
%     W = abs(s(N:L)).^2;
%     W = Periodic(W,alpha,fs);
%     W = W/mean(W);
%     XWX = CorrMatrix(x,W,N);
%     [h,kappa] = eigs(XWX,XX,1);
%     err(count) = abs(kappa-kappa_old)/abs(kappa_old);
%     if (err(count) < param.RE) || count >= param.iter 
%         test = 1;
%     end
%     count = count + 1;
%     kappa_old = kappa;
end
h_final = h;
count = count - 1;
s = filter(h,1,x);
s = s(N:end);
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
function p = Periodic(x,alpha,fs)
% Extract cyclic components with frequencies alpha in x

x = x(:);
L = length(x);
dt = 1/fs;
T = dt*L;
t = (0:dt:T-dt)';
K = length(alpha);

alpha(alpha==0) = [];
p = mean(x);
for k = 1:K
    c = mean(x.*exp(-2i*pi*alpha(k).*t));
    p = p + 2*real(c*exp(2i*pi*alpha(k).*t));
end

% % thresholding for improve weighting matrix
th = mean(p)+2*std(p);
% th = quantile(p,3)
% th = th(3)
p(p<th) = 0;
