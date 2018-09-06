function [h,s,kappa,W,count,err] = MaxCycloBD_SIMO(x,N,alpha,fs,param,p)
% [h,s,kappa,W,count,err,extr] = MaxCYCBD_SIMO(x,N,alpha,fs,param,h,p)
% SIMO BLIND DECONVOLUTION VIA MAXIMUM CYCLOSTATIONARITY
%  Simultaneous blind deconvolution of multiple signals x (single-input-multiple-output model)
%  by finding the optimal inverse filters (FIR form) that maximizes the cyclostationarity of the output.
%
%--------
% Inputs
%--------
%
% x.........observed signals of length L organized as a matrix LxK where K
%           is the number of observed signals
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
% h_final...matrix of N-long optimal inverse filters stacked in K columns
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

[L,K] = size(x);

for k = 1:K
    x(:,k) = x(:,k) - mean(x(:,k));
end

if nargin < 5
    param.RE = 1e-4;
    param.iter = 5;
end
% inverse FIR filters initialization
h = zeros(N*K,1);
for k = 1:K
    h((k-1)*N+1:k*N) = lpc(x(:,k),N-1).';
end

if nargin < 6 
    p = 2;
end

XX = CorrMatrix_SIMO(x,[],N);

test = 0;
count = 1;
kappa_old = 0;
err = zeros(param.iter,1);
while test == 0
    s = 0;
    for k = 1:K
        s = s + filter(h((k-1)*N+1:k*N),1,x(:,k));
    end
    W = abs(s(N:L)).^p;
    W = Periodic(W,alpha,fs);
    W = W/(mean(W).^(p/2));
    
    XWX = CorrMatrix_SIMO(x,W,N);
    [h,kappa] = eigs(XWX,XX,1);
    kappa = diag(kappa);
%     if param.n > 1
%         [kappa,I] = sort(kappa,'descend');
%         h = h(:,I);
%     end
    err(count) = abs(kappa(1)-kappa_old)/abs(kappa_old);
    if (err(count) < param.RE) || count >= param.iter
        test = 1;
    end
    count = count + 1;
    kappa_old = kappa(1);

end
count = count - 1;
h = reshape(h,N,K);
s = s(N:L,:);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = CorrMatrix_SIMO(x,W,N)
% Compute augmented correlation matrix (NKxNK) with weights W of K signals
% arranged in the colmuns of matrix x

K = size(x,2);

R = zeros(N*K,N*K);
for i = 1:K
    R((i-1)*N+1:i*N,(i-1)*N+1:i*N) = CorrMatrix(x(:,i),W,N);
    for j = i+1:K
        R((i-1)*N+1:i*N,(j-1)*N+1:j*N) = XCorrMatrix(x(:,i),x(:,j),W,N);
        R((j-1)*N+1:j*N,(i-1)*N+1:i*N) = R((i-1)*N+1:i*N,(j-1)*N+1:j*N)';
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = XCorrMatrix(x,y,W,N)
% Compute cross-correlation matrix (NxN) of signals x and y with weights W

L = length(x);
R = zeros(N);

if isempty(W)
    W = ones(L-N+1,1);
end

W = W(:);
x = x(:);
for i = 1:N
    for j = 1:N
        R(i,j) = mean(x(N+1-i:L+1-i).*conj(y(N+1-j:L+1-j)).*W);
    end
end
end

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = Periodic(W,alpha,fs)
% Extract periodic component with frequencies alpha in x

W = W(:);
L = length(W);
dt = 1/fs;
T = dt*L;
t = (0:dt:T-dt)';
K = length(alpha);

alpha(alpha==0) = [];
p = zeros(length(W),1);
p = mean(W);
for k = 1:K
    c = mean(W.*exp(-2i*pi*alpha(k).*t));
    p = p + 2*real(c*exp(2i*pi*alpha(k).*t));
end

% hard thresholding for improving results,
th = quantile(p,3); th = th(3);
p(p < th) = 0;

end