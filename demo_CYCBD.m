%% TEST OF MAXIMUM (SECOND-ORDER) CYCLOSTATIONARITY BLIND DECONVOLUTION
%  This example illustrates the detection of a cyclostationary source in
%  various harsh conditions
%
%------------------------------------------------------------------
% Code by and M. Buzzoni,  J. Antoni and G. D'Elia \\ December 2017
%------------------------------------------------------------------

clear
close all
clc
disp('-----------------------------------------------------------------------------')
disp('This example illustrates the detection of a cyclostationary train of impulses')
disp('in various harsh conditions')
disp('-----------------------------------------------------------------------------')

disp('case 1: train of impulses having Gaussian distributed amplitudes and additive')
disp('        Gaussian background noise')
disp('case 2: train of impulses having Gaussian distributed amplitudes, jitter')
disp('        effect and additive Gaussian background noise')
disp('case 3: a couple of impulse trains with different cyclic frequencies having')
disp('        Gaussian distributed amplitude and additive Gaussian background noise')
disp('case 4: train of impulses with Gaussian distributed amplitudes and additive')
disp('        Gaussian background noise with the addition of a single dominant impulse')
disp('case 5: train of impulses with Gaussian distributed amplitudes having')
disp('        fluctuating cycle and additive Gaussian background noise')
disp('case 6: train of impulses with Gaussian distributed amplitudes convoluted with two')
disp('        different IRFs (SIMO case) and additive Gaussian background noise')
fprintf('\n')
nExample = input('Select the case number: ');

if isempty(find(nExample == 1:6, 1))
    disp('Script aborted: you must select a case number between 1 and 6')
    return
end

disp(['Running case ' num2str(nExample) '...'])

switch nExample
    case 1
        load('test_signals/signal_case1.mat')
        
        % CYCBD
        N = 40;
        param.RE = 1e-4;
        param.iter = 100;
        [h_est,s0,kappa,W,count,err] = MaxCycloBD(y,N,alpha,fs,param,2);
        
        disp(['Number of iterations = ',num2str(count)])
        disp(['Estimated ICS2 = ',num2str(kappa')])
        
        % plots
        figure
        subplot(4,1,1)
        plot(s,'k')
        title('source signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,2)
        plot(x,'k')
        title('response signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,3)
        plot(y,'k')
        title('observed noisy signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,4)
        plot(s0,'k')
        title('estimated source'), xlabel('samples'), box off, set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(2,1,1)
        SE = abs(hilbert(abs(y).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        xline = alpha*length(x)/fs; xline = xline(1:20);
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (observed signal)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        subplot(2,1,2)
        SE = abs(hilbert(abs(s0).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (estimated source)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        
        pause(1)
        
        figure
        loglog(err,'k'),
        title('relative error on ICS_2'), xlabel('iterations'), set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(211),plot(h_est,'k'),xlabel('samples'),title('inverse impulse response'), box off, set(gca,'fontsize',8)
        subplot(212),plot(10*log10(abs(fft(h_est))),'k'), xlim([0 N/2])
        xlabel('non-dimensional frequency'), ylabel('dB'),title('inverse transfer function'), box off, set(gca,'fontsize',8)
        
        pause(1)
        
        figure, plot(W,'k')
        title('weights at convergence')
        xlabel('samples'), set(gca,'fontsize',8), box off
        
    case 2
        
        load('test_signals/signal_case2.mat')
        
        % CYCBD
        N = 40;
        param.RE = 1e-4;
        param.iter = 100;
        [h_est,s0,kappa,W,count,err] = MaxCycloBD(y,N,alpha,fs,param,2);
        
        disp(['Number of iterations = ',num2str(count)])
        disp(['Estimated ICS2 = ',num2str(kappa')])
        
        % plots
        figure
        subplot(4,1,1)
        plot(s,'k')
        title('source signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,2)
        plot(x,'k')
        title('response signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,3)
        plot(y,'k')
        title('observed noisy signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,4)
        plot(s0,'k')
        title('estimated source'), xlabel('samples'), box off, set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(2,1,1)
        SE = abs(hilbert(abs(y).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        xline = alpha*length(x)/fs; xline = xline(1:20);
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (observed signal)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        subplot(2,1,2)
        SE = abs(hilbert(abs(s0).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (estimated source)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        
        pause(1)
        
        figure
        loglog(err,'k'),
        title('relative error on ICS_2'), xlabel('iterations'), set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(211),plot(h_est,'k'),xlabel('samples'),title('inverse impulse response'), set(gca,'fontsize',8), box off
        subplot(212),plot(10*log10(abs(fft(h_est))),'k'), xlim([0 N/2])
        xlabel('non-dimensional frequency'), ylabel('dB'), title('inverse transfer function'), set(gca,'fontsize',8), box off
        
        pause(1)
        
        figure, plot(W,'k')
        title('weights at convergence')
        xlabel('samples'), set(gca,'fontsize',8), box off
        
    case 3
        
        load('test_signals/signal_case3.mat')
        
        % CYCBD (source 1)
        N = 40;
        param.RE = 1e-4;
        param.iter = 100;
        [h_est,s0,kappa,W,count,err] = MaxCycloBD(y,N,alpha,fs,param,2);
        
        disp(['Number of iterations (source 1) = ',num2str(count)])
        disp(['Estimated ICS2 (source 1) = ',num2str(kappa')])
        
        % plots (source 1)
        figure
        subplot(4,2,1)
        plot(s,'k')
        title('source signal 1'), xlabel('samples'), box off, set(gca,'fontsize',8),
        subplot(4,2,2)
        plot(s_second,'k')
        title('source signal 2'), xlabel('samples'), box off, set(gca,'fontsize',8),
        subplot(4,2,[3 4])
        plot(x + x_second,'k')
        title('response signal'), xlabel('samples'), box off, set(gca,'fontsize',8),
        subplot(4,2,[5 6])
        plot(y,'k')
        title('noisy observed signal'), xlabel('samples'), box off, set(gca,'fontsize',8),
        subplot(4,2,[7 8])
        plot(s0,'k')
        title('estimated source 1'), xlabel('samples'), box off, set(gca,'fontsize',8),
        
        figure
        subplot(2,1,1)
        SE = abs(hilbert(abs(y).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        xline = alpha*length(x)/fs; xline = xline(1:20);
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (observed signal)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        subplot(2,1,2)
        SE = abs(hilbert(abs(s0).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (estimated source 1)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        
        pause(1)
        
        figure
        loglog(err,'k')
        title('relative error on ICS_2'), xlabel('iterations'), set(gca,'fontsize',8), box off
        
        pause(1)
        
        figure
        subplot(211),plot(h_est,'k'), xlabel('samples'), title('inverse impulse response (source 1)'), set(gca,'fontsize',8)
        subplot(212),plot(10*log10(abs(fft(h_est))),'k'),xlim([0 N/2]), box off
        xlabel('Normalized frequency'),ylabel('dB'),title('inverse transfer function (source 1)'), set(gca,'fontsize',8), box off
        
        pause(1)
        
        figure, plot(W,'k')
        title('weights at convergence (source 1)'), xlabel('samples'), box off, set(gca,'fontsize',8)
        
        pause(1)
        
        % CYCBD (source 2)
        N = 40;
        param.RE = 1e-4;
        param.iter = 100;
        [h_est,s0,kappa,W,count,err] = MaxCycloBD(y,N,alpha_second,fs,param,2);
        
        disp(['Number of iterations (source 2) = ',num2str(count)])
        disp(['Estimated ICS2 (source 2) = ',num2str(kappa')])
        
        % plots (source 2)
        figure
        subplot(4,2,1)
        plot(s,'k')
        title('source signal 1'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,2)
        plot(s_second,'k')
        title('source signal 2'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,[3 4])
        plot(x + x_second,'k')
        title('response signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,[5 6])
        plot(y,'k')
        title('noisy observed signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,[7 8])
        plot(s0,'k')
        title('estimated source (source 2)'), xlabel('samples'), box off, set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(2,1,1)
        SE = abs(hilbert(abs(y).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        xline = alpha_second*length(x)/fs; xline = xline(1:40);
        l1 = line([xline(:) xline(:)]',[zeros(40,1) ones(40,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (observed signal)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        subplot(2,1,2)
        SE = abs(hilbert(abs(s0).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        l1 = line([xline(:) xline(:)]',[zeros(40,1) ones(40,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (estimated source 2)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        
        pause(1)
        
        figure
        loglog(err,'k')
        title('relative error on ICS_2 (source 2)'), xlabel('iterations'), set(gca,'fontsize',8), box off
        
        pause(1)
        
        figure
        subplot(211),plot(h_est,'k'),xlabel('samples'),title('inverse impulse response (source 2)'), set(gca,'fontsize',8)
        subplot(212),plot(10*log10(abs(fft(h_est))),'k'),xlim([0 N/2])
        xlabel('Normalized frequency'),ylabel('dB'),title('inverse transfer function (source 2)'), set(gca,'fontsize',8)
        
        pause(1)
        
        figure, plot(W,'k')
        title('weights at convergence (source 2)'), xlabel('samples'), box off, set(gca,'fontsize',8), box off
        
    case 4
        
        load('test_signals/signal_case4.mat')
        
        % CYCBD
        N = 40;
        param.RE = 1e-4;
        param.iter = 100;
        [h_est,s0,kappa,W,count,err] = MaxCycloBD(y,N,alpha,fs,param,2);
        
        disp(['Number of iterations = ',num2str(count)])
        disp(['Estimated ICS2 = ',num2str(kappa')])
        
        % plots
        figure
        subplot(4,2,1)
        plot(s,'k')
        title('source signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,2)
        plot(s_second,'k')
        title('dominant single impulse'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,[3 4])
        plot(x + x_second,'k')
        title('response signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,[5 6])
        plot(y,'k')
        title('noisy observed signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,[7 8])
        plot(s0,'k')
        title('estimated '), xlabel('samples'), box off, set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(2,1,1)
        SE = abs(hilbert(abs(y).^2)); SES = abs(fft(SE)).^2; SES(1) = 0;
        xline = alpha*length(x)/fs; xline = xline(1:20);
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (observed signal)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        subplot(2,1,2)
        SE = abs(hilbert(abs(s0).^2)); SES = abs(fft(SE)).^2; SES(1) = 0;
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (estimated source)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        
        pause(1)
        
        figure
        loglog(err,'k'), title('relative error on kappa'), xlabel('iterations'), set(gca,'fontsize',8), box off
        
        pause(1)
        
        figure
        subplot(211),plot(h_est,'k'),xlabel('samples'),title('Inverse impulse response')
        set(gca,'fontsize',8), box off
        subplot(212),plot(10*log10(abs(fft(h_est))),'k'),xlim([0 N/2])
        xlabel('Normalized frequency'),ylabel('dB'),title('Inverse transfer function'), set(gca,'fontsize',8), box off
        
        pause(1)
        
        figure, plot(W,'k')
        title('weights at convergence (source 1)'), xlabel('samples'), box off, set(gca,'fontsize',8)
        
    case 5
        
        load('test_signals/signal_case5.mat')

        % CYCBDang
        N = 250;
        param.RE = 1e-4;
        param.iter = 100;
        [h_est,s0,kappa,W,count,err] = MaxCycloBDangle(y,N,alpha,fs,tachoPulse,pulseXrev,param,2);
        
        disp(['Number of iterations = ',num2str(count)])
        disp(['Estimated ICS2 = ',num2str(kappa')])
        
        % plots
        figure
        subplot(4,1,1)
        plot(s,'k')
        title('source signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,2)
        plot(x,'k')
        title('response signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,3)
        plot(y,'k')
        title('noisy observed signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,1,4)
        plot(s0,'k')
        title('estimated '), xlabel('samples'), box off, set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(2,1,1)
        [angle_y] = TimeToAngle(y,fs,t_trig,20,z);
        [angle_s0] = TimeToAngle(s0,fs,t_trig,20,z); angle_s0 = angle_s0(1:end-20*z);
        SE = abs(hilbert(abs(angle_y).^2)); SES = abs(1/(length(angle_y)).*fft(SE)).^2; SES(1) = 0;
        xline = floor((1:8)*length(angle_y)/(20*z));
        l1 = line([xline(:) xline(:)]',[zeros(8,1) ones(8,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (observed signal)'), xlabel('non-dimensial order'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse order')
        subplot(2,1,2)
        SE = abs(hilbert(abs(angle_s0).^2)); SES = abs(fft(SE)).^2; SES(1) = 0;
        l1 = line([xline(:) xline(:)]',[zeros(8,1) ones(8,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(0:length(SES)-1,SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (estimated source)'), xlabel('non-dimensional order'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse order')
        
        pause(1)
        
        figure
        subplot(211),plot(h_est,'k'),xlabel('samples'),title('inverse impulse response')
        set(gca,'fontsize',8), box off
        subplot(212),plot(10*log10(abs(fft(h_est))),'k'),xlim([0 N/2])
        xlabel('Normalized frequency'),ylabel('dB'),title('inverse transfer function'), set(gca,'fontsize',8), box off
        
        pause(1)
        
        figure, plot(W,'k')
        title('weights at convergence'), xlabel('samples'), box off, set(gca,'fontsize',8)
        
    case 6
        
        load('test_signals/signal_case6.mat')
        load signaltest
        % CYCBD SIMO
        N = 40;
        param.RE = 1e-4;
        param.iter = 100;
        [h_est,s0,kappa,W,count,err] = MaxCycloBD_SIMO(Y,N,alpha,fs,param,2);
        
        disp(['Number of iterations = ',num2str(count)])
        disp(['Estimated ICS2 = ',num2str(kappa')])
        
         % plots
        figure
        subplot(4,2,[1 2])
        plot(s,'k')
        title('source signal'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,3)
        plot(x,'k')
        title('response signal 1'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,4)
        plot(x_second,'k')
        title('response signal 2'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,5)
        plot(Y(:,1),'k')
        title('observed noisy signal 1'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,6)
        plot(Y(:,2),'k')
        title('observed noisy signal 2'), xlabel('samples'), box off, set(gca,'fontsize',8)
        subplot(4,2,[7 8])
        plot(s0,'k')
        title('estimated source'), xlabel('samples'), box off, set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(2,1,1)
        SE = abs(hilbert(abs(y).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        xline = alpha*length(x)/fs; xline = xline(1:20)+1;
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (observed signal)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        subplot(2,1,2)
        SE = abs(hilbert(abs(s0).^2)); SES = abs(1/(length(SE)).*fft(SE)).^2; SES(1) = 0;
        l1 = line([xline(:) xline(:)]',[zeros(20,1) ones(20,1).*max(SES)*1.2]','Color','k','Linestyle',':'); hold on
        plot(SES,'k'), xlim([0 200]), ylim([0 max(SES)*1.2])
        title('SES (estimated source)'), xlabel('non-dimensional frequency'), set(gca,'fontsize',8), box off
        legend(l1(1,:),'impulse frequency')
        
        pause(1)
        
        figure
        loglog(err,'k'),
        title('relative error on ICS_2'), xlabel('iterations'), set(gca,'fontsize',8)
        
        pause(1)
        
        figure
        subplot(211),plot(h_est(:,1),'k'),hold on, plot(h_est(:,2),'color',[.5 .5 .5])
        xlabel('samples'),title('inverse impulse response'), box off, set(gca,'fontsize',8)
        legend('inv. IR 1','inv. IR 2')
        subplot(212),plot(10*log10(abs(fft(h_est(:,1)))),'k'), hold on, plot(10*log10(abs(fft(h_est(:,2)))),'color',[.5 .5 .5]), xlim([0 N/2])
        xlabel('non-dimensional frequency'), ylabel('dB'),title('inverse transfer function'), box off, set(gca,'fontsize',8)
        legend('inv. TF 1','inv. TF 2')
        
        pause(1)
        
        figure, plot(W,'k')
        title('weights at convergence')
        xlabel('samples'), set(gca,'fontsize',8), box off
        
        pause(1)
        
        % comparison between given filter and estimated filter
        Nw = 2^9;
        Nfft = 2*Nw;
        Hdir = zeros(Nfft/2+1,3);
        for k = 1:2
            Hdir(:,k) = tfestimate(s0,Y((N:length(s)),k),Nw,fix(2/3*Nw),Nfft);
        end
        
        figure
        subplot(211),plot(10*log10(abs(Hdir(:,1))),'k'), hold on, plot(10*log10(abs(Hdir(:,2))),'color',[.5 .5 .5]), xlim([0 length(Hdir)])
        title('direct filters through estimated source'), xlabel('non-dimensional frequency'), ylabel('dB'), set(gca,'fontsize',8)
        h = [himp; himp_second]';
        Hest = abs(fft(h));
        subplot(212), plot(linspace(0,2*length(Hdir),length(Hest)),10*log10(Hest(:,1)),'k'), hold on
        plot(linspace(0,2*length(Hdir),length(Hest)),10*log10(Hest(:,2)),'k'),
        xlim([0 length(Hdir)])
        title('given filters'), xlabel('non-dimensional frequency'), ylabel('dB'), ylim([-60 -10])
        set(gca,'fontsize',8)
end

%% % Angular resampling of a time signal
% References:
% [1]   K. R. Fyfe and E. D. S. Munck, “Analysis of Computed Order Tracking,” Mechanical Systems
%       and Signal Processing, vol. 11, no. 2, pp. 187–205, 1997.
%
% Code by G. D'Elia and M. Buzzoni

function [AngleSignal,indexNaN] = TimeToAngle(TimeSignal,fs,t_trig,Mpulse,pulsePerRev)
t = (0:length(TimeSignal)-1)/fs;

k = Mpulse/2:3*Mpulse/2-1;
deltatheta = 2*pi/(pulsePerRev*Mpulse);
theta = k*deltatheta;
index = 1:length(k);
tk = zeros(1,length(k)*(length(t_trig) - 2));
for ind = 1:length(t_trig) - 2;
    tMatrix = [1 t_trig(ind) t_trig(ind)^2
        1 t_trig(ind+1) t_trig(ind+1)^2
        1 t_trig(ind+2) t_trig(ind+2)^2];
    DeltaPhi = [0
        2*pi/pulsePerRev
        2*2*pi/pulsePerRev];
    b = pinv(tMatrix)*DeltaPhi;
    tk(1,index) = 1/(2*b(3)) .* (sqrt(b(2)^2 + 4*b(3)*(theta-b(1))) - b(2));
    index = index + length(k);
end
indexNaN = find(isnan(tk));

AngleSignal = interp1(t,TimeSignal,tk,'spline');
end


