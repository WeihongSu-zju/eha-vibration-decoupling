%% 1. Frequency extraction simulation
fs = 1e4;  T = 5;
t  = (1/fs:1/fs:T)';   Nt = numel(t);

baseF = 500;
mults = 1:5;   K = numel(mults);

% Simulation signal generation
sinTerm = 100*sin(pi*t);
fk   = t*(baseF*mults/T) + sinTerm;
phi  = 2*pi*cumtrapz(t, fk);
Yamp = t*(mults/T) + 0.1*sin(pi*t);
Y    = sum(Yamp .* cos(phi), 2);

% Adding noise to the signal at different SNR levels
SNRdB = [0 10 20];
Ps    = mean(Y.^2);
Pn    = Ps ./ (10.^(SNRdB/10));
xNoisy = Y + randn(Nt, numel(SNRdB)).*sqrt(Pn);

% Plotting time-domain signal for each SNR level
figure('Color','w');
t0 = 3001; L = 3000;
for i = 1:numel(SNRdB)
    subplot(numel(SNRdB),1,i);
    plot(t(t0:t0+L), xNoisy(t0:t0+L,i), 'LineWidth', 0.7);
    xlim([t(t0) t(t0+L)]); box on; grid off;
    xlabel('Time (s)'); ylabel('Amplitude');
end

% short-time Fourier transform (STFT)
stftOpt.win  = 1000;
stftOpt.ov   = 950;
stftOpt.nfft = 3000;
fr = [0 3000];

% Define tracking parameters
trk.startTime = 2;
trk.startFreq = 600;
trk.freqRange = 500;
trk.lambda    = 1;
trk.alpha     = 1500;
trk.rangeF    = 2;
trk.power     = 2;

% Time grid for frequency tracking
time0 = (0:Nt-1)'/fs;    % Time vector
dt    = 1/fs;            % Time step
winLen = max(3, min(round(0.01/dt), Nt));

rpm = nan(Nt, numel(SNRdB));

% Loop through each SNR level for frequency tracking
for i = 1:numel(SNRdB)
    y = xNoisy(:,i);

    % Compute the STFT for the noisy signal
    [S,F,Tspec] = stft(y, fs, ...
        Window=kaiser(stftOpt.win,5), ...
        OverlapLength=stftOpt.ov, ...
        FFTLength=stftOpt.nfft);

    mag = abs(S);   % Magnitude of the STFT
    FreqIdx = track_Freq(mag, F, Tspec, trk); 
    FreqF   = F(FreqIdx);

    % Plot the spectrogram with the tracked Freq
    figure('Color','w');
    plot_spectrogram(S,F,Tspec,stftOpt.win,fr);
    title(sprintf('Spectrogram (SNR = %d dB)', SNRdB(i)));
    xlim([0 T]); ylim(fr); colorbar; hold on;
    plot(Tspec, FreqF, 'g--', 'LineWidth', 1.5);

    % Smooth and interpolate the Freq frequencies to the original time grid
    FreqF = smoothdata(FreqF, 'sgolay', winLen);
    rpm(:,i) = interp1(Tspec, FreqF, time0, 'pchip', NaN);
end

% Plot RMSE results for frequency tracking
figure('Color','w'); hold on; box on; grid off;
colors = [0 0.45 0.74; 0.85 0.33 0.10; 0.47 0.67 0.19];

rmse = zeros(size(SNRdB));
refF = fk(:,3);   % Third-order harmonic as reference frequency

% Calculate RMSE for each SNR level
for i = 1:numel(SNRdB)
    valid = ~isnan(rpm(:,i)) & ~isnan(refF);
    err   = rpm(valid,i) - refF(valid);
    rmse(i) = sqrt(mean(err.^2));
    plot(t, rpm(:,i), '-.', 'Color', colors(i,:), 'LineWidth', 1.4);
end
plot(t, refF, 'k--', 'LineWidth', 1.6);

xlabel('Time (s)'); ylabel('Frequency (Hz)');
legend({'0 dB','10 dB','20 dB','Third-order harmonic'}, ...
    'Location','best','Box','off');
xlim([t(1) t(end)]);
ylim([min(refF)-100, max(refF)+100]);

% Display RMSE results
fprintf('\n=== RMSE Results ===\nSNR (dB)\tRMSE (Hz)\n-------------------------\n');
fprintf('%6.0f\t\t%8.3f\n', [SNRdB(:) rmse(:)]');

%% 2. Amplitude extraction simulation
% Perform VKF to estimate frequency components
[a,c] = vkf(xNoisy(:,1), fs, fk);
xEst  = real(a.*c);

% Plot VKF results for each frequency component
figure('Color','w','Units','pixels','Position',[100 100 1200 600]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

for k = 1:K
    ax = nexttile;
    [S,F,Tspec] = stft(xEst(:,k), fs, ...
        Window=kaiser(stftOpt.win,5), ...
        OverlapLength=stftOpt.ov, ...
        FFTLength=stftOpt.nfft);
    plot_spectrogram(S,F,Tspec,stftOpt.win,fr);
    ylim(fr); xlim([0 T]);
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    text(ax,0.02,0.95,sprintf('k=%d',k), ...
        'Units','normalized','FontName','Times New Roman', ...
        'FontSize',12,'FontWeight','bold','VerticalAlignment','top');
    ax.LooseInset = max(ax.TightInset,0.02);
end
% Plot the reconstructed spectrum
figure('Color','w');
[S,F,Tspec] = stft(sum(xEst,2), fs, ...
    Window=kaiser(stftOpt.win,5), ...
    OverlapLength=stftOpt.ov, ...
    FFTLength=stftOpt.nfft);
plot_spectrogram(S,F,Tspec,stftOpt.win,fr);
xlabel('Time (s)'); ylabel('Frequency (Hz)');
xlim([0 T]); ylim(fr); colorbar;

% Amplitude tracking evaluation (comparison between true and estimated amplitudes)
figure('Color','w'); hold on; grid off;
C = lines(K); 
for k = 1:K
    plot(t, Yamp(:,k), '-',  'Color', C(k,:), 'DisplayName', sprintf('A_%d',k));  % True amplitude
    plot(t, abs(a(:,k)), '--','Color', C(k,:), 'DisplayName', sprintf('B_%d',k)); % Estimated amplitude
end
ylim([0, 5.5]);
legend('Location','best','NumColumns',K,'Box','off');
xlabel('Time (s)'); ylabel('Amplitude');

for k = 1:K
    R = corrcoef(Yamp(:,k), abs(a(:,k)));
    fprintf('x_%d: r=%.4f, MSE=%.3e\n', k, R(1,2), mean((Yamp(:,k)-abs(a(:,k))).^2));  % Display results
end
