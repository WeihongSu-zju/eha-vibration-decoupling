function frequencyIdx = track_Freq(mag, F, Tspec, trk)
% Frequency tracking on STFT magnitude.
%
% INPUTS:
%   mag       - STFT magnitude matrix (nF x nT)
%   F         - Frequency vector (Hz)
%   Tspec     - Time vector for STFT
%   trk       - Tracking parameters struct
%
% OUTPUT:
%   frequencyIdx  - 1 x nT index of the frequency bin at each time frame

% Initialize parameters
magMax = max(mag(:)) + eps;
[nF, nT] = size(mag); 
df = F(2) - F(1);
freqBins = max(1, round(trk.freqRange / df)); 

% Initialize distance matrix and frequency index
D = inf(nF, nT);
frequencyIdx = nan(1, nT);

% Find start time and frequency
[~, t0] = min(abs(Tspec - trk.startTime));  
[~, f0] = min(abs(F - trk.startFreq));      

D(f0, t0) = -trk.alpha * (mag(f0, t0) / magMax);  % Initial cost
frequencyIdx(t0) = f0;  % Start tracking from initial frequency

% Forward tracking
for k = (t0 + 1):nT
    prev = frequencyIdx(k - 1);  
    r = max(1, prev - freqBins) : min(nF, prev + freqBins);
    fd = (r - prev)' * df;
    pen = trk.lambda * abs(fd / df - trk.rangeF) .^ trk.power; 
    cost = D(prev, k - 1) + pen - trk.alpha * (mag(r, k) / magMax);
    [~, id] = min(cost);
    frequencyIdx(k) = r(id);
    D(r, k) = min(D(r, k), cost);
end

% Backward tracking
for k = (t0 - 1):-1:1
    nxt = frequencyIdx(k + 1);  
    r = max(1, nxt - freqBins) : min(nF, nxt + freqBins);
    fd = (r - nxt)' * df;
    pen = trk.lambda * abs(fd / df + trk.rangeF) .^ trk.power;
    cost = D(nxt, k + 1) + pen - trk.alpha * (mag(r, k) / magMax);
    [~, id] = min(cost);
    frequencyIdx(k) = r(id);
    D(r, k) = min(D(r, k), cost);
end
end
