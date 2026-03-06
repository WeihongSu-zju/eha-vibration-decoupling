function varargout = vkf(y, fs, f) 
%==========================================================================
% VKF Vector Kalman Filter with embedded defaults 
%
% INPUTS: 
%   y           - Time-domain signal 
%   fs          - Sampling frequency (Hz) 
%   f           - Frequency components matrix
%
% OUTPUTS: 
%   1. Filtered signal with frequency components.
%   2. Filtered signal and frequency components.
%   3. Filtered signal, frequency components, and scaling factor (r).
%==========================================================================

% Define bandwidth as a fraction of the sampling rate (kHz)
bw = fs / 1000; 

% Force column vector for input signal y
if size(y, 2) > size(y, 1) 
    y = y.'; 
end 

% Number of time samples and frequency orders
n_t = length(y); 
n_ord = size(f, 2); 

% Phase factor calculation
c = exp(2i * pi * cumsum(f, 1) / fs);

% Set parameters for phase factor and Q matrix
phi = (pi / fs) * bw; 
s = 0:2; 
sgn = (-1) .^ s; 
Q = zeros(3, 3); 
Q(1, :) = 1; 

% Construct Q matrix
for ii = 1:2 
    Q(ii + 1, :) = (s .^ (2 * (ii - 1))) .* sgn; 
end 

% Solve for q using Q matrix
bq = zeros(3, 1); 
bq(1) = 2^3; 
q = round(Q \ bq).'.*sgn; 

% Calculate r and check denominator conditions
k = 0:2; 
den = 2 * sum(q(:) .* cos(k(:) * phi)); 
num = sqrt(2) - 1; 
r = sqrt(num / den); 

% Error handling for ill-conditioned cases
if den <= 0 || r > sqrt(1 / (2 * q(1) * eps)) 
    error(generatemsgid('IllConditioned'), 'Ill-conditioned B: bandwidth too small.'); 
end 

% Construct difference equation matrix A and coefficients
coeffs = pascal(3, 1); 
coeffs = coeffs(end, :); 
A = spdiags(ones(n_t - 2, 1) * coeffs, 0:2, n_t - 2, n_t); 

% Matrix R based on r value
R = spdiags(r * ones(n_t, 1), 0, n_t, n_t); 

% AR matrix multiplication
AR = A * R; 

% Matrix B for system equations
B = AR' * AR + speye(n_t); 

% Handle multi-frequency orders
nn = n_t * n_ord; 
BB = kron(speye(n_ord), B); 

if n_ord > 1 
    % Construct upper triangular matrix for multi-frequency components
    blU = n_ord * (n_ord - 1) / 2; 
    I = zeros(n_t * blU, 1); 
    J = zeros(n_t * blU, 1); 
    V = complex(zeros(n_t * blU, 1), 0); 
    idx = 0; 
    t = (1:n_t).'; 
    
    for ki = 1:n_ord 
        for kj = ki + 1:n_ord 
            % Define row and column indices for off-diagonal matrix elements
            rows = (ki - 1) * n_t + t; 
            cols = (kj - 1) * n_t + t; 
            vv = conj(c(:, ki)) .* c(:, kj); 
            I(idx + 1:idx + n_t) = rows; 
            J(idx + 1:idx + n_t) = cols; 
            V(idx + 1:idx + n_t) = vv; 
            idx = idx + n_t; 
        end 
    end 
    
    % Construct final sparse matrix BB
    BBu = sparse(I, J, V, nn, nn); 
    BB = BB + BBu + BBu'; 
end 

% Compute the complex conjugate of the frequency components
cy = conj(reshape(c, nn, 1)) .* repmat(y, n_ord, 1); 

% Solve for the output signal
x = 2 * reshape(BB \ cy, n_t, n_ord); 

% Return desired outputs
switch nargout 
    case 1 
        varargout{1} = real(x .* c); % Filtered signal with frequency components 
    case 2 
        varargout{1} = x; 
        varargout{2} = c; % Filtered signal and frequency components 
    case 3 
        varargout{1} = x; 
        varargout{2} = c; 
        varargout{3} = r; % Filtered signal, frequency components, and scaling factor r 
    otherwise 
        varargout{1} = real(x .* c); % Default filtered signal with frequency components 
end 
end
