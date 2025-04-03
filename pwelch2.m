%% Function that implements the WOSA estimator in 2 dimensions
% Input parameters:
% - image: image of which to estimate the PSD.
% - window_size: window size, can be an integer or a
% vector, if window_size is an integer number the window used is square
% window_size * window_size, otherwise if window_size is a
% row vector the window is rectangular with width window_size(1, 2) and
% height window_size(1, 1)
% - overlap: overlap level of windows (0 <= overlap < 1)
% - window_handle: is the handle to the function for creating the window
% used to select the segments on which to estimate the DSP using the
% periodogram, an example is:
% eg.: window_handle = @(x) rectwin(x)
% to use a rectangular window
% eg.: window_handle = @(x) hamming(x)
% to use a Hamming window
%
% Output:
% - Pwosa: power spectral density estimated by WOSA method in 2
% dimensions
function Pwosa = pwelch2(image, window_size, overlap, window_handle)
    if overlap < 0 || overlap >= 1
        disp('error: overlap must be in range [0; 1)');
        return;
    end
    size_ = size(image);
    W = size_(2); % image width
    H = size_(1); % image height
    size_ = size(window_size); % window size
    if size_(2) == 1 % square window
        window_size_w = window_size; % window width
        window_size_h = window_size_w; % window width equals to window height
    elseif size_(2) == 2 % rectangulare window
        window_size_w = window_size(1, 2); % window width 
        window_size_h = window_size(1, 1); % window height
    end
    n_window_w = ceil(W / ceil((window_size_w * (1 - overlap)))); % number of windows along x direction
    n_window_h = ceil(H / ceil((window_size_h * (1 - overlap)))); % number of windows along y direction
    total_window = n_window_w * n_window_h; % number of total windows
    step_x = ceil(window_size_w * (1 - overlap)); % horizontal step
    step_y = ceil(window_size_h * (1 - overlap)); % vertical step
    k = 1.0 / total_window; % weight of a single segment 
    Pwosa = zeros(H, W); % 2D WOSA PSD result
    for i = 1:n_window_w
        for j = 1:n_window_h
            x_start = (i-1)*step_x + 1; % horizontal index where selection starts
            x_end   = x_start + window_size_w -1; % horizontal selection ends
            y_start = (j-1)*step_y + 1; % vertical index where selection starts
            y_end   = y_start + window_size_h -1; % horizontal selection ends
            % check if indexes exceeds image dimension
            if x_end <= W && y_end <= H 
                w = window_handle(window_size_h)*window_handle(window_size_w)'; % 2D window
                segment = image(y_start:y_end, x_start:x_end).*w; % segment multiplied window
            elseif x_end <= W && y_end > H 
                w = window_handle(window_size_h - (y_end - H))*window_handle(window_size_w)';
                segment = image(y_start:end, x_start:x_end).*w;
            elseif x_end > W && y_end <= H 
                w = window_handle(window_size_h)*window_handle(window_size_w - (x_end - W))';
                segment = image(y_start:y_end, x_start:end).*w;
            elseif x_end > W && y_end > H 
                w = window_handle(window_size_h - (y_end - H))*window_handle(window_size_w - (x_end - W))';
                segment = image(y_start:end, x_start:end).*w;
            end
            Ew    = sum(sum(w.^2)); % window's energy
            Pper  = 1/(W*H)*abs(fft2(segment, H, W)).^2; % periodogram
            Pwosa = Pwosa + Ew * k * Pper; % adding the current periodogram in Pwosa weighted per k
        end
    end
end