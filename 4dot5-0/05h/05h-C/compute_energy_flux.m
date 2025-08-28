function [Pi, eps] = compute_energy_flux(u, v, dx, dy, filter_size)
% Compute SGS energy flux from 2D PIV data using Gaussian filtering.
% Inputs:
%   u, v: 2D velocity fields (size [Ny, Nx])
%   dx, dy: Grid spacing in x and y directions
%   filter_size: Standard deviation of Gaussian filter (in pixels)
% Output:
%   Pi: SGS energy flux field (size [Ny, Nx])
%   eps: non filtered field dissipation rate field (size [Ny, 1])

% Step 1: Apply Gaussian filter to velocities
u_filtered = imgaussfilt(u, filter_size);
v_filtered = imgaussfilt(v, filter_size);

% Step 2: Compute filtered strain rate tensor (S_ij)
[du_dx_filtered, du_dy_filtered] = gradient(u_filtered, dx, dy);
[dv_dx_filtered, dv_dy_filtered] = gradient(v_filtered, dx, dy);

S11_filtered = du_dx_filtered;
S12_filtered = 0.5 * (du_dy_filtered + dv_dx_filtered);
S22_filtered = dv_dy_filtered;

% Step 3: Compute tau_ij = filtered(u_i u_j) - filtered(u_i)*filtered(u_j)
uu = u .* u;
uv = u .* v;
vv = v .* v;

uu_filtered = imgaussfilt(uu, filter_size);
uv_filtered = imgaussfilt(uv, filter_size);
vv_filtered = imgaussfilt(vv, filter_size);

tau11_filtered = uu_filtered - u_filtered .* u_filtered;
tau12_filtered = uv_filtered - u_filtered .* v_filtered;
tau22_filtered = vv_filtered - v_filtered .* v_filtered;

% Step 4: Compute energy flux Pi = -tau_ij * S_ij
Pi = -(tau11_filtered .* S11_filtered + 2 * tau12_filtered .* S12_filtered + tau22_filtered .* S22_filtered);

% Step 5: Compute non filted dissipation rate epsilon
du_dy = gradient(mean(u, 2), dy);

S12 = du_dy;

tau12 = mean(uv, 2) - mean(u, 2) .* mean(v, 2);

eps = -1 * tau12 .* S12;
end