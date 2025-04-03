% This script processes `mat` from PIVLab toolbox 
% to calculate turbulence statistic parameters
clc; clear; close all;

%------------------------------------------------------------------------%
%% read data
%------------------------------------------------------------------------%

% before running, add mat file path here!!!
matPath = '.\01h-L\';
data = load(strcat(matPath, 'PIVlab.mat'));
% vector, 2d plane coordinate values
Y = data.y{1, 1}(:, 1);
X = data.x{1, 1}(1, :);
% 3d matrix
u_original = data.u_original;
v_original = data.v_original;
if isempty(data.u_filtered{1}) && isempty(data.v_filtered{1})
    u_filtered = u_original;
    v_filtered = v_original;
else
    u_filtered = data.u_filtered;
    v_filtered = data.v_filtered;
end
plane_nomask = data.typevector_original{1, 1};
%    o--------o
%   /        /|
%  /        / |
% o--------o  |
% |        |  o
% |     my | /
% |        |/ lt
% o--------o
%     nx
% scalar, 3d matrix dimension length
[my, nx] = size(u_original{1}); % FIXME
lt = length(u_original); % snapshot number

tCnt = zeros(my, nx); % valid cell numbers in snapshots
uSum = zeros(my, nx);
vSum = zeros(my, nx);

%------------------------------------------------------------------------%
%% clean data
%------------------------------------------------------------------------%

for k = 1:lt
	for i = 1:my
		for j = 1:nx
			if plane_nomask(i, j)
				tCnt(i, j) = tCnt(i, j) + 1; % snapshot numbers of each valid cell
			else
				% remove masked cells
				u_filtered{k}(i, j) = 0;
				v_filtered{k}(i, j) = 0;
			end
		end
	end
end

%------------------------------------------------------------------------%
%% statistic process
% READ FIRST: https://www.mit.edu/course/1/1.061/www/dream/SEVEN/SEVENTHEORY.PDF
%------------------------------------------------------------------------%

for k = 1:lt
	uSum = uSum + u_filtered{k};
	vSum = vSum + v_filtered{k};
end

% time average velocity of each cell
% 2d matrix
U_t = uSum ./ tCnt;
V_t = vSum ./ tCnt;

% double average velocity
% vector
U_xt = mean(U_t, 2); % average in x direction
V_xt = mean(V_t, 2); % average in x direction

% double average component
% FIXME: this should be same as U_xt or V_xt, but not?
% U_xt = sum(uSum, 2) ./ sum(tCnt, 2);
% V_xt = sum(vSum, 2) ./ sum(tCnt, 2);

% turbulent velocity
u_pri = cell(lt, 1); % u'
v_pri = cell(lt, 1); % v'
%% cross- or self- correlation, second moment
uv = zeros(my, nx); % u'v' Reynold shear stress
uu = zeros(my, nx); % u'u' TKE of u
vv = zeros(my, nx); % v'v' TKE of v
%% third moment
uuu = zeros(my, nx); % u'u'u'

for k = 1:lt
	u_pri{k} = u_filtered{k} - U_t;
	v_pri{k} = v_filtered{k} - V_t;
	for i = 1:my
		for j = 1:nx
			if plane_nomask(i, j) == 0
				u_pri{k}(i, j) = 0;
				v_pri{k}(i, j) = 0;
			end
		end
	end
	uv = uv + u_pri{k} .* v_pri{k};
	uu = uu + u_pri{k}.^2;
	vv = vv + v_pri{k}.^2;
	uuu = uuu + u_pri{k}.^3;
end

%% double average: average in time direction then in x direction
uv_xt = mean(uv ./ tCnt, 2); % u'v' Reynold shear stress
uu_xt = mean(uu ./ tCnt, 2); % u'u' TKE of u
vv_xt = mean(vv ./ tCnt, 2); % v'v' TKE of v
%% rms means space average of stdandard deviation
u_rms = mean(sqrt(uu ./ tCnt), 2); % turbulence strength of u
v_rms = mean(sqrt(vv ./ tCnt), 2); % turbulence strength of v
uuu_xt = mean(uuu ./ tCnt, 2); % u'u'u' third moment of u

% assemable = [Y, U_xt, V_xt, uv_xt, uu_xt, vv_xt, u_rms, v_rms];
%------------------------------------------------------------------------%
%% spanwise correlation
% corr_ref_inx = [floor(nx / 2) +  1, floor(my / 2) + 1];
corr_ref_inx = [floor(nx / 2) +  1, 1];
corr_mat = cellfun(@(x) x .* x(corr_ref_inx(2), corr_ref_inx(1)), u_pri, 'UniformOutput', false);
self_corr_mat = cellfun(@(x) x .* x, u_pri, 'UniformOutput', false);
% ref: https://stackoverflow.com/a/42328979/18736354
corr_mat = mean(cat(3, corr_mat{:}), 3); 
self_corr_mat = mean(cat(3, self_corr_mat{:}), 3);
% normalized correlation matrix
corr_mat = corr_mat ./ sqrt(self_corr_mat .* self_corr_mat(corr_ref_inx(2), corr_ref_inx(1)));

figure();
imagesc(corr_mat'); colorbar;

figure();
plot(corr_mat(:, 4));

%% prepare data row in y direction, column in time direction
u_pri_yt = cellfun(@(x) x(:, floor(nx/2) + 1), u_pri, 'UniformOutput', false);
u_pri_yt = cat(2, u_pri_yt{:}); % convert to 2D u fluctuation matrixt, row is spanwise position, column is time

if contains(matPath, 'C')
    Fs = 150;
    % window_len = [];
elseif contains(matPath, 'L')
    Fs = 24;
    % window_len = 16384;
end

%% 1d spectrum analysis using direct method
dt = 1 / Fs;

f_ts = (1:lt-1) / lt / dt; % frequency domain vector
X = fft(u_pri_yt')';
GT = (1/lt/dt)*X(:,2:lt/2).* conj(X(:,2:lt/2));
f_ts = f_ts(2:lt/2);

figure();
loglog(f_ts, GT(10,:))
xlabel('f [Hz]');
ylabel('G_x(k_x: y=10 m)');

%% 1d spectrum analysis using welch method
% When x is a matrix, the PSD is computed independently 
% for each column and stored in the corresponding column of pxx. 
[pxx, fs] = pwelch(u_pri_yt', [], [], [], Fs); % column first
pxx = pxx';

figure();
loglog(fs, pxx(10, :));
xlabel('f [Hz]');
ylabel('G_x(k_x: y=10 m)');

%% 2d spectrum analysis using direct method
% Ref: https://ocean-physics.seos.uvic.ca/~jklymak/Phy580/html/Fft2d.html
%------------------------------------------------------------------------%
dt = 1 / Fs;
dy = abs(mean(diff(Y)));

f_ts = (1:lt-1) / lt / dt; % frequency domain vector
k_ys = (1:my-1) / my / dy; % wavenumber domain vector

YT = fft(fft(u_pri_yt')');
GYT = (1/lt/dt) * (1/my/dy) * YT(2:my/2,2:lt/2) .* conj(YT(2:my/2,2:lt/2));
k_ys = k_ys(2:my/2);
f_ts = f_ts(2:lt/2);

figure();

pcolor(f_ts,k_ys,real(GYT)); shading flat;
colorbar;
set(gca,'xscale','log','yscale','log')
xlabel('f [Hz]');
ylabel('k_y [cpm]');


xgrid = logspace(log10(min(f_ts)), log10(max(f_ts)), length(f_ts)*0.5);
ygrid = logspace(log10(min(k_ys)), log10(max(k_ys)), length(k_ys)*10);
[Xmesh,Ymesh] = meshgrid(xgrid, ygrid);
% Ref outer product: https://ww2.mathworks.cn/matlabcentral/answers/67282-function-which-returns-the-outer-product-of-two-vectors#comment_2623325
Zdata = real(GYT) .* (k_ys.' * f_ts);
Zmesh = griddata(f_ts, double(k_ys), double(Zdata), Xmesh, double(Ymesh));

con = figure('Position', [10 10 1000 618]);
% SEE: https://stackoverflow.com/a/44817243/18736354
% contourf(Xmesh, Ymesh, Zmesh, 20, 'LineStyle', 'none');
contour(Xmesh, Ymesh, Zmesh, 20);
% colormap(jet);
col = colorbar();
set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); set(gca, 'YScale', 'log'); 
set(xlabel("$f (\rm Hz)$"), 'Interpreter', 'latex'); 
set(ylabel("$k_y(\rm m^{-1})$"), 'Interpreter', 'latex');
set(ylabel(col, "$k_y fS_{uu} (\rm m^2/s^2)$"), 'Interpreter', 'latex');

%% 2d spectrum analysis using welch method (WOSA)
% Ref: https://ww2.mathworks.cn/matlabcentral/fileexchange/108894-pwelch2
dt = 1 / Fs;
dy = abs(mean(diff(Y)));

f_ts = (1:lt-1) / lt / dt; % frequency domain vector
k_ys = (1:my-1) / my / dy; % wavenumber domain vector

Pwosa = pwelch2(u_pri_yt, [floor(my/2)+1, floor(lt/2)+1], 0.6, @(x) hamming(x));

Pwosa_half = Pwosa(2:my/2, 2:lt/2);
k_ys = k_ys(2:my/2);
f_ts = f_ts(2:lt/2);

% figure();clf
% subplot(1,2,1);
% 
% pcolor(f_ts,k_ys,real(Pwosa_half));shading flat;
% set(gca,'xscale','log','yscale','log')
% xlabel('f [Hz]');
% ylabel('k_y [cpm]');
% 
% subplot(1,2,2);
% 
% surface(f_ts,k_ys,real(Pwosa_half));shading interp;
% set(gca,'xscale','log','yscale','log')
% xlabel('f [Hz]');
% ylabel('k_y [cpm]');
% view(30,30);
%% denoise
neps = 10;
% nmode = 8;
% smooth_window = {"gaussian", 200};
flim = 1;
noise_f = 0.025:0.025:flim; % 1.5:0转速比

valid = zeros(size(f_ts));
for i = 2:length(f_ts)
    eps = abs(f_ts(i) - f_ts(i-1)) * neps; % tolerance
    valid(i) = min(abs(noise_f - f_ts(i))) > eps;
end
valid = logical(valid);
f_denoise = f_ts(valid);
Pwosa_denoise = Pwosa_half(:, valid);

k_xs = 2*pi*f_denoise/mean(U_xt);
% [Xmesh,Ymesh] = meshgrid(k_xs,k_ys);
% pre_Pwosa = Pwosa_denoise .* Xmesh .* Ymesh; 
% Ref outer product: https://ww2.mathworks.cn/matlabcentral/answers/67282-function-which-returns-the-outer-product-of-two-vectors#comment_2623325
pre_Pwosa = Pwosa_denoise .* (k_ys.' * k_xs); % premultiplied spectrum

xgrid = logspace(log10(min(k_xs)), log10(max(k_xs)), length(k_xs)*0.1);
ygrid = logspace(log10(min(k_ys)), log10(max(k_ys)), length(k_ys)*100);
[Xmesh,Ymesh] = meshgrid(xgrid, ygrid);
Zmesh = griddata(double(k_xs), double(k_ys), double(pre_Pwosa), double(Xmesh), double(Ymesh));

%%
figure()

% contour(Xmesh, Ymesh, pre_Pwosa, 20);
contourf(Xmesh, Ymesh, Zmesh, 30, 'LineStyle','none');
colormap(jet);
col = colorbar();
clim([0, 30000]);
xlim([1, 500]);
ylim([0, 100])
set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); set(gca, 'YScale', 'log'); 
set(xlabel("$k_x (\rm m^{-1})$"), 'Interpreter', 'latex'); 
set(ylabel("$k_y (\rm m^{-1})$"), 'Interpreter', 'latex');
set(ylabel(col, "$k_x k_y S_{uu} (\rm m^2/s^2)$"), 'Interpreter', 'latex');