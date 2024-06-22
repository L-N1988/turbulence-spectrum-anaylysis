% This script processes `mat` from PIVLab toolbox 
% to calculate turbulence statistic parameters
clc; clear; close all;

%------------------------------------------------------------------------%
% read data
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
% clean data
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
% statistic process
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
% cross- or self- correlation, second moment
uv = zeros(my, nx); % u'v' Reynold shear stress
uu = zeros(my, nx); % u'u' TKE of u
vv = zeros(my, nx); % v'v' TKE of v
% third moment
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

% double average: average in time direction then in x direction
uv_xt = mean(uv ./ tCnt, 2); % u'v' Reynold shear stress
uu_xt = mean(uu ./ tCnt, 2); % u'u' TKE of u
vv_xt = mean(vv ./ tCnt, 2); % v'v' TKE of v
% rms means space average of stdandard deviation
u_rms = mean(sqrt(uu ./ tCnt), 2); % turbulence strength of u
v_rms = mean(sqrt(vv ./ tCnt), 2); % turbulence strength of v
uuu_xt = mean(uuu ./ tCnt, 2); % u'u'u' third moment of u

%------------------------------------------------------------------------%
% spectrum analysis
%------------------------------------------------------------------------%
center = [floor(my / 2) + 1, floor(nx / 2) + 1];
if contains(matPath, 'C')
    Fs = 250;
    window_len = [];
    spectrum_mat = 'pxx_f-C.mat';
elseif contains(matPath, 'L')
    Fs = 24;
    window_len = 16384;
    spectrum_mat = 'pxx_f-L.mat';
end
% Fs = 24;
X = zeros(lt, 1);
pxx = 0; f = 0;
for i = -1:1:1
	for j = -1:1:1
		for k = 1:lt
			X(k) = u_pri{k}(center(1) + i, center(2) + j);
		end
		[pxx_, f_] = pwelch(X, window_len, [], [], Fs);
        pxx = pxx + pxx_; f = f + f_;
	end
end
pxx = pxx ./ 9; f = f ./ 9;

save(fullfile(matPath, spectrum_mat), "f", "pxx", "u_rms", "U_xt");