clear; close all; clc;

folderPath = genpath('BSpline_Bruno');
addpath(folderPath);

%% Read data. Get Anterior surface. Visualize.
data_path = "data/Human_20 yo RIEB15-1632_OD_data.xls";

M = readmatrix(data_path, 'Sheet', 'Centered and Aligned', 'Range', 'A:B');
figure(1); scatter(M(:,2), M(:,1)); title("Original plot of data"); xlim([-5,5]); ylim([-5,5]); % Plot to confirm

% Bottom in original data, is indented from suture, so I
% replicate from top side to smoothen out data and make consistent
ant = M(M(:, 2) > 0, :); % Filter out anterior
figure(2); scatter(ant(:,1), ant(:,2)); title("Anterior Raw"); % Plot to check

% Process data for algorithm (places anterior on top, optic axis is x-axis)
ant_new = sortrows(ant, 2);
X_data = ant_new(:,2); Y_data = -1*ant_new(:,1);

% Get unique data
[Y_unique, idx] = unique(X_data);
X_unique = Y_data(idx);

figure(3); scatter(X_unique, Y_unique, 'LineWidth', 2); title("Right = Anterior. Left = Posterior");

%% Fit splines
figure(4); hold on;

% Variables
k=7; % Order
nknots = 6; % Number of knots - automatically removes knots
fixknots = [];

% Function forced to have: 
%   s(min) = 0; s(max) = 0;
%   s(min) = infSlope; s(max) = -infSlope
infSlope = 50;

pntcon(1) = struct('p', 0, 'x', [min(X_unique) max(X_unique)], 'v', [0 0]);
pntcon(2) = struct('p', 1, 'x', [min(X_unique) max(X_unique)], 'v', [infSlope -infSlope]);

options = struct('animation', 1, ...
    'figure', 4, ...
    'waitbar', 1, ...
    'display', 1, ...
    'KnotRemoval', 'none', ...
    'd', 1, 'lambda', 1e-3, 'regmethod', 'discrete', ...
    'qpengine', '', ...
    'sigma', [], ...
    'pntcon', pntcon);

[pp ier] = BSFK(X_unique, Y_unique, k, nknots, fixknots, options)

%% Take pp and construct symbolic
L = pp.pieces;
K = pp.order;
symbolicArr = sym(zeros(L,1));
syms x

for i = 1:L
    poly_coef = pp.coefs(i,:);
    poly_args = (x-pp.breaks(i)) .^ (fliplr(0:K-1));

    poly = poly_args * poly_coef';
    symbolicArr(i) = poly;
end

% Automatically find bounds
boundArr = sym([]);
for i = 1:L
    % Check if symbolic piecewise crosses over 3

    eqn = symbolicArr(i) == 3;
    bound = double(vpa(solve(eqn, x)));

    realBound = bound(imag(bound) == 0); % Take real solutions
    realBound = realBound(realBound > pp.breaks(i) & realBound < pp.breaks(i+1)); % Make sure bound is within range of piecewise

       % Loop over all solutions
       for j = 1:length(realBound)
           bound = realBound(j);
           if bound < 0
               negBound = bound
           else
               posBound = bound
           end
       end
    end

% Manually find bound
% negBound = -1.94;
% posBound = 0.85041;

%% Plot symbolic arr
figure(5); hold on;
for i=1:L
    if (pp.breaks(i) < negBound && pp.breaks(i+1) > negBound)
        fplot(symbolicArr(i), [pp.breaks(i), negBound])
    end
    if (pp.breaks(i) < posBound && pp.breaks(i+1) > posBound)
        fplot(symbolicArr(i), [posBound, pp.breaks(i+1)])
    end
    
    if (pp.breaks(i) < negBound && pp.breaks(i+1) < negBound) || (pp.breaks(i) > posBound && pp.breaks(i+1) > posBound)
        fplot(symbolicArr(i), [pp.breaks(i), pp.breaks(i+1)])
    end
end

%% Find metrics

% Get curvature
curvatureArr = sym(zeros(L,1));
for i = 1:L
    % pre-calculating first/second order derivatives of x(t), y(t). x(t) =
    % t, so x'(t) = 1, x''(t) = 0
    y_p = 1;
    y_pp = 0;
    x_p = diff(symbolicArr(i), 1);
    x_pp = diff(symbolicArr(i), 2);

    % formula for a curvature of a curve parameterized by x(t), y(t) -
    % can be found on wikipedia
    curvatureArr(i) = (x_p*y_pp - y_p*x_pp) / (x_p^2 + y_p^2)^(3/2);
end

% Plot curvature
figure(6); hold on;
for i=1:L
    fplot(curvatureArr(i), [pp.breaks(i), pp.breaks(i+1)])
end

% Get smoothing energy
smth_energy = 0;
for i=1:L

    if (pp.breaks(i) < negBound && pp.breaks(i+1) > negBound)      
        smth_energy = smth_energy + double(vpa(vpaintegral(diff(curvatureArr(i), 1) ^ 2, pp.breaks(i), negBound)));
    end

    if (pp.breaks(i) < posBound && pp.breaks(i+1) > posBound)
        smth_energy = smth_energy + double(vpa(vpaintegral(diff(curvatureArr(i), 1) ^ 2, posBound, pp.breaks(i+1))));
    end
    
    if (pp.breaks(i) < negBound && pp.breaks(i+1) < negBound) || (pp.breaks(i) > posBound && pp.breaks(i+1) > posBound)
        smth_energy = smth_energy + double(vpa(vpaintegral(diff(curvatureArr(i), 1) ^ 2, pp.breaks(i), pp.breaks(i+1))));
    end
end


% Get bending energy
bend_energy = 0;
for i=1:L
    % First derivative
    y_p = 1;
    x_p = diff(symbolicArr(i), 1);

    % Numerator
    num = diff(y_p/x_p, 1);

    % Final Integral Expression
    y_pp_x = num / x_p;
    expr = y_pp_x ^ 2;

    % Evaluate integral from bounds
    if (pp.breaks(i) < negBound && pp.breaks(i+1) > negBound)      
        bend_energy = bend_energy + eval(vpaintegral(expr, pp.breaks(i), negBound));
    end

    if (pp.breaks(i) < posBound && pp.breaks(i+1) > posBound)
        bend_energy = bend_energy + eval(vpaintegral(expr, posBound, pp.breaks(i+1)));
    end
    
    if (pp.breaks(i) < negBound && pp.breaks(i+1) < negBound) || (pp.breaks(i) > posBound && pp.breaks(i+1) > posBound)
        bend_energy = bend_energy + eval(vpaintegral(expr, pp.breaks(i), pp.breaks(i+1)));
    end
end

% Get arclength
arclength = 0;
for i=1:L    
    temp_len = double(vpaintegral(sqrt(1^2 + diff(symbolicArr(i), 1)^2), pp.breaks(i), pp.breaks(i+1)));

    arclength = arclength + temp_len
end

% Get fit
y_fit = ppval(pp, X_unique);

fit = sqrt(mean((y_fit - Y_unique) .^ 2)) * 10^3;

% Print final results
format long g;
results = [smth_energy, bend_energy, arclength, fit] * 2
