clear; clc; close all;

if ~exist('../results', 'dir')
    mkdir('../results');
end

% ================================
% QUESTION 1
% ================================

dt = 2.5e-3;
N  = 1632;
dx = 1e-3;
K = 100;
tauList = (-K:K)*dt;
nTau = length(tauList);
i0 = 37;
iList = 1:71;
nPoints = length(iList);
rVector = (iList - i0)*dx;

refFile = sprintf('../data/signaux/signal%03d-026.dat', i0);
[uRef, ~] = load_velocity(refFile);
uRefFluc = uRef - mean(uRef);

if length(uRef) ~= N
    error('Le fichier de référence ne contient pas %d échantillons.', N);
end

CorrMatrix = zeros(nPoints, nTau);

for idx = 1:nPoints
    iCur = iList(idx);
    targetFile = sprintf('../data/signaux/signal%03d-026.dat', iCur);
    [uTarget, ~] = load_velocity(targetFile);
    if length(uTarget) ~= N
        error('Le fichier %s ne contient pas %d échantillons.', targetFile, N);
    end
    uTargetFluc = uTarget - mean(uTarget);
    for k = -K:K
        col = k + K + 1;
        if k >= 0
            nMax = N - k - 1;
            num = sum(uRefFluc(1:nMax) .* uTargetFluc(1+k:nMax+k));
            den = sqrt(sum(uRefFluc(1:nMax).^2) * sum(uTargetFluc(1+k:nMax+k).^2));
        else
            kk = abs(k);
            nMax = N - kk - 1;
            num = sum(uRefFluc(1+kk:nMax+kk) .* uTargetFluc(1:nMax));
            den = sqrt(sum(uRefFluc(1+kk:nMax+kk).^2) * sum(uTargetFluc(1:nMax).^2));
        end
        if den > 1e-14
            CorrMatrix(idx, col) = num / den;
        else
            CorrMatrix(idx, col) = 0;
        end
    end
end

pointsToPlot = [10, 26, 37, 60];
figure;
set(gcf, 'Visible', 'off')
hold on;
for iVal = pointsToPlot
    idx = find(iList == iVal);
    plot(tauList, CorrMatrix(idx, :), 'LineWidth', 0.5, 'DisplayName', sprintf('i = %d', iVal));
end
hold off;
xlabel('\tau (s)');
ylabel('Coefficient de corrélation R_{0i}(\tau)');
title('Courbes de corrélation pour quelques points (Question 1)');
legend('show'); grid on;
saveas(gcf, '../results/courbes_correlation.png');

% ================================
% QUESTION 2
% ================================

figure;
set(gcf, 'Visible', 'off')
contourf(rVector, tauList, CorrMatrix', 20, 'LineColor','none');
colorbar;
xlabel('Distance r (m)');
ylabel('\tau (s)');
title('Isocontours de R_{0i}(r,\tau) (Question 2)');
saveas(gcf, '../results/isocontours_corr.png');

tauMaxCenter = zeros(nPoints,1);
for iRow = 1:nPoints
    [~, idxMax] = max(CorrMatrix(iRow, :));
    tauMaxCenter(iRow) = tauList(idxMax);
end

Uc_center = sum(rVector(:).*tauMaxCenter(:)) / sum(tauMaxCenter(:).^2);
fprintf('Vitesse de convection estimée (raie centrale): U_c = %.10f m/s\n', Uc_center);

figure;
set(gcf, 'Visible', 'off')
contourf(rVector, tauList, CorrMatrix', 20, 'LineColor','none'); 
colorbar;
hold on;
rFit = linspace(min(rVector), max(rVector), 100);
tauTrend = rFit / Uc_center;
plot(rFit, tauTrend, 'r-', 'LineWidth', 2);
hold off;
xlabel('Distance r (m)');
ylabel('\tau (s)');
title(sprintf('Isocontours avec tendance centrale : τ = r/U_c, U_c = %.4f m/s', Uc_center));
legend('Tendance centrale', 'Location', 'Best');
grid on;
saveas(gcf, '../results/isocontours_corr_with_trend.png');

% ================================
% QUESTION 3
% ================================

I = 71;
J = 39;
Uc_center = 0.02;

data = readmatrix('../data/champs/champ0020.dat', ...
    'FileType', 'text', 'NumHeaderLines', 2);

X_raw = data(:,1);
Y_raw = data(:,2);
u_raw = data(:,3);
v_raw = data(:,4);

X = reshape(X_raw, [I, J])';
Y = reshape(Y_raw, [I, J])';
u = reshape(u_raw, [I, J])';
v = reshape(v_raw, [I, J])';

u_moving = u - Uc_center; 
v_moving = v;

nbPointsDomain = 20;
xRange = linspace(min(X(:)), max(X(:)), nbPointsDomain);
yRange = linspace(min(Y(:)), max(Y(:)), nbPointsDomain);
[Xstart_domain, Ystart_domain] = meshgrid(xRange, yRange);

nbPointsInlet = 50;
xLeft = min(X(:));
yStart = linspace(min(Y(:)), max(Y(:)), nbPointsInlet);
xStart_inlet = xLeft * ones(size(yStart));

figure;
set(gcf, 'Visible', 'off')
quiver(X, Y, u, v);
xlabel('X (m)'); ylabel('Y (m)');
title('Champ de vecteurs - Référentiel fixe');
axis equal; grid on;
saveas(gcf, '../results/champ_vect_ref_mes.png');

figure;
set(gcf, 'Visible', 'off')
verts_fixed_domain = stream2(X, Y, u, v, Xstart_domain, Ystart_domain);
hlines_fixed_domain = streamline(verts_fixed_domain);
set(hlines_fixed_domain, 'LineWidth', 0.5);
xlabel('X (m)'); ylabel('Y (m)');
title('Lignes de courant - Référentiel fixe');
axis equal; grid on;
saveas(gcf, '../results/streamlines_ref_mes_domain.png');

figure;
set(gcf, 'Visible', 'off')
verts_fixed_inlet = stream2(X, Y, u, v, xStart_inlet, yStart);
hlines_fixed_inlet = streamline(verts_fixed_inlet);
set(hlines_fixed_inlet, 'LineWidth', 0.5);
xlabel('X (m)'); ylabel('Y (m)');
title('Lignes de courant - Référentiel fixe (Inlet Only)');
axis equal; grid on;
saveas(gcf, '../results/streamlines_ref_mes_inlet.png');

figure;
set(gcf, 'Visible', 'off')
quiver(X, Y, u_moving, v_moving);
xlabel('X (m)'); ylabel('Y (m)');
title('Champ de vecteurs - Référentiel en translation');
axis equal; grid on;
saveas(gcf, '../results/champ_vect_ref_moving.png');

figure;
set(gcf, 'Visible', 'off')
verts_moving_domain = stream2(X, Y, u_moving, v_moving, ...
                              Xstart_domain, Ystart_domain);
hlines_moving_domain = streamline(verts_moving_domain);
set(hlines_moving_domain, 'LineWidth', 0.5);
xlabel('X (m)'); ylabel('Y (m)');
title('Lignes de courant - Référentiel en translation');
axis equal; grid on;
saveas(gcf, '../results/streamlines_ref_moving_domain.png');

figure;
set(gcf, 'Visible', 'off')
verts_moving_inlet = stream2(X, Y, u_moving, v_moving, ...
                             xStart_inlet, yStart);
hlines_moving_inlet = streamline(verts_moving_inlet);
set(hlines_moving_inlet, 'LineWidth', 0.5);
xlabel('X (m)'); ylabel('Y (m)');
title('Lignes de courant - Référentiel en translation (Inlet Only)');
axis equal; grid on;
saveas(gcf, '../results/streamlines_ref_moving_inlet.png');

% ================================
% QUESTION 4
% ================================

dx = mean(diff(X(1, :)));
dy = mean(diff(Y(:, 1)));

[dVdY, dVdX] = gradient(v, dy, dx);
[dUdY, dUdX] = gradient(u, dy, dx);

omega_z = dVdX - dUdY;

figure;
set(gcf, 'Visible', 'off')
contourf(X, Y, omega_z, 20, 'LineColor','none');
colorbar;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Composante de la vorticité ω_z');
saveas(gcf, '../results/vorticite.png');

nbPointsDomain = 20;
xRange = linspace(min(X(:)), max(X(:)), nbPointsDomain);
yRange = linspace(min(Y(:)), max(Y(:)), nbPointsDomain);
[Xstart_domain, Ystart_domain] = meshgrid(xRange, yRange);
verts_fixed = stream2(X, Y, u, v, Xstart_domain, Ystart_domain);

figure;
set(gcf, 'Visible', 'off')
contourf(X, Y, omega_z, 20, 'LineColor', 'none'); 
colorbar;
hold on;
hlines_fixed = streamline(verts_fixed);
set(hlines_fixed, 'LineWidth', 0.5, 'Color', [0.5 0.5 0.5]); 
hold off;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Superposition de la vorticité \omega_z et des lignes de courant (référentiel fixe)');
saveas(gcf, '../results/vorticite_streamlines.png');

% ================================
% QUESTION 5
% ================================

Q2D = (1/4) * (dUdY - dVdX).^2 ...
     - (1/2) * (dUdX.^2 + dVdY.^2 + (1/2)*(dUdY + dVdX).^2);

figure;
set(gcf, 'Visible', 'off')
contourf(X, Y, Q2D, 20, 'LineColor','none');
colorbar;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Critère Q_{2D}');
saveas(gcf, '../results/Q2D.png');
