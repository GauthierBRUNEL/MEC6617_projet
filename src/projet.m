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
% Crée un fond neutre blanc pour "encadrer" le champ comme un contourf vide
contourf(X, Y, zeros(size(X)), 1, 'LineColor', 'none');
colormap([1 1 1]);  % fond blanc
hold on;
quiver(X, Y, u, v);
hold off;
xlabel('X (m)');
ylabel('Y (m)');
title('Champ de vecteurs - Référentiel fixe');
axis equal; grid on;
ylim([-10, 10]);  % borne verticale comme demandé
saveas(gcf, '../results/champ_vect_box_ref_mes.png');

figure;
set(gcf, 'Visible', 'off')
% On fait une "fausse" carte de fond blanche uniforme
contourf(X, Y, zeros(size(X)), 1, 'LineColor', 'none'); 
colormap([1 1 1]);  % blanc
hold on;
% Tracé des streamlines
streamslice(X, Y, u, v, 2);  % Densité contrôlée
hold off;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Lignes de courant - Référentiel fixe');
ylim([-10, 10]);
saveas(gcf, '../results/streamlines_box_ref_fixe.png');

figure;
set(gcf, 'Visible', 'off')
% Crée un fond neutre blanc pour "encadrer" le champ comme un contourf vide
contourf(X, Y, zeros(size(X)), 1, 'LineColor', 'none');
colormap([1 1 1]);  % fond blanc
hold on;
quiver(X, Y, u_moving, v_moving);
hold off;
xlabel('X (m)');
ylabel('Y (m)');
title('Champ de vecteurs - Référentiel en tranlation');
axis equal; grid on;
ylim([-10, 10]);  % borne verticale comme demandé
saveas(gcf, '../results/champ_vect_box_ref_moving.png');

figure;
set(gcf, 'Visible', 'off')
% On fait une "fausse" carte de fond blanche uniforme
contourf(X, Y, zeros(size(X)), 1, 'LineColor', 'none'); 
colormap([1 1 1]);  % blanc
hold on;
% Tracé des streamlines
streamslice(X, Y, u_moving, v_moving, 2);  % Densité contrôlée
hold off;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Lignes de courant - Référentiel en translation');
ylim([-10, 10]);
saveas(gcf, '../results/streamlines_box_ref_moving.png');

% ================================
% QUESTION 4
% ================================

dx = mean(diff(X(1, :)));
dy = mean(diff(Y(:, 1)));

% Méthode de différences centrées avec circshift
dvdx = (circshift(v, [0, -1]) - circshift(v, [0, 1])) / (2 * dx);
dudy = (circshift(u, [-1, 0]) - circshift(u, [1, 0])) / (2 * dy);
omega_z = dvdx - dudy;


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

figure;
set(gcf, 'Visible', 'off')
contourf(X, Y, omega_z, 20, 'LineColor', 'none'); 
colorbar;
hold on;
streamslice(X, Y, u, v, 2);  % 2 = densité modérée
hold off;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Superposition de la vorticité \omega_z et des lignes de courant (référentiel fixe)');
saveas(gcf, '../results/vorticite_streamlines.png');

% === Superposition vorticité + lignes de courant (référentiel en translation) ===
figure;
set(gcf, 'Visible', 'off')
contourf(X, Y, omega_z, 20, 'LineColor', 'none');
colorbar;
hold on;
streamslice(X, Y, u_moving, v_moving, 2);  % Densité contrôlée
hold off;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Superposition de la vorticité \omega_z et des lignes de courant (référentiel translation)');
saveas(gcf, '../results/vorticite_streamlines_moving.png');


% ================================
% QUESTION 5
% ================================

% Dérivées centrées
dudx = (circshift(u, [0, -1]) - circshift(u, [0, 1])) / (2 * dx);
dudy = (circshift(u, [-1, 0]) - circshift(u, [1, 0])) / (2 * dy);
dvdx = (circshift(v, [0, -1]) - circshift(v, [0, 1])) / (2 * dx);
dvdy = (circshift(v, [-1, 0]) - circshift(v, [1, 0])) / (2 * dy);

% Calcul du terme rotationnel (vorticité / 2)
Omega = 0.5 * (dvdx - dudy);

% Calcul du tenseur symétrique S^2
S2 = dudx.^2 + dvdy.^2 + 0.5 * (dudy + dvdx).^2;

% Formule de Q_{2D}
Q2D = Omega.^2 - S2;

% Affichage
figure;
set(gcf, 'Visible', 'off')
contourf(X, Y, Q2D, 20, 'LineColor','none');
colorbar;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Critère Q_{2D}');
saveas(gcf, '../results/Q2D.png');

figure;
set(gcf, 'Visible', 'off')
contourf(X, Y, Q2D, 20, 'LineColor','none');
colormap(parula); colorbar;
hold on;
streamslice(X, Y, u_moving, v_moving, 2);
hold off;
axis equal; grid on;
xlabel('X (m)');
ylabel('Y (m)');
title('Superposition de Q_{2D} et des lignes de courant (référentiel translation)');
saveas(gcf, '../results/Q2D_streamlines_moving.png');
