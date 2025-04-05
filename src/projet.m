%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fichier : projet.m
% Objet   : Traitement des données de turbulence – Questions 1 et 2
%
% Hypothèses :
% - Les fichiers de données sont nommés "signalXXX-026.dat" avec XXX allant
%   de 001 à 071, correspondant à la ligne j = 26.
% - dt = 2.5e-3 s, N = 1632 échantillons, dx = 1e-3 m.
% - K = 100, donc τ ∈ [-K*dt, K*dt].
%
% Énoncé :
% Q1 : Calculer les coefficients de corrélation temporelle R_{0i}(τ) entre
%      le point de référence (i0 = 37, j = 26) et tous les points de la ligne
%      (1 ≤ i ≤ 71) pour τ ∈ [-K*dt, K*dt] et tracer quelques courbes.
%
% Q2 : Tracer les isocontours de R_{0i}(r,τ) (avec r = (i-i0)*dx en abscisse et τ en ordonnée),
%      extraire pour chaque point le décalage τ_max pour lequel R_{0i}(τ) est maximal
%      (pour τ ≥ 0) et déterminer la vitesse de convection U_c via un ajustement linéaire
%      de la forme r = U_c * τ_max. Ici, RANSAC est utilisé pour extraire la tendance dominante.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% Paramètres
dt = 2.5e-3;         % Pas temporel (2.5 ms)
N  = 1632;           % Nombre d'échantillons par série
dx = 1e-3;           % Pas spatial (1 mm)
K = 100;             % Nombre de pas pour τ (donc τ ∈ [-K*dt, K*dt])
tauList = (-K:K)*dt; % Vecteur des décalages (en s)
nTau = length(tauList);

% Points spatiaux sur la ligne j = 26
i0 = 37;             % Point de référence
iList = 1:71;        % Indices des points sur la ligne
nPoints = length(iList);
rVector = (iList - i0)*dx;  % r = (i - i0)*dx (en m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 1 : Calcul et tracé des courbes de corrélation R_{0i}(τ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Chargement du signal de référence (i0 = 37)
refFile = sprintf('../data/signaux/signal%03d-026.dat', i0);
[uRef, ~] = load_velocity(refFile);
uRefFluc = uRef - mean(uRef);  % u'_0 = u0 - ̄u₀

if length(uRef) ~= N
    error('Le fichier de référence ne contient pas %d échantillons.', N);
end

% 2) Calcul des coefficients de corrélation R_{0i}(τ)
CorrMatrix = zeros(nPoints, nTau);  % Matrice (nPoints x nTau)

for idx = 1:nPoints
    iCur = iList(idx);
    
    % Chargement du signal du point iCur
    targetFile = sprintf('../data/signaux/signal%03d-026.dat', iCur);
    [uTarget, ~] = load_velocity(targetFile);
    
    if length(uTarget) ~= N
        error('Le fichier %s ne contient pas %d échantillons.', targetFile, N);
    end
    
    % Calcul de la fluctuation u'_i = u_i - ̄uᵢ
    uTargetFluc = uTarget - mean(uTarget);
    
    % Calcul pour chaque décalage τ = k*dt, pour k = -K ... K.
    for k = -K:K
        col = k + K + 1;  % Transformation de k = -K...K en indice 1 ... 2K+1
        if k >= 0
            nMax = N - k - 1;  % Borne supérieure corrigée
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

% 3) Tracé des courbes de corrélation pour quelques points choisis
pointsToPlot = [10, 26, 37, 60];  % Exemples d'indices parmi iList
figure;
hold on;
for iVal = pointsToPlot
    idx = find(iList == iVal);
    plot(tauList, CorrMatrix(idx, :), 'LineWidth', 1.5, 'DisplayName', sprintf('i = %d', iVal));
end
hold off;
xlabel('\tau (s)');
ylabel('Coefficient de corrélation R_{0i}(\tau)');
title('Courbes de corrélation pour quelques points (Question 1)');
legend('show'); grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 2 : Tracé des isocontours, extraction de la tendance dominante et estimation de U_c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Tracé des isocontours du champ de corrélation R_{0i}(r,τ)
figure;
contourf(rVector, tauList, CorrMatrix', 20, 'LineColor','none');
colorbar;
xlabel('Distance r (m)');
ylabel('\tau (s)');
title('Isocontours de R_{0i}(r,\tau) (Question 2)');

% 2) Extraction de τ_max pour chaque point (considérant τ ≥ 0)
tauPos = tauList(K+1:end);  % On ne prend que τ ≥ 0 (indice K+1 correspond à τ = 0)
tauMax = zeros(nPoints,1);
for idx = 1:nPoints
    corrPos = CorrMatrix(idx, K+1:end);  % Partie pour τ ≥ 0
    [~, imax] = max(corrPos);
    tauMax(idx) = tauPos(imax);
end

% 3) Tracé du nuage de points (τ_max vs. r)
figure;
scatter(tauMax, rVector, 50, 'filled');
xlabel('\tau_{max} (s)');
ylabel('Distance r (m)');
title('Nuage de points (τ_{max} vs. r) (Question 2)');
grid on;

% 4) Utilisation de RANSAC pour extraire la tendance linéaire dominante
maxIter = 500;     % Nombre d'itérations
threshold = 0.002;   % Seuil de distance (à ajuster selon l'échelle de r)
minInliers = 5;      % Nombre minimal d'inliers pour valider un modèle

[bestM, bestB, inliersIdx] = ransacLineFit(tauMax, rVector, maxIter, threshold, minInliers);
fprintf('RANSAC: U_c = %.4f m/s (pente)\n', bestM);

% 5) Tracé du nuage de points avec la droite RANSAC
figure;
scatter(tauMax, rVector, 50, 'filled'); hold on;
scatter(tauMax(inliersIdx), rVector(inliersIdx), 50, 'g', 'filled'); % Inliers en vert
tauFit = linspace(min(tauMax), max(tauMax), 100);
rFit = bestM * tauFit + bestB;
plot(tauFit, rFit, 'r-', 'LineWidth', 2);
hold off;
xlabel('\tau_{max} (s)');
ylabel('Distance r (m)');
title('Ajustement RANSAC pour l''estimation de U_c (Question 2)');
legend('Tous les points', 'Cluter le plus long', sprintf('Droite RANSAC: U_c = %.4f m/s', bestM), 'Location', 'Best');
grid on;
