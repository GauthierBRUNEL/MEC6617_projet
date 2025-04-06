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
%      (pour τ ≥ 0, en se limitant à une fenêtre primaire) et déterminer la vitesse de 
%      convection U_c via un ajustement linéaire forcé à passer par l'origine : r = U_c * τ_max.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

% Créer le dossier results s'il n'existe pas
if ~exist('../results', 'dir')
    mkdir('../results');
end

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
uRefFluc = uRef - mean(uRef);  % u'_0 = u0 - \bar{u}_0

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
    
    % Calcul de la fluctuation u'_i = u_i - \bar{u}_i
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
saveas(gcf, '../results/courbes_correlation.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 2 : Tracé des isocontours, extraction du maximum primaire et estimation de U_c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Tracé des isocontours du champ de corrélation R_{0i}(r,τ)
figure;
contourf(rVector, tauList, CorrMatrix', 20, 'LineColor','none');
colorbar;
xlabel('Distance r (m)');
ylabel('\tau (s)');
title('Isocontours de R_{0i}(r,\tau) (Question 2)');
saveas(gcf, '../results/isocontours_corr.png');

% 2) Extraction de τ_max pour chaque point, en se limitant à une fenêtre primaire
%    Ici, nous sélectionnons uniquement les τ ≤ tau_lim pour privilégier le maximum primaire.
tau_lim = 0.02;  % Limite pour le maximum primaire (par exemple, 20 ms)
tauPos = tauList(K+1:end);  % Valeurs pour τ ≥ 0
valid_idx = find(tauPos <= tau_lim);
if isempty(valid_idx)
    valid_idx = 1:length(tauPos);
end

tauMax = zeros(nPoints,1);
for idx = 1:nPoints
    corrPos = CorrMatrix(idx, K+1:end);  % Partie pour τ ≥ 0
    corrPrimary = corrPos(valid_idx);    % Restriction à la fenêtre primaire
    [~, imax] = max(corrPrimary);
    tauMax(idx) = tauPos(valid_idx(imax));
end

% 3) Tracé du nuage de points (τ_max vs. r)
figure;
scatter(tauMax, rVector, 50, 'filled');
xlabel('\tau_{max} (s)');
ylabel('Distance r (m)');
title('Nuage de points (τ_{max} vs. r) (Question 2)');
grid on;
saveas(gcf, '../results/nuage_points.png');

% 4) Ajustement linéaire forcé à passer par l'origine pour estimer U_c
%    On impose que r = U_c * τ_max (c'est-à-dire b = 0).
Uc = sum(rVector(:) .* tauMax(:)) / sum(tauMax(:).^2);
fprintf('Estimated convection speed U_c = %.4f m/s\n', Uc);

% 5) Tracé du nuage de points avec la droite d'ajustement forcée par l'origine
figure;
scatter(tauMax, rVector, 50, 'filled'); hold on;
tauFit = linspace(min(tauMax), max(tauMax), 100);
rFit = Uc * tauFit;  % Ajustement linéaire forcé: r = U_c * τ
plot(tauFit, rFit, 'r-', 'LineWidth', 2);
hold off;
xlabel('\tau_{max} (s)');
ylabel('Distance r (m)');
title('Ajustement linéaire (passant par 0) pour l''estimation de U_c (Question 2)');
legend('Données', sprintf('Droite: U_c = %.4f m/s', Uc), 'Location', 'Best');
grid on;
saveas(gcf, '../results/ajustement_lineaire.png');