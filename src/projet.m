% Définition des fichiers
ref_file = '../data/signaux/signal037-026.dat';
target_file = '../data/signaux/signal026-026.dat';

% Fonction pour charger les données
function [u, v] = load_velocity(file)
    fid = fopen(file, 'r');
    if fid == -1
        error('Impossible d''ouvrir le fichier %s', file);
    end
    
    % Sauter les métadonnées (4 lignes + 1 ligne de variables)
    for i = 1:5
        fgetl(fid);
    end
    
    % Lire les données numériques
    data = fscanf(fid, '%f %f', [2, Inf]);
    fclose(fid);
    
    % Extraire les vitesses u et v
    u = data(1, :);
    v = data(2, :); 
end

% Charger les séries temporelles
[uref, vref] = load_velocity(ref_file);
[utarget, vtarget] = load_velocity(target_file);

% Nombre d'échantillons
N = length(uref);

% Méthode 1 : Calcul classique avec l'équation (1)
mean_uref = mean(uref);
mean_utarget = mean(utarget);
mean_vref = mean(vref);
mean_vtarget = mean(vtarget);

num_u_eq1 = sum((uref - mean_uref) .* (utarget - mean_utarget));
num_v_eq1 = sum((vref - mean_vref) .* (vtarget - mean_vtarget));

den_u_eq1 = sqrt(sum((uref - mean_uref).^2) * sum((utarget - mean_utarget).^2));
den_v_eq1 = sqrt(sum((vref - mean_vref).^2) * sum((vtarget - mean_vtarget).^2));

corr_u_eq1 = num_u_eq1 / den_u_eq1;
corr_v_eq1 = num_v_eq1 / den_v_eq1;

% Méthode 2 : Calcul via les équations (2), (3) et (4)
num_u_eq234 = sum(uref .* utarget) - (1/N) * sum(uref) * sum(utarget);
num_v_eq234 = sum(vref .* vtarget) - (1/N) * sum(vref) * sum(vtarget);

den_u_eq234 = sqrt((sum(uref.^2) - (1/N) * sum(uref)^2) * (sum(utarget.^2) - (1/N) * sum(utarget)^2));
den_v_eq234 = sqrt((sum(vref.^2) - (1/N) * sum(vref)^2) * (sum(vtarget.^2) - (1/N) * sum(vtarget)^2));

corr_u_eq234 = num_u_eq234 / den_u_eq234;
corr_v_eq234 = num_v_eq234 / den_v_eq234;

% Méthode 3 : Utilisation de corrcoef
corr_u_builtin = corrcoef(uref, utarget);
corr_v_builtin = corrcoef(vref, vtarget);

% Affichage des résultats
disp('Coefficient de corrélation temporelle pour u entre les points 37 et 26:');
disp(['Calcul avec corrcoef: ', num2str(corr_u_builtin(1,2))]);
disp(['Calcul avec équation (1): ', num2str(corr_u_eq1)]);
disp(['Calcul avec équations (2), (3), (4): ', num2str(corr_u_eq234)]);

disp('Coefficient de corrélation temporelle pour v entre les points 37 et 26:');
disp(['Calcul avec corrcoef: ', num2str(corr_v_builtin(1,2))]);
disp(['Calcul avec équation (1): ', num2str(corr_v_eq1)]);
disp(['Calcul avec équations (2), (3), (4): ', num2str(corr_v_eq234)]);

% 1. Calcul de la matrice de corrélation
corr_matrix = zeros(N, N);
for tau = 1:N
    for r = 1:N-tau
        if (r+tau) <= N
            num = sum((uref(r:N-tau) - mean(uref)) .* (utarget(r+tau:N) - mean(utarget)));
            den = sqrt(sum((uref(r:N-tau) - mean(uref)).^2) * sum((utarget(r+tau:N) - mean(utarget)).^2));
            if den ~= 0
                corr_matrix(r, tau) = num / den;
            else
                corr_matrix(r, tau) = 0;
            end
        end
    end
end

% 2. Calcul et correction de la vitesse de convection Uc
[~, max_indices] = max(corr_matrix, [], 2); 
valid_indices = find(max_indices > 0 & max_indices < N);  

if ~isempty(valid_indices)
    r_values = valid_indices;
    tau_values = max_indices(valid_indices);
    
    % Inversion du signe de tau pour obtenir une vitesse positive
    tau_eff = -tau_values;
    
    % Ajustement linéaire r = Uc * tau_eff
    p = polyfit(tau_eff, r_values, 1);
    Uc = p(1);
    
    disp(['Valeur estimée de la vitesse de convection Uc: ', num2str(Uc)]);
else
    disp('Impossible de calculer Uc : données insuffisantes.');
end

% 3. Tracé des isocontours
figure;
contourf(1:N, 1:N, corr_matrix, 20, 'LineColor', 'none');
colorbar;
xlabel('Distance r');
ylabel('Décalage temporel \tau');
title('Isocontours du coefficient de corrélation R(r, \tau)');
