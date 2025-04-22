clear; clc;

% Constantes physiques
kappa = 0.41;
y_max = 500;
N = 1000; %Nombre de point de maillage
y_plus = linspace(0.01, y_max, N);
dy = y_plus(2) - y_plus(1);

% Initialisation
U_euler = zeros(1, N);
U_rk2   = zeros(1, N);
U_rk4   = zeros(1, N);
v_values = zeros(1, N); % pour tracer v(y^+)

% Fonction de résolution implicite de v = dU+/dy+
solve_v = @(y) fzero(@(v) (1 + kappa^2 * y.^2 .* abs(v)) .* v - 1, 0.5);

% Pré-calculer v(y^+) pour tous les y^+
for i = 1:N
    v_values(i) = solve_v(y_plus(i));
end

%% --- Méthode Euler ---
for i = 2:N
    U_euler(i) = U_euler(i-1) + v_values(i-1) * dy;
end

%% --- Runge-Kutta ordre 2 ---
for i = 2:N
    k1 = v_values(i-1);
    k2 = v_values(i);
    U_rk2(i) = U_rk2(i-1) + (dy / 2) * (k1 + k2);
end

%% --- Runge-Kutta ordre 4 ---
for i = 2:N
    % Interpolation pour valeurs intermédiaires
    y1 = y_plus(i-1);
    y2 = y1 + dy/2;
    y3 = y1 + dy;

    k1 = solve_v(y1);
    k2 = solve_v(y2);
    k3 = solve_v(y2);
    k4 = solve_v(y3);

    U_rk4(i) = U_rk4(i-1) + (dy / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

%% --- Solution exacte analytique ---
U_exact = (1/kappa) * ( ...
    (1 - sqrt(1 + 4 * (kappa * y_plus).^2)) ./ (2 * kappa * y_plus) + ...
    log(2 * kappa * y_plus + sqrt(1 + 4 * (kappa * y_plus).^2)) );

%% --- Tracé de v(y⁺) ---
figure;
plot(y_plus, v_values, 'k', 'LineWidth', 1.5);
xlabel('$y^+$','Interpreter','latex','FontSize',14);
ylabel('$\frac{dU^+}{dy^+}$','Interpreter','latex','FontSize',14);
title('Profil de $v(y^+) = \frac{dU^+}{dy^+}$','Interpreter','latex','FontSize',16);
grid on;

%% --- Tracé des profils U⁺ avec la méthode de Prandtl ---
figure;
plot(y_plus, U_exact, 'k', 'LineWidth', 2); hold on;
plot(y_plus, U_euler, '--b', 'LineWidth', 1.4);hold on;
plot(y_plus, U_rk2, '-.g', 'LineWidth', 1.6);
plot(y_plus, U_rk4, '-r', 'LineWidth', 2);
xlabel('y^+'); ylabel('U^+');
title('Tracé des profils U+, avec comparaison des différentes méthodes numériques');
legend('Exacte','Euler','RK2','RK4','Location','southeast');
grid on;

%% --- Tracé de l’écart relatif entre méthodes pour Prandtl ---
rel_err_euler = abs((U_rk4 - U_euler) ./ U_rk4) * 100;
rel_err_rk2   = abs((U_rk4 - U_rk2) ./ U_rk4) * 100;

figure;
semilogy(y_plus, rel_err_euler, '--b', 'LineWidth', 1.4); hold on;
semilogy(y_plus, rel_err_rk2, '-.g', 'LineWidth', 1.6);
xlabel('$y^+$','Interpreter','latex','FontSize',14);
ylabel('Erreur relative (%)','Interpreter','latex','FontSize',14);
title('Erreur relative par rapport à RK4','Interpreter','latex','FontSize',16);
legend({'Euler vs RK4','RK2 vs RK4'}, 'Location','northeast');
grid on;

%% --- Tracé des erreurs relatives par rapport à la solution analytique---
rel_err_euler = abs((U_euler - U_exact) ./ U_exact) * 100;
rel_err_rk2   = abs((U_rk2 - U_exact) ./ U_exact) * 100;
rel_err_rk4   = abs((U_rk4 - U_exact) ./ U_exact) * 100;

figure;
semilogy(y_plus, rel_err_euler, '--b', 'LineWidth', 1.4); hold on;
semilogy(y_plus, rel_err_rk2, '-.g', 'LineWidth', 1.6);
semilogy(y_plus, rel_err_rk4, '-r', 'LineWidth', 2);
xlabel('y^+'); ylabel('Erreur relative (%)');
title('Écarts relatifs par rapport à la solution exacte');
legend('Euler','RK2','RK4','Location','northeast');
grid on;

%% --- Tableau de comparaison des méthodes avec erreurs relatives pour Prandtl vs Analytique ---
valeurs_y = [50, 100, 200, 300, 400, 500];
fprintf('\n%6s | %10s | %10s | %10s | %10s || %12s | %12s | %12s\n', ...
        'y^+', 'U^+ exacte', 'Euler', 'RK2', 'RK4', ...
        'Err Euler (%)', 'Err RK2 (%)', 'Err RK4 (%)');
fprintf(repmat('-', 1, 96)); fprintf('\n');

for y_target = valeurs_y
    [~, idx] = min(abs(y_plus - y_target));
    
    ue = U_exact(idx);
    u1 = U_euler(idx);
    u2 = U_rk2(idx);
    u4 = U_rk4(idx);
    
    err1 = abs((u1 - ue)/ue)*100;
    err2 = abs((u2 - ue)/ue)*100;
    err4 = abs((u4 - ue)/ue)*100;
    
    fprintf('%6.0f | %10.4f | %10.4f | %10.4f | %10.4f || %12.4f | %12.4f | %12.4f\n', ...
        y_target, ue, u1, u2, u4, err1, err2, err4);
end

%% --- Estimation de B pour le modèle de Prandtl et solution analytique ---

valeurs_y = [200, 300, 400, 500];
fprintf('\n--- Estimation de la constante B ---\n');
fprintf('%6s | %12s | %12s || %10s | %10s\n', ...
        'y^+', 'U^+ RK4', 'U^+ exacte', 'B (RK4)', 'B (exacte)');
fprintf(repmat('-', 1, 64)); fprintf('\n');

for y_target = valeurs_y
    [~, idx] = min(abs(y_plus - y_target));
    
    U4 = U_rk4(idx);
    Ue = U_exact(idx);
    ln_y = log(y_plus(idx));
    
    B4 = U4 - (1/kappa) * ln_y;
    Be = Ue - (1/kappa) * ln_y;
    
    fprintf('%6.0f | %12.4f | %12.4f || %10.4f | %10.4f\n', ...
            y_target, U4, Ue, B4, Be);
end

%% --- Analyse complète de B(y⁺) pour RK4 et solution exacte (Prandtl et solution analytique) ---

% Sélection des points y⁺ > 30
indices_B = find(y_plus > 30);
y_B = y_plus(indices_B);
U4_B = U_rk4(indices_B);
Ue_B = U_exact(indices_B);

% Calcul de B(y⁺)
B_RK4 = U4_B - (1/kappa) * log(y_B);
B_exact = Ue_B - (1/kappa) * log(y_B);

% Tracé de B(y⁺)
figure;
plot(y_B, B_exact, '-k', 'LineWidth', 1.8); hold on;
plot(y_B, B_RK4, '--r', 'LineWidth', 1.8);
xlabel('y^+'); ylabel('B estimé');
title('Variation de la constante B selon y^+');
legend('Solution exacte', 'Méthode RK4', 'Location','best');
grid on;

% Calcul des statistiques
moy_B_exact = mean(B_exact);
std_B_exact = std(B_exact);

moy_B_RK4 = mean(B_RK4);
std_B_RK4 = std(B_RK4);

% Affichage des résultats statistiques
fprintf('\n--- Statistiques sur B estimé (y^+ > 30) ---\n');
fprintf('Méthode       | Moyenne(B) | Écart-type(B)\n');
fprintf('--------------|------------|---------------\n');
fprintf('Solution exacte |   %.4f    |     %.4f\n', moy_B_exact, std_B_exact);
fprintf('RK4            |   %.4f    |     %.4f\n', moy_B_RK4, std_B_RK4);

%% Résolution avec modèle de Van Driest (Euler, RK2, RK4) ---

% Paramètre Van Driest
A0 = 26;

% Fonction solve_v pour Van Driest
solve_v_vd = @(y) fzero(@(v) ...
    (1 + (kappa * y * (1 - exp(-y / A0)))^2 * abs(v)) * v - 1, 0.5);

% Initialisation
U_euler_vd = zeros(1, N);
U_rk2_vd   = zeros(1, N);
U_rk4_vd   = zeros(1, N);
v_values_vd = zeros(1, N);

% Calcul de v(y⁺) pour Van Driest
for i = 1:N
    v_values_vd(i) = solve_v_vd(y_plus(i));
end

%% --- Méthode Euler (Van Driest) ---
for i = 2:N
    U_euler_vd(i) = U_euler_vd(i-1) + v_values_vd(i-1) * dy;
end

%% --- Méthode RK2 (Van Driest) ---
for i = 2:N
    k1 = v_values_vd(i-1);
    k2 = v_values_vd(i);
    U_rk2_vd(i) = U_rk2_vd(i-1) + (dy/2) * (k1 + k2);
end

%% --- Méthode RK4 (Van Driest) ---
for i = 2:N
    y1 = y_plus(i-1);
    y2 = y1 + dy/2;
    y3 = y1 + dy;

    k1 = solve_v_vd(y1);
    k2 = solve_v_vd(y2);
    k3 = solve_v_vd(y2);
    k4 = solve_v_vd(y3);

    U_rk4_vd(i) = U_rk4_vd(i-1) + (dy / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

%% --- Tracé des profils U⁺ avec la méthode de Van Driest ---
figure;
plot(y_plus, U_exact, 'k', 'LineWidth', 2); hold on;
plot(y_plus, U_euler_vd, '--b', 'LineWidth', 1.4);hold on;
plot(y_plus, U_rk2_vd, '-.g', 'LineWidth', 1.6);
plot(y_plus, U_rk4_vd, '-r', 'LineWidth', 2);
xlabel('y^+'); ylabel('U^+');
title('Tracé des profils U+ avec la méthode de Van Driest, comparaison des différentes méthodes numériques');
legend('Exacte','Euler','RK2','RK4','Location','southeast');
grid on;

%% --- Tracé de l’écart relatif entre méthodes pour VD ---
rel_err_euler_vd = abs((U_rk4_vd - U_euler_vd) ./ U_rk4_vd) * 100;
rel_err_rk2_vd   = abs((U_rk4_vd - U_rk2_vd) ./ U_rk4_vd) * 100;

figure;
semilogy(y_plus, rel_err_euler_vd, '--b', 'LineWidth', 1.4); hold on;
semilogy(y_plus, rel_err_rk2_vd, '-.g', 'LineWidth', 1.6);
xlabel('$y^+$','Interpreter','latex','FontSize',14);
ylabel('Erreur relative (%)','Interpreter','latex','FontSize',14);
title('Erreur relative par rapport à RK4','Interpreter','latex','FontSize',16);
legend({'Euler vs RK4','RK2 vs RK4'}, 'Location','northeast');
grid on;

%% --- Tracé comparatif : Prandtl vs Van Driest ---
figure;
plot(y_plus, U_rk4, 'r--', 'LineWidth', 1.5); hold on;
plot(y_plus, U_rk4_vd, 'b-', 'LineWidth', 1.8);
xlabel('y^+'); ylabel('U^+');
title('Comparaison des profils U^+ (modèle de Van Driest vs Prandtl)');
legend('Prandtl (RK4)', 'Van Driest (RK4)', 'Location','southeast');
grid on;

%% --- Tracé des deux lois logarithmiques associées (B différents pour Prandtl et Van Driest) ---

% Calcul des constantes B pour chaque modèle
B_prandtl = mean(U_rk4(y_plus > 200) - (1/kappa) * log(y_plus(y_plus > 200)));
B_vandriest = mean(U_rk4_vd(y_plus > 200) - (1/kappa) * log(y_plus(y_plus > 200)));

% Lois logarithmiques
U_log_prandtl = (1/kappa) * log(y_plus) + B_prandtl;
U_log_vandriest = (1/kappa) * log(y_plus) + B_vandriest;

% Tracé complet : profils RK4 + deux lois logarithmiques + loi linéaire
U_linear = y_plus;
U_linear(y_plus > 30) = NaN;  % coupe la loi linéaire après y⁺ = 30

figure;
plot(y_plus, U_rk4, 'r-', 'LineWidth', 2); hold on;
plot(y_plus, U_rk4_vd, 'b-', 'LineWidth', 2);
plot(y_plus, U_log_prandtl, 'm-.', 'LineWidth', 1.6);
plot(y_plus, U_log_vandriest, 'k:', 'LineWidth', 1.8);
plot(y_plus, U_linear, 'g--', 'LineWidth', 1.5);
xlabel('y^+'); ylabel('U^+');
title('Profils RK4 et lois asymptotiques (linéaire + logarithmiques)');
legend('Prandtl (RK4)', 'Van Driest (RK4)', ...
       'Loi log (Prandtl)', 'Loi log (Van Driest)', ...
       'Loi linéaire', ...
       'Location','southeast');
grid on;

% Affichage des constantes B
fprintf('\n--- Constantes B pour les lois logarithmiques ---\n');
fprintf('B (Prandtl RK4)     = %.4f\n', B_prandtl);
fprintf('B (Van Driest RK4)  = %.4f\n', B_vandriest);


%% --- Tracé des longueurs de mélange : Prandtl vs Van Driest ---

lm_prandtl = kappa * y_plus;
lm_vandriest = kappa * y_plus .* (1 - exp(-y_plus / A0));

figure;
plot(y_plus, lm_prandtl, '--r', 'LineWidth', 1.5); hold on;
plot(y_plus, lm_vandriest, '-b', 'LineWidth', 2);
xlabel('y^+'); ylabel('\ell_m^+');
title('Longueurs de mélange : Prandtl vs Van Driest');
legend('Prandtl', 'Van Driest', 'Location','northwest');
grid on;

%% Modèle de Van Driest avec déficit logarithmique ---

% Données physiques
delta = 0.077;         % m
u_tau = 0.0021;        % m/s
nu = 1e-6;             % m²/s
Re_tau = u_tau * delta / nu;  % ≈ 161.7

% Grille en y⁺
N_deficit = 1000;
y_plus_def = linspace(0.01, Re_tau, N_deficit);
dy_def = y_plus_def(2) - y_plus_def(1);

% Fonction longueur de mélange Van Driest
ellm_vd = @(y) kappa * y .* (1 - exp(-y / A0));

% Initialisation
v_deficit = zeros(1, N_deficit);
U_deficit = zeros(1, N_deficit);

% Résolution de v(y⁺) avec déficit (point par point)
for i = 1:N_deficit
    y = y_plus_def(i);
    RHS = 1 - y / Re_tau;
    solve_v = @(v) (1 + ellm_vd(y)^2 * abs(v)) * v - RHS;
    v_deficit(i) = fzero(solve_v, 1);
end

% Intégration avec RK4
for i = 2:N_deficit
    y1 = y_plus_def(i-1);
    y2 = y1 + dy_def/2;
    y3 = y1 + dy_def;

    RHS1 = 1 - y1 / Re_tau;
    RHS2 = 1 - y2 / Re_tau;
    RHS3 = 1 - y3 / Re_tau;

    k1 = fzero(@(v) (1 + ellm_vd(y1)^2 * abs(v)) * v - RHS1, 1);
    k2 = fzero(@(v) (1 + ellm_vd(y2)^2 * abs(v)) * v - RHS2, 1);
    k3 = fzero(@(v) (1 + ellm_vd(y2)^2 * abs(v)) * v - RHS2, 1);
    k4 = fzero(@(v) (1 + ellm_vd(y3)^2 * abs(v)) * v - RHS3, 1);

    U_deficit(i) = U_deficit(i-1) + (dy_def / 6) * (k1 + 2*k2 + 2*k3 + k4);
end

%% --- Tracé seul du profil avec déficit ---
figure;
plot(y_plus_def, U_deficit, 'b', 'LineWidth', 2);
xlabel('y^+'); ylabel('U^+');
title('Profil U^+ avec modèle de Van Driest et déficit logarithmique');
grid on;

%% --- Comparaison avec les deux autres profils RK4 ---
figure;
plot(y_plus, U_rk4, '--r', 'LineWidth', 1.5); hold on;
plot(y_plus, U_rk4_vd, '-g', 'LineWidth', 1.5);
plot(y_plus_def, U_deficit, '-b', 'LineWidth', 2);
plot(y_plus, U_log_vandriest, 'k:', 'LineWidth', 1.8);
plot(y_plus, U_linear, 'g--', 'LineWidth', 1.5);
xlabel('y^+'); ylabel('U^+');
title('Profils RK4 et lois asymptotiques (linéaire + logarithmiques)');

xlabel('y^+'); ylabel('U^+');
title('Comparaison des profils RK4 : Prandtl, Van Driest et Van Driest avec déficit');
legend('Prandtl (RK4)', 'Van Driest (RK4)', 'Van Driest + déficit (RK4)', ...
    'Loi log (Prandtl)', 'Loi log (Van Driest)', 'Loi linéaire',...
       'Location', 'southeast');
grid on;