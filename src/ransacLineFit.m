function [bestM, bestB, inliersIdx] = ransacLineFit(tau, r, maxIter, threshold, minInliers)
% ransacLineFit ajuste une droite de la forme r = m*tau + b via RANSAC.
% Entrées :
%   tau, r      : vecteurs (de même taille) représentant les données (en s et en m)
%   maxIter     : nombre maximal d'itérations (ex : 500)
%   threshold   : seuil de distance pour considérer un point comme inlier (ex : 0.002)
%   minInliers  : nombre minimal d'inliers pour valider le modèle (ex : 5)
%
% Sorties :
%   bestM, bestB : pente et intercept de la meilleure droite (U_c = bestM)
%   inliersIdx   : indices des points considérés comme inliers

    % Forcer tau et r en vecteurs colonnes
    tau = tau(:);
    r = r(:);
    N = length(tau);
    
    bestM = 0; bestB = 0; bestInliersCount = 0; inliersIdx = [];
    
    for iter = 1:maxIter
        % Choix aléatoire de 2 points distincts
        pts = randperm(N,2);
        t1 = tau(pts(1)); r1 = r(pts(1));
        t2 = tau(pts(2)); r2 = r(pts(2));
        
        if abs(t2 - t1) < 1e-14
            continue;
        end
        
        m = (r2 - r1) / (t2 - t1);
        b = r1 - m*t1;
        
        % Calcul de la distance verticale pour chaque point à la droite
        distAll = abs(r - (m*tau + b));
        inliers = find(distAll < threshold);
        inliersCount = numel(inliers);
        
        if inliersCount > bestInliersCount && inliersCount >= minInliers
            bestInliersCount = inliersCount;
            bestM = m;
            bestB = b;
            inliersIdx = inliers;
        end
    end
    
    % Réajustement sur les inliers pour affiner le modèle
    if bestInliersCount >= minInliers
        tauIn = tau(inliersIdx);
        rIn   = r(inliersIdx);
        p = polyfit(tauIn, rIn, 1);
        bestM = p(1);
        bestB = p(2);
    end
end
