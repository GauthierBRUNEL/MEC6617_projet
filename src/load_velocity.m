function [u, v] = load_velocity(filename)
% load_velocity lit un fichier de données de type "signalXXX-026.dat".
% On suppose que le fichier comporte 5 lignes d'en-tête à ignorer,
% puis 2 colonnes de données numériques (u et v).

    fid = fopen(filename, 'r');
    if fid == -1
        error('Impossible d''ouvrir le fichier %s', filename);
    end
    % Ignorer 5 lignes d'en-tête
    for i = 1:5
        fgetl(fid);
    end
    data = fscanf(fid, '%f %f', [2, Inf]);
    fclose(fid);
    data = data';
    u = data(:,1);
    v = data(:,2);
end
