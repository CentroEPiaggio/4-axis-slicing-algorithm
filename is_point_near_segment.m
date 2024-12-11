function is_near = is_point_near_segment(A, B, C, epsilon)
    % Calcola il vettore AB e AC
    AB = B - A;
    AC = C - A;
    
    % Calcola il parametro t per la proiezione di C su AB
    t = dot(AC, AB) / norm(AB)^2;
    
    % Verifica se la proiezione cade sul segmento
    % if t < 0
    %     % Punto più vicino è A
    %     closest_point = A;
    % elseif t > 1
    %     % Punto più vicino è B
    %     closest_point = B;
    is_near = 0;
    if (t>0+epsilon && t<1-epsilon)
        % Proiezione sul segmento
        closest_point = A + t * AB;
        % Calcola la distanza tra C e il punto più vicino
        distance = norm(C - closest_point);
        % Verifica se la distanza è inferiore alla tolleranza
        is_near = (distance <= epsilon);
    end
end
