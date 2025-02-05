
close all
clear 
clc


%% Oggetto STL 

TR = stlread("cilindro.stl");

fv.vertices = TR.Points;
fv.faces = TR.ConnectivityList;
figure;
trisurf(TR, 'FaceColor', 'cyan', 'EdgeColor', 'none'); 
camlight; lighting phong; title('Original Mesh');
axis equal
grid off
axis off



%% Distensione del volume cilindrico nel piano

x = fv.vertices(:,1);
y = fv.vertices(:,2);
z = fv.vertices(:,3);

theta = atan2(z, y); % Angolo polare
r = sqrt(y.^2 + z.^2); % Raggio cilindrico

% Aggiornamento dei vertici
r_medio = (max(r)+min(r))/2;
Points_flat = [x,theta*r_medio,(r-min(r))];
excursion = max(Points_flat(:,3))-min(Points_flat(:,3));

% Aggiornamento delle facce
new_faces = [];
for i = 1:size(fv.faces, 1)
    face = fv.faces(i, :);
    theta_face = theta(face);
    % Verifica se la faccia attraversa la linea di taglio
    angoli_sorted = sort(unique(theta));
    offset = max(0.05,angoli_sorted(end)-angoli_sorted(end-1));
    cutAng = pi-offset;
    % Condizione per non entrare nella linea di taglio
    condition = max(theta_face) > -cutAng &&  min(theta_face) > -cutAng && max(theta_face) < cutAng &&  min(theta_face) < cutAng;
    if condition
        new_faces = [new_faces; face];
    end
end

% Rimozione dei vertici che non sono più referenziati ad alcun triangolo
referencedVertices = unique(new_faces(:));
newVertices = Points_flat(referencedVertices, :);

% Crea una mappa per aggiornare gli indici delle facce
indexMap = zeros(size(Points_flat, 1), 1);
indexMap(referencedVertices) = 1:length(referencedVertices);
newFaces = indexMap(new_faces);
% Verifico se qualche faccia ha vertici fuori dai limiti
invalid_faces = find(any(newFaces == 0, 2));
if ~isempty(invalid_faces)
    warning('Facce con indici di vertici non validi trovate e rimosse.');
    newFaces(invalid_faces, :) = [];  % Rimuovi le facce non valide
end

fv_flat.faces = newFaces;
fv_flat.vertices = newVertices;


%% Visualizzazione del volume disteso (eventualmente con normali)

figure;
trisurf(triangulation(fv_flat.faces, fv_flat.vertices), 'FaceColor', [0.2, 0.6, 0.2], 'EdgeColor', 'black');
title('Volume disteso');
axis equal;
axis off
camlight; lighting phong;

% % Visualizza le relative normali
% TR_flat = triangulation(fv_flat.faces, fv_flat.vertices);
% faceCenters1 = incenter(TR_flat);
% faceNormals1 = -faceNormal(TR_flat);
% figure;
% trisurf(TR_flat, 'FaceColor', [0.2, 0.6, 0.2], 'EdgeColor', 'black');
% hold on;
% quiver3(faceCenters1(:, 1), faceCenters1(:, 2), faceCenters1(:, 3),faceNormals1(:, 1), faceNormals1(:, 2), faceNormals1(:, 3), 2, 'Color', 'b', 'LineWidth', 2); 
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% axis equal;
% grid off;
% axis off;
% hold off;

%% CHIUSURA DIRETTA DEI BUCHI - se fai questo passare direttamente a %% Salva e visualizza il nuovo file STL ricostruito

% Identifico i bordi ------------------------------------------------------------------------------------
faces = fv_flat.faces;  
edges_all = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
edges_canonical = sort(edges_all, 2);
[~, unique_idx, idx] = unique(edges_canonical, 'rows', 'stable');
edges_unique = edges_all(unique_idx, :);
edge_counts = accumarray(idx, 1);
bordi = edges_unique(edge_counts == 1, :);

% Passo 1: Raggruppa i bordi in loop chiusi -----------------------------------------
vertici_bordo = bordi(:); % Elenco di tutti i vertici nei bordi
loop_chiusi = {}; % Celle che conterranno i loop chiusi

% Ricostruisci i loop
while ~isempty(bordi)
    % Inizia un nuovo loop
    loop = bordi(1, :); % Prendi il primo bordo
    bordi(1, :) = []; % Rimuovilo dai bordi
    while true
        % Trova un bordo che condivida un vertice con l'ultimo vertice del loop
        idx = find(bordi(:, 1) == loop(end) | bordi(:, 2) == loop(end), 1);
        if isempty(idx)
            break; % Loop completato
        end
        % Aggiungi il bordo al loop
        next_edge = bordi(idx, :);
        bordi(idx, :) = []; % Rimuovi il bordo trovato
        % Determina l'ordine del vertice
        if next_edge(1) == loop(end)
            loop = [loop, next_edge(2)];
        else
            loop = [loop, next_edge(1)];
        end
    end
    loop_chiusi{end+1} = loop; % Salva il loop trovato
end

% Passo 2: Crea triangoli per chiudere i loop ----------------------------------------
nuove_facce = []; % Array per le nuove facce
for i = 1:length(loop_chiusi)
    loop = loop_chiusi{i};
    n = length(loop);
    % Fissa un vertice del loop e crea triangoli con i successivi
    for j = 2:n-1
        nuove_facce = [nuove_facce; loop(1), loop(j), loop(j+1)];
    end
end

% Passo 3: Aggiungi le nuove facce al volume originale ------------------------------
fv_flat.faces = [fv_flat.faces; nuove_facce];

TR_flat_closed = triangulation(fv_flat.faces,fv_flat.vertices);

% Visualizza il volume chiuso con le relative normali
faceCenters = incenter(TR_flat_closed);
faceNormals = faceNormal(TR_flat_closed);
figure;
patch('Vertices', TR_flat_closed.Points, 'Faces', TR_flat_closed.ConnectivityList,'FaceColor', [0.8, 0.8, 1.0], 'EdgeColor', 'k', 'FaceAlpha', 0.5);
hold on;
quiver3(faceCenters(:, 1), faceCenters(:, 2), faceCenters(:, 3),faceNormals(:, 1), faceNormals(:, 2), faceNormals(:, 3), 1.5, 'Color', 'r', 'LineWidth', 1.5); 
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
grid on;
hold off;

SOLID_FV.vertices = fv_flat.vertices;
SOLID_FV.faces = fv_flat.faces;


%% Isolamento della superficie superiore MODO1 (elimino facce di base)

% Calcolo delle normali delle facce
normali = faceNormal(triangulation(fv_flat.faces, fv_flat.vertices));
tolleranza = 0.01;
facce_Superiori = normali(:,3)>0-tolleranza;
facce_superiori = ~facce_Superiori;
%indicizzo
new_faces_sup = fv_flat.faces(facce_superiori, :);
vertici_referenziati = unique(new_faces_sup(:));
new_vertices_sup = fv_flat.vertices(vertici_referenziati, :);
indexMap = zeros(size(fv_flat.vertices, 1), 1);
indexMap(vertici_referenziati) = 1:length(vertici_referenziati);
new_faces_sup = indexMap(new_faces_sup);

%Aggiorno la struttura fv_flat
fv_flat.faces = new_faces_sup;
fv_flat.vertices = new_vertices_sup;


%% Isolamento della superficie superiore MODO2 (mantiene solo le facce superiori orizzontali)

% Calcolo delle normali delle facce
normali = faceNormal(triangulation(fv_flat.faces, fv_flat.vertices));
tolleranza = 0.01;
facce_nonSuperiori = normali(:,3)<1-tolleranza;
facce_superiori = ~facce_nonSuperiori;
%indicizzo
new_faces_sup = fv_flat.faces(facce_superiori, :);
vertici_referenziati = unique(new_faces_sup(:));
new_vertices_sup = fv_flat.vertices(vertici_referenziati, :);
indexMap = zeros(size(fv_flat.vertices, 1), 1);
indexMap(vertici_referenziati) = 1:length(vertici_referenziati);
new_faces_sup = indexMap(new_faces_sup);

%Aggiorno la struttura fv_flat
fv_flat.faces = new_faces_sup;
fv_flat.vertices = new_vertices_sup;

%aggiorno l'altezza (altimenti mi sposta all'altezza di base)
fv_flat.vertices(:,3) = fv_flat.vertices(:,3)+excursion;


%% Visualizzazione della superficie isolata

figure;
trisurf(triangulation(fv_flat.faces, fv_flat.vertices), 'FaceColor', [0.2, 0.6, 0.2], 'EdgeColor', 'black');
title('Superficie superiore isolata');
axis equal;
grid on
axis on
camlight; lighting phong;


%% Ricostruzione del volume parendo dalla superficie superiore MODO1 (per superfici semplici, senza troppi buchi)

SOLID_FV = surf2solid(fv_flat,'elevation',0);
TR_flat_closed = triangulation(SOLID_FV.faces, SOLID_FV.vertices);
stlwrite(TR_flat_closed,'cilindro_flat.stl');

% Visualizza il volume chiuso
figure;
trisurf(TR_flat_closed, 'FaceColor', 'cyan','EdgeColor', 'k');
title('Volume chiuso ricostruito');
axis equal;
grid off
axis off
camlight; lighting phong;


%% Ricostruzione del volume parendo dalla superficie superiore MODO2 (per superfici complesse, con trama)

% Identifico i bordi ------------------------------------------------------------------------------------
faces = fv_flat.faces;  
edges_all = [faces(:, [1, 2]); faces(:, [2, 3]); faces(:, [3, 1])];
edges_canonical = sort(edges_all, 2);
[~, unique_idx, idx] = unique(edges_canonical, 'rows', 'stable');
edges_unique = edges_all(unique_idx, :);
edge_counts = accumarray(idx, 1);
bordi = edges_unique(edge_counts == 1, :);

% Identifico se 2 edges consecutivi appartenenti a 2 faces costituiscono 1 edges di una terza face. 
lista_vertici_bordo=union(bordi(:,1),bordi(:,2));
bordi_intersecati=[];
tolleranza = 5e-2;  % tolleranza per la collinearità

for i=1:size(bordi,1)
    bordo=bordi(i,:);
    %coordinate
    A=fv_flat.vertices(bordo(1),:);
    B=fv_flat.vertices(bordo(2),:);
    altri_vertici=lista_vertici_bordo;
    index=[find(altri_vertici==bordo(1));find(altri_vertici==bordo(2))];
    altri_vertici(index)=[];
    for k=1:size(altri_vertici,1)
        V=fv_flat.vertices(altri_vertici(k),:);
        if (is_point_near_segment(A, B, V, tolleranza))
            bordi_intersecati=[bordi_intersecati;i];
            altri_bordiA = find((bordi(:,1)==bordo(1)  & bordi(:,2)==altri_vertici(k)) | (bordi(:,1)==altri_vertici(k) & bordi(:,2)==bordo(1)));
            altri_bordiB = find((bordi(:,1)==bordo(2)  & bordi(:,2)==altri_vertici(k)) | (bordi(:,1)==altri_vertici(k) & bordi(:,2)==bordo(2)));
            altri_bordi = union(altri_bordiA,altri_bordiB)';
            bordi_intersecati=[bordi_intersecati;altri_bordi];
            break;
        end
    end
end

bordi_intersecati=unique(bordi_intersecati);
bordi(bordi_intersecati,:)=[];

% % plotto i bordi
figure;
trisurf(triangulation(fv_flat.faces, fv_flat.vertices), 'FaceColor', [0.2, 0.6, 0.2], 'EdgeColor', 'black');
title('Superficie superiore isolata');
axis equal;
grid off
axis off
camlight; lighting phong;
hold on
%bordi1=bordi(bordi_intersecati,:);
% for i = 1:size(bordi, 1)
%     v1 = fv_flat.vertices(bordi(i, 1), :);
%     v2 = fv_flat.vertices(bordi(i, 2), :);
%     plot3([v1(1), v2(1)], [v1(2), v2(2)], [v1(3), v2(3)], 'r-', 'LineWidth', 2.5);
%     hold on;
%     grid off
%     axis off
% end
%-----------------------------------------------------------------------------------------------------------------

% 1. Parametro di estrusione
distanza_estrusione = -excursion; 

% 2. Calcolo delle normali per ciascun vertice e vertici estrusi
normali_vere = faceNormal(triangulation(fv_flat.faces, fv_flat.vertices));
vertici_estrusi = fv_flat.vertices + distanza_estrusione * (normali_vere(1,:));

% 3. Costruzione delle facce laterali
facce_laterali = [];
num_vertici = size(fv_flat.vertices, 1);
for i = 1:size(bordi, 1)
    % Vertici del bordo inferiore
    v1 = bordi(i, 1);
    v2 = bordi(i, 2);
    
    % Vertici del bordo estruso
    v1_estruso = v1 + num_vertici;
    v2_estruso = v2 + num_vertici;
    
    % Crea due triangoli per ciascun lato del bordo
    facce_laterali = [facce_laterali; 
                      v1, v1_estruso, v2; 
                      v2, v1_estruso, v2_estruso];
end

% 4. Creazione della geometria estrusa completa
SOLID_FV.vertices = [fv_flat.vertices; vertici_estrusi];
facce_nuove = fv_flat.faces + num_vertici;
facce_nuove = flip(facce_nuove, 2); % inverti l'ordine dei vertici per invertire le normali
SOLID_FV.faces = [fv_flat.faces; facce_nuove;facce_laterali];
TR_flat_closed = triangulation(SOLID_FV.faces, SOLID_FV.vertices);

% Visualizza il volume chiuso
figure;
trisurf(TR_flat_closed, 'FaceColor', 'cyan', 'EdgeColor', 'black');
title('Volume chiuso ricostruito');
axis equal;
axis off
grid off
camlight; lighting phong;

% Visualizza il volume chiuso con le relative normali
faceCenters = incenter(TR_flat_closed);
faceNormals = faceNormal(TR_flat_closed);
figure;
patch('Vertices', TR_flat_closed.Points, 'Faces', TR_flat_closed.ConnectivityList,'FaceColor', [0.8, 0.8, 1.0], 'EdgeColor', 'k', 'FaceAlpha', 0.5);
hold on;
quiver3(faceCenters(:, 1), faceCenters(:, 2), faceCenters(:, 3),faceNormals(:, 1), faceNormals(:, 2), faceNormals(:, 3), 1.5, 'Color', 'r', 'LineWidth', 1.5); 
xlabel('X');
ylabel('Y');
zlabel('Z');
axis equal;
grid on;
hold off;


%% Salva e visualizza il nuovo file STL ricostruito

subplot(1, 2, 1);
trisurf(TR, 'FaceColor', 'cyan', 'EdgeColor', 'none'); 
camlight; lighting phong; title('Original Mesh');
axis equal
subplot(1, 2, 2);
trisurf(TR_flat_closed, 'FaceColor', 'magenta', 'EdgeColor', 'none');
camlight; lighting phong; title('Unrolled Mesh');
axis equal
grid off

stlwrite(TR_flat_closed,'volume_disteso_finale.stl');


%% Slicer estratto da IMAGObot

faces = fv_flat.faces;
vertex = fv_flat.vertices;
wallLineCountTop = 1;
TopLayers = 1;
WallLineCountBottom = 1;
BottomLayers = 1;
LayerHeightmm = 0.4;
LineWidthmm = 1; %0.4;
InfillDensity_2 = 70;
InfillDensity = 70;
NonPlanarOffsetmm = -0.1;
PlanarOffsetmm = 0;
InfillOverlap_2 = 50;
InfillOrientation_2 = 45;
RelativeInfillOrientation_2 = 00;
SpatialResolutionmm = 0.5;
WallLineCount = 1;
InfillOrientation = 45;
RelativeInfillOrientation = 90;
SafetyHeightmm = 12;
SamplesDensityDropDown = "Uniform";
NonPlanarPrintCheckBox = 0;

[TOT_number,TOP_number,TOP_nan,BOTTOM_number,BOTTOM_nan,PLANAR_number, PLANAR_nan,A_true, B_true, C_true]=IMAGOslicer(faces, vertex, ...
    wallLineCountTop, TopLayers,WallLineCountBottom,BottomLayers,LayerHeightmm,LineWidthmm, ...
    InfillDensity_2,InfillDensity,NonPlanarOffsetmm,PlanarOffsetmm,InfillOverlap_2, ...
    InfillOrientation_2,RelativeInfillOrientation_2,SpatialResolutionmm,WallLineCount,...
    InfillOrientation,RelativeInfillOrientation,SafetyHeightmm,SamplesDensityDropDown,NonPlanarPrintCheckBox);

if ~isempty(A_true)
    plot3(A_true(:,1),A_true(:,2),A_true(:,3),'Color',[0.9290, 0.6940, 0.1250],'LineWidth',1);
    hold on
    axis equal
    axis off
    grid off
end
if ~isempty(B_true)
    plot3(B_true(:,1),B_true(:,2),B_true(:,3),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1);
end
if ~isempty(C_true)
    plot3(C_true(:,1),C_true(:,2),C_true(:,3),'Color',[0, 0.4470, 0.7410],'LineWidth',1);
end

TopLayer=TOP_number;
PlanarLayer=PLANAR_number;
BottomLayer=BOTTOM_number;
TotLayer=TOT_number;

%% Creazione del file txt con i punti della traiettoria

[file,path] = uiputfile(strcat("*.txt"),strcat("Save .",string("txt")," file as:"));

                if(file == 0 & path == 0)
                    return;
                else

                    fid = fopen(strcat(path,file), 'wt');

                    TotLayer=[BottomLayer; PlanarLayer; TopLayer];

                    for i = 1:length(TotLayer)
                        line(i,:) = strcat(num2str(TotLayer(i,1)),...
                            ", ",num2str(TotLayer(i,2)),...
                            ", ",num2str(TotLayer(i,3)));
                    end
                    fprintf(fid, '%s\n', line);
                    fclose(fid);
                end


%% Importare il file traiettoria

filename = 'ramoebb.txt'; 
trajectory = readmatrix(filename, 'Delimiter', ',');

figure
plot3(trajectory(:,1),trajectory(:,2),trajectory(:,3))
axis equal
xlabel('x')
ylabel('y')
zlabel('z')


%% Punti da dare alla stampante ROTAZIONE IN MM

x = trajectory(:,1);
r_mandrino = 4; %[mm]
r = trajectory(:,3);
z = r;
r_perCalcolo = (r-min(r))+r_mandrino; %Normalizzo per risalire allo spazio percorso
theta_rad = trajectory(:,2)/r_medio; %così facendo, però, i punti hanno una fessura nella zona in cui ho fatto
                                     %il taglio, quindi vado a "tirarli" per chiudere il taglio. Scalo i valori 
                                     %da un intervallo [a,b] a un nuovo intervallo [c,d] con la seguente formula:
                                     % new_vect = [(d-c)/(b-a)]*(vect-a)+c
safe = (LineWidthmm/2)/r_medio*0.05;
c=-pi+safe;
d=pi-safe;
theta_rad_chiuso = (( d-c )/( max(theta_rad)-min(theta_rad) ))*(theta_rad-min(theta_rad))+(-pi);

for i=1:length(theta_rad_chiuso)
    theta_mm(i) = r_perCalcolo(i)*theta_rad_chiuso(i);
end

Points_stampa_mm = [x, theta_mm', z];
plot3(Points_stampa_mm(:,1),Points_stampa_mm(:,2),Points_stampa_mm(:,3))
axis equal


%% Punti da dare alla stampante ROTAZIONE IN DEG

r_mandrino = 4; %[mm]
x = trajectory(:,1);
r = trajectory(:,3); 
z = r;
r_perCalcolo = (r-min(r))+r_mandrino; %Normalizzo per risalire allo spazio percorso
theta_rad = trajectory(:,2)/r_medio; %così facendo, però, i punti hanno una fessura nella zona in cui ho fatto
                                     %il taglio, quindi vado a "tirarli" per chiudere il taglio. Scalo i valori 
                                     %da un intervallo [a,b] a un nuovo intervallo [c,d] con la seguente formula:
                                     % new_vect = [(d-c)/(b-a)]*(vect-a)+c
safe = (LineWidthmm/2)/r_medio*0.5;
c=-pi+safe;%+safe;
d=pi-safe;%-safe;
theta_rad_chiuso = (( d-c )/( max(theta_rad)-min(theta_rad) ))*(theta_rad-min(theta_rad))+(-pi);

for i=1:length(theta_rad_chiuso)
    theta_deg(i,:) = rad2deg(theta_rad_chiuso(i));
end

Points_stampa_deg = [x, theta_deg, z];
plot3(Points_stampa_deg(:,1),Points_stampa_deg(:,2),Points_stampa_deg(:,3))
axis equal


%% Avvolgere i punti planari attorno al cilindro (serve solo a visualizzare la geometria)

% riconverto in coordinate cartesiane originali
x_fin = trajectory(:,1);
theta_fin = trajectory(:,2)/(min(r)+4); 

safe = (LineWidthmm/2)/max(r)*0.05;
c=-pi+safe;
d=pi-safe;
theta_fin_chiuso = (( d-c )/( max(theta_fin)-min(theta_fin) ))*(theta_fin-min(theta_fin))+(-pi);

r_fin = trajectory(:,3)+8;
y_fin = 0.5*r_fin.*cos(theta_fin_chiuso);
z_fin = 0.5*r_fin.*sin(theta_fin_chiuso);
Points_cylindrical = [x_fin, y_fin, z_fin]; 
figure;
plot3(Points_cylindrical(:,1),Points_cylindrical(:,2),Points_cylindrical(:,3),'LineWidth',0.5)
grid off
axis equal

%Strato esterno non planare-----------------------------------------------
figure;
x_ext = B_true(:,1);
theta_ext = B_true(:,2)/(min(r)+4); 

safe = (LineWidthmm/2)/max(r)*0.05;
c=-pi+safe;
d=pi-safe;
theta_fin_chiuso = (( d-c )/( max(theta_ext)-min(theta_ext) ))*(theta_ext-min(theta_ext))+(-pi);

r_ext = B_true(:,3)+4;
y_ext = r_ext.*cos(theta_fin_chiuso);
z_ext = r_ext.*sin(theta_fin_chiuso);
Points_cylindrical = [x_ext, y_ext, z_ext]; 
hold on
axis off
axis equal
grid off
plot3(Points_cylindrical(:,1),Points_cylindrical(:,2),Points_cylindrical(:,3),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);


%Strato interno non planare-----------------------------------------------
x_ext = C_true(:,1);
theta_ext = C_true(:,2)/(min(r)+4); 

safe = (LineWidthmm/2)/max(r)*0.05;
c=-pi+safe*90;
d=pi-safe*90;
theta_fin_chiuso = (( d-c )/( max(theta_ext)-min(theta_ext) ))*(theta_ext-min(theta_ext))+(-pi);

r_ext = C_true(:,3)+4;
y_ext = r_ext.*cos(theta_fin_chiuso);
z_ext = r_ext.*sin(theta_fin_chiuso);
Points_cylindrical = [x_ext, y_ext, z_ext]; 
hold on
axis equal
axis off
grid off
plot3(Points_cylindrical(:,1),Points_cylindrical(:,2),Points_cylindrical(:,3),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1);


%Infill planare-----------------------------------------------------------
x_ext = A_true(:,1);
theta_ext = A_true(:,2)/min(r)+4; 

safe = (LineWidthmm/2)/max(r)*0.05;
c=-pi+safe;
d=pi-safe;
theta_fin_chiuso = (( d-c )/( max(theta_ext)-min(theta_ext) ))*(theta_ext-min(theta_ext))+(-pi);

r_ext = A_true(:,3)+4;
y_ext = r_ext.*cos(theta_fin_chiuso);
z_ext = r_ext.*sin(theta_fin_chiuso);
Points_cylindrical = [x_ext, y_ext, z_ext]; 
hold on
plot3(Points_cylindrical(:,1),Points_cylindrical(:,2),Points_cylindrical(:,3),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',0.5);
axis equal
axis off 
grid off


%% Seleziono tipo di controllo asse B (mm o deg)

Points_stampa = Points_stampa_deg;


%% Coordinate per G-code ESTRUSIONE A PISTONE

Points_stampa_conA = [Points_stampa(:,1),Points_stampa(:,3),zeros(length(Points_stampa(:,1)),1),Points_stampa(:,2)];  % X,Z,A,B
%Coordinata_spaziale_asseB = r_perCalcolo;
r_max = max(Points_stampa_conA(:,2)); %identifichiamo la safety heigh

% Per cambiare la safety heigh
%index=find(Points_stampa_conA(:,2)==r_max);
%Points_stampa_conA(index,2)=r_max*0.7; % 70% dell'altezza massima

referenceVolume = LayerHeightmm*LineWidthmm; %Sezione filamento depositato
volumeRef_siringa = pi*((12)^2)/4; %ml  
extrusorMoltiplicator = 1; %Parametro per ridurre o accentuare la deposizione totale di materiale
Poits_before_stopExtrusion = 5;
Poits_after_stopExtrusion = 2;

p=2;

while p < length(Points_stampa_conA(:,1))-5
    if Points_stampa(p+Poits_before_stopExtrusion,3) > (r_max-1)
        distance = sqrt( (Points_stampa_conA(p,1)-Points_stampa_conA(p-1,1))^2 + ( deg2rad((Points_stampa_conA(p,4)-Points_stampa_conA(p-1,4)))*r_perCalcolo(p) )^2 );
        extrusion_increment = (distance*referenceVolume/volumeRef_siringa)*extrusorMoltiplicator;
        Points_stampa_conA(p,3) = Points_stampa_conA(p-1,3) + extrusion_increment;
        for i=1:(Poits_before_stopExtrusion+Poits_after_stopExtrusion)
            Points_stampa_conA(p+i,3) = Points_stampa_conA(p+(i-1),3); %non estrudo altro materiale nei punti precedenti al 
                                                                       %distacco nè nei due punti successivi al distacco
        end
        p=p+(Poits_before_stopExtrusion+Poits_after_stopExtrusion);
    else
        distance = sqrt( (Points_stampa_conA(p,1)-Points_stampa_conA(p-1,1))^2 + ( deg2rad((Points_stampa_conA(p,4)-Points_stampa_conA(p-1,4)))*r_perCalcolo(p) )^2 );
        extrusion_increment = (distance*referenceVolume/volumeRef_siringa)*extrusorMoltiplicator;
        Points_stampa_conA(p,3) = Points_stampa_conA(p-1,3) + extrusion_increment;
    end
    p=p+1;
end
%Points_stampa_conA(end-2:end,:) = [];


%% Scrittura G-code ESTRUSIONE A PISTONE

vel = " F200";
vel_retraction = " F800";

[file,path] = uiputfile(strcat("*.ngc"),strcat("Save .","ngc"," file as:"));
fid = fopen(strcat(path,file), 'wt');

point_txt = [";Cylindrical mandrel 3D printing - EXTRUSION";
    ""];
point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(1,1)), " Z", num2str(Points_stampa_conA(1,2)), " A", num2str(Points_stampa_conA(1,3)), " B", num2str(Points_stampa_conA(1,4)), vel)];
i=2;
while i < size(Points_stampa_conA, 1)-3
    if Points_stampa_conA(i+1,2) > (r_max-1)
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i,1)), " Z", num2str(Points_stampa_conA(i,2)), " A", num2str(Points_stampa_conA(i,3)), " B", num2str(Points_stampa_conA(i,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+1,1)), " Z", num2str(Points_stampa_conA(i+1,2)), " A", num2str(Points_stampa_conA(i+1,3)), " B", num2str(Points_stampa_conA(i+1,4)),vel_retraction)];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+2,1)), " Z", num2str(Points_stampa_conA(i+2,2)), " A", num2str(Points_stampa_conA(i+2,3)), " B", num2str(Points_stampa_conA(i+2,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+3,1)), " Z", num2str(Points_stampa_conA(i+3,2)), " A", num2str(Points_stampa_conA(i+3,3)), " B", num2str(Points_stampa_conA(i+3,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+4,1)), " Z", num2str(Points_stampa_conA(i+4,2)), " A", num2str(Points_stampa_conA(i+4,3)), " B", num2str(Points_stampa_conA(i+4,4)), vel)];
        i=i+4;
    else
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i,1)), " Z", num2str(Points_stampa_conA(i,2)), " A", num2str(Points_stampa_conA(i,3)), " B", num2str(Points_stampa_conA(i,4)))];
    end
    i=i+1;
end
point_txt = [point_txt; "M2"];

for i = 1:size(point_txt,1)
    line(i,:) = point_txt(i);
end
fprintf(fid, '%s\n', line);
fclose(fid);


%% Coordinate per G-code FDM

Points_stampa_conA = [Points_stampa(:,1),Points_stampa(:,3),zeros(length(Points_stampa(:,1)),1),Points_stampa(:,2)];  % X,Z,A,B
%Coordinata_spaziale_asseB = r_perCalcolo;
r_max = max(Points_stampa_conA(:,2)); %identifichiamo la safety heigh

% Per cambiare la safety heigh
%index=find(Points_stampa_conA(:,2)==r_max);
%Points_stampa_conA(index,2)=r_max*0.7; % 70% dell'altezza massima

referenceSezione = LayerHeightmm*LineWidthmm; %Sezione della linea stampata [mm]^2
diametro_filamento = 1.75; %[mm]
sezione_filamento = pi*(diametro_filamento/2)^2; %[mm]^2
extrusorMoltiplicator = 0.95; % 3 %Parametro per ridurre o accentuare la deposizione totale di materiale (VALUTATO SPERIMENTALMENTE DA CURA)
retraction = 2;

p=2;

while p < length(Points_stampa_conA(:,1))-3
    if Points_stampa(p+1,3) > (r_max-1)
        distance = sqrt( (Points_stampa_conA(p,1)-Points_stampa_conA(p-1,1))^2 + ( deg2rad((Points_stampa_conA(p,4)-Points_stampa_conA(p-1,4)))*r_perCalcolo(p) )^2 );
        extrusion_increment = extrusorMoltiplicator*distance*referenceSezione/sezione_filamento;
        Points_stampa_conA(p,3) = Points_stampa_conA(p-1,3) + extrusion_increment;
        Points_stampa_conA(p+1,3) = Points_stampa_conA(p,3) - retraction;
        Points_stampa_conA(p+2,3) = Points_stampa_conA(p+1,3);
        Points_stampa_conA(p+3,3) = Points_stampa_conA(p+2,3) + retraction; % + extrusion_increment; %ripristino il valore e spingo
        p=p+3;
    else
        distance = sqrt( (Points_stampa_conA(p,1)-Points_stampa_conA(p-1,1))^2 + ( deg2rad((Points_stampa_conA(p,4)-Points_stampa_conA(p-1,4)))*r_perCalcolo(p) )^2 );
        extrusion_increment = extrusorMoltiplicator*distance*referenceSezione/sezione_filamento;
        Points_stampa_conA(p,3) = Points_stampa_conA(p-1,3) + extrusion_increment;
    end
    p=p+1;
end
%Points_stampa_conA(end-3:end,:) = [];


%% Scrittura G-code FDM

vel = " F300";
vel_retraction = " F900";

[file,path] = uiputfile(strcat("*.ngc"),strcat("Save .","ngc"," file as:"));
fid = fopen(strcat(path,file), 'wt');

point_txt = [";Cylindrical mandrel 3D printing - FDM";
    ""];
point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(1,1)), " Z", num2str(Points_stampa_conA(1,2)), " A", num2str(Points_stampa_conA(1,3)), " B", num2str(Points_stampa_conA(1,4)), vel)];
i=2;
while i < size(Points_stampa_conA, 1)-3
    if Points_stampa_conA(i+1,4)~=Points_stampa_conA(i,4) && Points_stampa_conA(i+1,1)==Points_stampa_conA(i,1)
        counter = 1;
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i,1)), " Z", num2str(Points_stampa_conA(i,2)), " A", num2str(Points_stampa_conA(i,3)), " B", num2str(Points_stampa_conA(i,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+1,1)), " Z", num2str(Points_stampa_conA(i+1,2)), " A", num2str(Points_stampa_conA(i+1,3)), " B", num2str(Points_stampa_conA(i+1,4)),vel_retraction)];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+1,1)), " Z", num2str(Points_stampa_conA(i+1,2)), " A", num2str(Points_stampa_conA(i+1,3)), " B", num2str(Points_stampa_conA(i+1,4)),vel)];
        i=i+1;
        % while Points_stampa_conA(i+counter+1,1)==Points_stampa_conA(i+counter,1) && Points_stampa_conA(i+1,1)==Points_stampa_conA(i,1)
        %     counter = counter+1;
        %     point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+counter,1)), " Z", num2str(Points_stampa_conA(i+counter,2)), " A", num2str(Points_stampa_conA(i+counter,3)), " B", num2str(Points_stampa_conA(i+counter,4)))];
        % end
        % counter = counter +1;
        % point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+counter,1)), " Z", num2str(Points_stampa_conA(i+counter,2)), " A", num2str(Points_stampa_conA(i+counter,3)), " B", num2str(Points_stampa_conA(i+counter,4)),vel)];
        % i=i+counter;
    end
    if Points_stampa_conA(i+1,2) > (r_max-1)
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i,1)), " Z", num2str(Points_stampa_conA(i,2)), " A", num2str(Points_stampa_conA(i,3)), " B", num2str(Points_stampa_conA(i,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+1,1)), " Z", num2str(Points_stampa_conA(i+1,2)), " A", num2str(Points_stampa_conA(i+1,3)), " B", num2str(Points_stampa_conA(i+1,4)),vel_retraction)];
        point_txt = [point_txt; strcat("G0 F6000 X", num2str(Points_stampa_conA(i+2,1)), " Z", num2str(Points_stampa_conA(i+2,2)), " B", num2str(Points_stampa_conA(i+2,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+3,1)), " Z", num2str(Points_stampa_conA(i+3,2)), " A", num2str(Points_stampa_conA(i+3,3)), " B", num2str(Points_stampa_conA(i+3,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+3,1)), " Z", num2str(Points_stampa_conA(i+3,2)), " A", num2str(Points_stampa_conA(i+3,3)), " B", num2str(Points_stampa_conA(i+3,4)), vel)];
        i=i+3;
    else
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i,1)), " Z", num2str(Points_stampa_conA(i,2)), " A", num2str(Points_stampa_conA(i,3)), " B", num2str(Points_stampa_conA(i,4)))];
    end
    i=i+1;
end
point_txt = [point_txt; "M2"];

for i = 1:size(point_txt,1)
    line(i,:) = point_txt(i);
end
fprintf(fid, '%s\n', line);
fclose(fid);



%% Scrittura G-code FDM

vel = " F300";
vel_retraction = " F900";

[file,path] = uiputfile(strcat("*.ngc"),strcat("Save .","ngc"," file as:"));
fid = fopen(strcat(path,file), 'wt');

point_txt = [";Cylindrical mandrel 3D printing - FDM";
    ""];
point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(1,1)), " Z", num2str(Points_stampa_conA(1,2)), " A", num2str(Points_stampa_conA(1,3)), " B", num2str(Points_stampa_conA(1,4)), vel)];
i=2;
while i < size(Points_stampa_conA, 1)-3
    if Points_stampa_conA(i+1,4)~=Points_stampa_conA(i,4) && Points_stampa_conA(i+1,1)==Points_stampa_conA(i,1)
        counter = 1;
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i,1)), " Z", num2str(Points_stampa_conA(i,2)), " A", num2str(Points_stampa_conA(i,3)), " B", num2str(Points_stampa_conA(i,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+1,1)), " Z", num2str(Points_stampa_conA(i+1,2)), " A", num2str(Points_stampa_conA(i+1,3)), " B", num2str(Points_stampa_conA(i+1,4)),vel_retraction)];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+1,1)), " Z", num2str(Points_stampa_conA(i+1,2)), " A", num2str(Points_stampa_conA(i+1,3)), " B", num2str(Points_stampa_conA(i+1,4)),vel)];
        i=i+1;
        % while Points_stampa_conA(i+counter+1,1)==Points_stampa_conA(i+counter,1) && Points_stampa_conA(i+1,1)==Points_stampa_conA(i,1)
        %     counter = counter+1;
        %     point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+counter,1)), " Z", num2str(Points_stampa_conA(i+counter,2)), " A", num2str(Points_stampa_conA(i+counter,3)), " B", num2str(Points_stampa_conA(i+counter,4)))];
        % end
        % counter = counter +1;
        % point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+counter,1)), " Z", num2str(Points_stampa_conA(i+counter,2)), " A", num2str(Points_stampa_conA(i+counter,3)), " B", num2str(Points_stampa_conA(i+counter,4)),vel)];
        % i=i+counter;
    end
    if Points_stampa_conA(i+1,2) > (r_max-1)
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i,1)), " Z", num2str(Points_stampa_conA(i,2)), " A", num2str(Points_stampa_conA(i,3)), " B", num2str(Points_stampa_conA(i,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+1,1)), " Z", num2str(Points_stampa_conA(i+1,2)), " A", num2str(Points_stampa_conA(i+1,3)), " B", num2str(Points_stampa_conA(i+1,4)),vel_retraction)];
        point_txt = [point_txt; strcat("G0 F6000 X", num2str(Points_stampa_conA(i+2,1)), " Z", num2str(Points_stampa_conA(i+2,2)), " B", num2str(Points_stampa_conA(i+2,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+3,1)), " Z", num2str(Points_stampa_conA(i+3,2)), " A", num2str(Points_stampa_conA(i+3,3)), " B", num2str(Points_stampa_conA(i+3,4)))];
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i+3,1)), " Z", num2str(Points_stampa_conA(i+3,2)), " A", num2str(Points_stampa_conA(i+3,3)), " B", num2str(Points_stampa_conA(i+3,4)), vel)];
        i=i+3;
    else
        point_txt = [point_txt; strcat("G1 X", num2str(Points_stampa_conA(i,1)), " Z", num2str(Points_stampa_conA(i,2)), " A", num2str(Points_stampa_conA(i,3)), " B", num2str(Points_stampa_conA(i,4)))];
    end
    i=i+1;
end
point_txt = [point_txt; "M2"];

for i = 1:size(point_txt,1)
    line(i,:) = point_txt(i);
end
fprintf(fid, '%s\n', line);
fclose(fid);
