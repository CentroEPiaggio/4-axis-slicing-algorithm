function [TOT_number,TOP_number,TOP_nan,BOTTOM_number,BOTTOM_nan,PLANAR_number, PLANAR_nan,A_true, B_true, C_true]=IMAGOslicer(faces,points,nNonPlanarPerimetersUp,...
    nNonPlanarLayersUp,...
    nNonPlanarPerimetersDown,...
    nNonPlanarLayersDown,...
    layer_height,...
    extrusion_width,...
    percNonPlanar,...
    percPlanar,...
    off_nonplanar_small,...
    off_planar_small,...
    infill_overlap,...
    theta1NonPlanar,...
    theta_relative_nonPlanar,...
    spatial_resolution,...
    nPlanarPerimeters,...
    theta1Planar,...
    theta_relative_planar,...
    zRet,...
    sample_density,...
    non_planar)


    TOT_number=[];
    TOP_number=[];
    BOTTOM_number=[];
    PLANAR_number=[];
    A_true=[];
    B_true=[];
    C_true=[];
    temp_bottom=[];
    temp_planar=[];
    temp_top=[];
    
    
    %% Non Planar Up
    
    close all
    warning('off')
    
    %% Carico il file stl
    
    object=triangulation(faces,points);
    
    mediaX=mean(points(:,1));
    mediaY=mean(points(:,2));
    mediaZ=mean(points(:,3));
    
    %% Parametri di input
    
    thresUp=85; % soglia normali di top
    thresDown=95; % soglia normali di bottom
    
    off_nonplanar=-extrusion_width/2; %offset dei perimetri interni
    off_planar=-extrusion_width/2;
    
    dx=extrusion_width; % distanza tra i punti di controllo
    
    angle=180;
    
    
    %% Creo triangoli e vettore delle altezze dei layer
    
    index_tri = object.ConnectivityList; % indici dei triangoli della mesh
    tri = zeros(length(index_tri),9); % matrice dei triangoli della mesh
    
    for h=1:length(index_tri)
        k = 1;
        for j=1:3:size(tri,2)
            tri(h,j:j+2) = points(index_tri(h,k), :);
            k = k+1;
        end
    end
    
    z = ((min(points(:,3))):layer_height:max(points(:,3))); % divisione in layers
    
    TR=triangulation(index_tri,points); % mesh originale
    
    %% Head-tail
    
    percNonPlanar=percNonPlanar/100;
    percPlanar=percPlanar/100;
    infill_overlap=infill_overlap/100;
    
    [shell,triTot,triThresUp,triThresDown]=head_tail_true(tri,z,thresUp,thresDown,TR); % slicing
    
    percTriUp=triThresUp./triTot; % percentuale di triangoli superano l'angolo soglia
    percTriDown=triThresDown./triTot;
    
    polyShell=cell2poly(shell);
    
    for i=1:size(polyShell,2)
        if isempty(polyShell{1,i})
            polyShell{1,i}=polyshape(0,0);
        end
    end
    
    if size(polyShell{1,end}.Vertices,1)<3
        polyShell(:,end)=[];
        shell(:,end)=[];
        z(:,end)=[];
    end
    
    maximum_height=max(points(:,3));
    minimum_height=min(points(:,3));
    
    maximum_z=z(1,end);
    shift=maximum_height-maximum_z;
    height=maximum_height-minimum_height;
    
    %% Selezione del primo layer non planare e generazione dello spazio per le non planar surface
    
    if non_planar==1
    
        percThresUp=0;
        percThresDown=0;
    
        zFirst=find(percTriUp>=percThresUp,1);
        zEnd=size(polyShell,2);
    
        zThres=1;
        roomUp=cell(nNonPlanarLayersUp,size(shell,2));
    
        for n=1:nNonPlanarLayersUp
            roomUp{n,zEnd-(n-1)}=polyShell{1,zEnd-(n-1)};
            for i=zFirst:zEnd-n
                roomUp{n,i}=subtract(polyShell{1,i},polyShell{1,i+1});
                for j=i+1:zEnd
                    inters=intersect(roomUp{n,i},polyShell{1,j});
                    if size(inters.Vertices)>1
                        roomUp{n,i}=subtract(roomUp{n,i},inters);
                    end
                end
                polyShell{1,i}=subtract(polyShell{1,i},roomUp{n,i});
            end
            polyShell{1,zEnd-(n-1)}=subtract(polyShell{1,zEnd-(n-1)},roomUp{n,zEnd-(n-1)});
        end
    
        shell1=cell(size(shell,1),size(shell,2));
    
        for i=1:size(polyShell,2)
            if ~isempty(polyShell{1,i})
                [xt,yt]=boundary(polyShell{1,i});
                shell1{1,i}=[xt yt z(i)*ones(size(boundary(polyShell{1,i}),1),1)];
            end
        end
    
        shell=shell1;
    
        %% Crop della mesh dal layer selezionato in su
    
        % per ogni layer non planare faccio un crop della mesh e proietto la
        % silohuette in 2D per ottenere il perimetro non planare fluttuante
    
        polyFloat=cell(1,nNonPlanarLayersUp);
    
        for n=1:nNonPlanarLayersUp
            for i=1:size(roomUp,2)
                if ~isempty(roomUp{1, i})
                    if ~isempty(roomUp{n, i})
                        if ~isempty(roomUp{n, i}.Vertices)
                            zThres=i+n-1; % primo layer
                            break
                        end
                    end
                end
            end
    
            if zThres<=1
                TR_new=TR;
                indexTri=[];
            else
    
                % costruisco una nuova mesh con solo i triangoli sopra il primo layer in
                % modo da calcolarmi il perimetro della surface non planare
    
                index_new=index_tri;
    
                ind1=find(points(:,3)<z(zThres));
                temp=[];
    
                for i=1:size(ind1,1)
                    [i1,~]=find(index_new==ind1(i));
                    temp=[temp;i1];
                end
    
                for i=1:size(temp,1)
                    index_new(temp(i),:)= [nan nan nan];
                end
    
                index_new(any(isnan(index_new), 2), :) = [];
    
                TR_new=triangulation(index_new,points); % mesh croppata con i triangoli INTERAMENTE sopra il primo layer
                indexTri=find(tri(:,3)<z(zThres) & tri(:,6)>z(zThres) | tri(:,3)>z(zThres) & tri(:,6)<z(zThres) | tri(:,6)<z(zThres) & tri(:,9)>z(zThres)|tri(:,6)>z(zThres) & tri(:,9)< z(zThres) |tri(:,3)>z(zThres) & tri(:,9)< z(zThres)|tri(:,3)<z(zThres) & tri(:,9)>z(zThres));
    
            end
    
            if ~isempty(indexTri)
    
                trUp=tri(indexTri,:); % triangoli che sono a cavallo del primo layer
    
                U=zeros(size(indexTri,1),6); % punti di intersezione del piano con i triangoli
    
                for i=1:size(indexTri,1)
                    [~,temp]=plane_tri_inter([0 0 z(zThres)], [0 0 1], tri(indexTri(i),1:3), tri(indexTri(i),4:6), tri(indexTri(i),7:9));
                    U(i,:)=temp;
                end
    
                terzaColonna=trUp(:,3);
                sestaColonna=trUp(:,6);
                nonaColonna=trUp(:,9);
    
                R=cell(1,size(terzaColonna,1));
    
                for i=1:size(terzaColonna,1)
                    i3=find(terzaColonna>=z(zThres));
                    i6=find(sestaColonna>=z(zThres));
                    i9=find(nonaColonna>=z(zThres));
                end
    
                for i=1:size(i3,1)
                    R{1,i3(i)}=[trUp(i3(i),1) trUp(i3(i),2) trUp(i3(i),3)];
                end
    
                for i=1:size(i6,1)
                    R{2,i6(i)}=[trUp(i6(i),4) trUp(i6(i),5) trUp(i6(i),6)];
                end
    
                for i=1:size(i9,1)
                    R{3,i9(i)}=[trUp(i9(i),7) trUp(i9(i),8) trUp(i9(i),9)];
                end
    
                triMid=nan*ones(size(indexTri,1),9); % triangoli creati dall'intersezione con il layer
                quadMid=nan*ones(size(indexTri,1),12); % quadrilateri creati dall'intersezione con il layer
    
                for i=1:size(R,2)
                    count=0;
                    for j=1:3
                        if ~isempty(R{j,i})
                            count=count+1;
                        end
                    end
    
                    R{4,i}=count;
    
                    if count==1
                        triMid(i,:)=[U(i,1:3) U(i,4:6) cell2mat(R(1:3,i))];
                    elseif count==2
                        quad1=nan*ones(2,3);
                        for j=1:3
                            if ~isempty(R{j,i})
                                quad1(j,:)=R{j,i};
                            end
                        end
                        quad1(any(isnan(quad1), 2), :) = [];
                        quadMid(i,:)=[U(i,1:3) U(i,4:6) quad1(1,1) quad1(1,2) quad1(1,3) quad1(2,1) quad1(2,2) quad1(2,3)];
                    end
                end
    
                quadMid(any(isnan(quadMid), 2), :) = [];
                triMid(any(isnan(triMid), 2), :) = [];
    
                %ordino i punti dei quadrilateri
    
                for i=1:size(quadMid,1)
                    xyz =[quadMid(i,1:3);quadMid(i,4:6);quadMid(i,7:9);quadMid(i,10:12)];
                    xyzc = mean(xyz,1);
                    P = xyz - xyzc;
                    [~,~,V] = svd(P,0);
                    [~,is] = sort(atan2(P*V(:,1),P*V(:,2)));
                    xyz = xyz(is([1:end 1]),:);
                    quadMid(i,:)=[xyz(1,:) xyz(2,:) xyz(3,:) xyz(4,:)];
                end
    
                %sommo tutti i triangoli e quadrilateri originati dall'intersezione con
                %il primo layer
    
                pout1=polyshape([triMid(1,1),triMid(1,4),triMid(1,7)],[triMid(1,2),triMid(1,5),triMid(1,8)]);
    
                for i=2:size(triMid,1)
                    poly=polyshape([triMid(i,1),triMid(i,4),triMid(i,7)],[triMid(i,2),triMid(i,5),triMid(i,8)]);
                    poly=polybuffer(poly,0.001,'JointType','miter');
                    pout1=union(pout1,poly);
                end
    
                for i=1:size(quadMid,1)
                    poly=polyshape([quadMid(i,1),quadMid(i,4),quadMid(i,7),quadMid(i,10)],[quadMid(i,2),quadMid(i,5),quadMid(i,8),quadMid(i,11)]);
                    poly=polybuffer(poly,0.001,'JointType','miter');
                    pout1=union(pout1,poly);
                end
    
                % aggiungo anche tutti i triangoli con tutti i punti superiori al primo
                % layer, ottenendo cosi la proiezione della mesh 2D croppata
    
                pout=polyshape(TR_new.Points(TR_new.ConnectivityList(1,:),1),TR_new.Points(TR_new.ConnectivityList(1,:),2));
    
                for i = 2:size(TR_new.ConnectivityList,1)
                    p1 = polyshape(TR_new.Points(TR_new.ConnectivityList(i,:),1),TR_new.Points(TR_new.ConnectivityList(i,:),2));
                    pout=union(pout,p1);
                end
    
                pout=union(pout,pout1);
            else
                pout=polyshape(TR_new.Points(TR_new.ConnectivityList(1,:),1),TR_new.Points(TR_new.ConnectivityList(1,:),2));
                for i = 2:size(TR_new.ConnectivityList,1)
                    p1 = polyshape(TR_new.Points(TR_new.ConnectivityList(i,:),1),TR_new.Points(TR_new.ConnectivityList(i,:),2));
                    pout=union(pout,p1);
                end
            end
            pout=polybuffer(pout,0.00001,'JointType','miter');
    
            %pout=rmholes(pout); % rimuovo eventuali buchi nel perimetro
            polyFloat{1,n}=pout;
    
            % rimuovo i boundary con area infinitesima generati erronamente dal
            % calcolatore
    
            a2=polyFloat{1,n}.NumHoles+polyFloat{1,n}.NumRegions;
    
            polyArea=[];
            flag=0;
            j=1;
    
            while j<=a2
                polyArea=area(polyFloat{1,n},j);
                if abs(polyArea)<1
                    polyFloat{1,n}=rmboundary(polyFloat{1,n},j);
                    a2=polyFloat{1,n}.NumHoles+polyFloat{1,n}.NumRegions;
                else
                    j=j+1;
                end
            end
            polyFloat{1,n}=polybuffer(polyFloat{1,n},-off_nonplanar/1000,'JointType','miter');
        end
    
        %% Offset negativo sul perimetro più esterno per poterlo proiettare correttamente
    
        polyInnerFloat=cell(1,size(polyFloat,2));
    
        for i=1:size(polyFloat,2)
            if nNonPlanarPerimetersUp>=2
                for n=2:nNonPlanarPerimetersUp
                    polyFloat{n,i}=polybuffer(polyFloat{1,i},off_nonplanar*(n-1),'JointType','miter');
                end
            end
            polyInnerFloat{1,i}=polyFloat{end,i};
        end
    
        perimFloat=cell(1,size(polyFloat,2));
        innerFloat=cell(1,size(polyFloat,2));
    
        for i=1:size(polyFloat,2)
            w=1;
            polyInnerFloat{1,i}=polybuffer(polyInnerFloat{1,i},off_nonplanar_small,'JointType','miter');
            for n=1:size(polyFloat,1)
                polyFloat{n,i}=polybuffer(polyFloat{n,i},off_nonplanar_small,'JointType','miter');
                temp=polyFloat{n,i}.NumHoles+polyFloat{n,i}.NumRegions;
                for j=1:temp
                    [perimFloat{w,i}(:,1),perimFloat{w,i}(:,2)]=boundary(polyFloat{n,i},j);
                    perimFloat{w,i}(:,3)=(max(points(:,3))+5*layer_height+i)*ones(size(perimFloat{w,i},1),1);
                    w=w+1;
                end
            end
        end
    
        %% Perimetri interni (se ci sono)
    
        off_infill=-extrusion_width*(1-infill_overlap);
    
        for i=1:size(polyInnerFloat,2)
            w=1;
            polyInnerFloat{1,i}=polybuffer(polyInnerFloat{1,i},off_infill,'JointType','miter');
            temp=polyInnerFloat{1,i}.NumHoles+polyInnerFloat{1,i}.NumRegions;
            for j=1:temp
                [innerFloat{w,i}(:,1),innerFloat{w,i}(:,2)]=boundary(polyInnerFloat{1,i},j);
                innerFloat{w,i}(:,3)=(max(points(:,3))+layer_height+i)*ones(size(innerFloat{w,i},1),1);
                w=w+1;
            end
        end
    
        %% Infill
    
        infillFloat=lines1(innerFloat,dx,percNonPlanar,theta1NonPlanar,theta_relative_nonPlanar,points);
    
    
        %% Remesh
    
        theta_up=90; % non planar surface superiore (angolo soglia 90 gradi)
    
        F = faceNormal(TR); %normali alle face dei triangoli
    
        u=[0 0 1]; %normale all'asse z
        ThetaInDegrees=zeros(length(F),1); % calcolo gli angoli tra le normali e l'asse z
        k=1;
    
        for i=1:length(F)
            CosTheta = max(min(dot(u,F(i,:))/(norm(u)*norm(F(i,:))),1),-1);
            ThetaInDegrees(k) = real(acosd(CosTheta));
            k=k+1;
        end
    
        indexUP=find(ThetaInDegrees<theta_up); %scarto tutte le normali che formano un angolo maggiore di 90 con l'asse z
    
        tri_up=tri(indexUP,:);
        index_up=index_tri(indexUP,:);
    
        %% Resample
    
        resPerim=resample_cell(perimFloat,spatial_resolution);
        resInfill=resample_cell(infillFloat,spatial_resolution);
    
        %% Projection
    
        % con findtria trovo gli inidici dei triangoli in cui sono contenuti i
        % punti da proiettare
        % con trianglerayint calcolo la proiezione dei punti sui triangoli trovati
        % prima
    
        projectPerim=cell(size(resPerim,1),size(resPerim,2));
        in1=cell(1,size(resPerim,2));
        projectInfill=cell(size(resInfill,1),size(resInfill,2));
        in2=cell(1,size(resInfill,2));
        listPerim=cell(size(resInfill,1),size(resInfill,2));
    
        for w=1:size(resPerim,2)
            for i=1:size(resPerim,1)
                k=1;
                if ~isempty(resPerim{i,w})
                    [tp,ti] = findtria(points(:,1:2),index_up,resPerim{i,w}(:,1:2));
                    for j=1:size(resPerim{i,w},1)
                        IDperim=ti(tp(j,1):tp(j,2));
    
                        if ~isempty(IDperim)
                            [flag,~,~,~,inters] = TriangleRayIntersection(resPerim{i,w}(j,:), [0 0 -1], tri_up(IDperim,1:3), tri_up(IDperim,4:6), tri_up(IDperim,7:9),'planeType','one sided','border','inclusive');
                            intersIDperim=[inters IDperim];
                            intersIDperim= intersIDperim(~any(isnan(intersIDperim),2),:);
                            intersIDperim=sortrows(intersIDperim,3,'descend');
                            if ~isempty(intersIDperim)
                                projectPerim{i,w}(k,:)=intersIDperim(1,:);
                                k=k+1;
                            end
                        end
                    end
                end
    
                if sample_density=="Mesh" % rimuovo i punti ridondanti, ovvero ottengo 2 punti per ogni triangolo tranne per i punti piu esterni
                    if ~isempty(projectPerim{i,w})
                        for t=2:size(projectPerim{i,w},1)-1
                            if projectPerim{i,w}(t,4)==projectPerim{i,w}(t-1,4) && projectPerim{i,w}(t,4)==projectPerim{i,w}(t+1,4)
                                projectPerim{i,w}(t,1:3)=nan;
                            end
                        end
                    end
                end
    
                if ~isempty(projectPerim{i,w})
                    projectPerim{i,w}= projectPerim{i,w}(~any(isnan(projectPerim{i,w}),2),:);
                end
    
                projectPerim{i,w}=projectPerim{i,w}(:,1:3)-[0 0 layer_height*(w-1)];
    
                %check if all points are inside the mesh
                %if w~=1
                %in1{1,w} = intriangulation(points,index_tri,projectPerim{i,w});
                %             projectPerim{i,w}=[projectPerim{i,w}(in1{1,w},1) projectPerim{i,w}(in1{1,w},2) projectPerim{i,w}(in1{1,w},3)];
                %outPoints=in1{1,w}(in1{1,w}==0);
                %             check=size(outPoints,1)/size(in1{1,w},1);
                %end
            end
        end
    
        for w=1:size(resInfill,2)
            tempp=[];
            for i=1:size(resInfill,1)
                k=1;
                if ~isempty(resInfill{i,w})
                    [tp,ti] = findtria(points(:,1:2),index_up,resInfill{i,w}(:,1:2));
                    for j=1:size(resInfill{i,w},1)
                        IDinfill=ti(tp(j,1):tp(j,2));
                        if ~isempty(IDinfill)
                            [~,~,~,~,inters1] = TriangleRayIntersection(resInfill{i,w}(j,:), [0 0 -1], tri_up(IDinfill,1:3), tri_up(IDinfill,4:6), tri_up(IDinfill,7:9),'planeType','one sided','border','inclusive');
                            intersIDinfill=[inters1 IDinfill];
                            intersIDinfill= intersIDinfill(~any(isnan(intersIDinfill),2),:);
                            intersIDinfill=sortrows(intersIDinfill,3,'descend');
                            if ~isempty(intersIDinfill)
                                projectInfill{i,w}(k,:)=intersIDinfill(1,:);
                                checkk=intersect(infillFloat{i,w},resInfill{i,w}(j,:),'rows');
                                if ~isempty(checkk)
                                    tempp=[tempp;projectInfill{i,w}(k,:)];
                                end
                                k=k+1;
                            end
                        end
                    end
    
                    listPerim{i,w}=tempp;
    
                    if sample_density=="Mesh"
                        if ~isempty(projectInfill{i,w})
                            for t=2:size(projectInfill{i,w},1)-1
                                checkk2=intersect(listPerim{i,w},projectInfill{i,w}(t,:),'rows');
                                if projectInfill{i,w}(t,4)==projectInfill{i,w}(t-1,4) && projectInfill{i,w}(t,4)==projectInfill{i,w}(t+1,4) && isempty(checkk2)
                                    projectInfill{i,w}(t,1:3)=nan;
                                end
                            end
                        end
                    end
    
                    if ~isempty(projectInfill{i,w})
                        projectInfill{i,w}= projectInfill{i,w}(~any(isnan(projectInfill{i,w}),2),:);
                    end
    
                    projectInfill{i,w}=projectInfill{i,w}(:,1:3)-[0 0 layer_height*(w-1)];
    
                    %check if all points are inside the mesh
    
                    %if w~=1
                    %in2{1,w} = intriangulation(points,index_tri,projectInfill{i,w});
                    %                     projectInfill{i,w}=[projectInfill{i,w}(in2{1,w},1) projectInfill{i,w}(in2{1,w},2) projectInfill{i,w}(in2{1,w},3)];
                    %end
                end
            end
        end
    
    
        warning('off')
    
        projectPerim=flip(projectPerim,2);
        projectInfill=flip(projectInfill,2);
    
        %%
        clear polyFloat
    
        %% Non Planar Down
    
        shell=flip(shell,2); % passo alle surface non planari di bottom, faccio un flip dei layer precedentemente calcolati
        percTriDown=flip(percTriDown); %flip del vettore delle percentuali di triangoli associate ai vari layer
    
        % ruoto i punti di ogni layer
    
        for i=1:size(shell,2)
            for j=1:size(shell,1)
                if ~isempty(shell{j,i})
                    shell{j,i}(:,1)=shell{j,i}(:,1)-mediaX;
                    shell{j,i}(:,2)=shell{j,i}(:,2)-mediaY;
                    shell{j,i}(:,3)=shell{j,i}(:,3)-mediaZ;
    
                    shell{j,i}=(roty(angle)*shell{j,i}(:,:)')';
    
                    shell{j,i}(:,1)=shell{j,i}(:,1)+mediaX;
                    shell{j,i}(:,2)=shell{j,i}(:,2)+mediaY;
                    shell{j,i}(:,3)=shell{j,i}(:,3)+mediaZ;
    
                end
            end
        end
    
        layerz=[ones(length(z),1),ones(length(z),1),z'];
    
        layerz(:,1)=layerz(:,1)-mediaX;
        layerz(:,2)=layerz(:,2)-mediaY;
        layerz(:,3)=layerz(:,3)-mediaZ;
    
        layerz=roty(angle)*layerz';
        layerz=layerz';
    
        layerz(:,1)=layerz(:,1)+mediaX;
        layerz(:,2)=layerz(:,2)+mediaY;
        layerz(:,3)=layerz(:,3)+mediaZ;
    
    
        z_new=layerz(:,3);
        z_new=flip(z_new);
    
        points(:,1)=points(:,1)-mediaX;
        points(:,2)=points(:,2)-mediaY;
        points(:,3)=points(:,3)-mediaZ;
    
        points=roty(angle)*points';
        points=points';
    
        points(:,1)=points(:,1)+mediaX;
        points(:,2)=points(:,2)+mediaY;
        points(:,3)=points(:,3)+mediaZ;
    
        index_tri = object.ConnectivityList; % gli indici dei triangoli non cambiano perchè sono associati ai medesimi punti che sono stati ruotati
        tri = zeros(length(index_tri),9); % triangoli della nuova mesh
    
        for h=1:length(index_tri)
            k = 1;
            for j=1:3:size(tri,2)
                tri(h,j:j+2) = points(index_tri(h,k), :);
                k = k+1;
            end
        end
    
        TR=triangulation(index_tri,points); % nuova mesh ruotata
    
        polyShell=cell2poly(shell);
    
        for i=1:size(polyShell,2)
            if isempty(polyShell{1,i})
                polyShell{1,i}=polyshape(0,0);
            end
        end
    
    
        %%
    
        percThresUp=0;
        percThresDown=0;
    
        zFirstDown=find(percTriDown>=percThresDown,1);
        zEndDown=size(polyShell,2);
    
    
        %%
        zThres=1;
        roomDown=cell(nNonPlanarLayersDown,size(shell,2)); % spazio da lasciare per le non planar surface
    
        for n=1:nNonPlanarLayersDown
            roomDown{n,zEndDown-(n-1)}=polyShell{1,zEndDown-(n-1)};
            for i=zFirstDown:zEndDown-n
                roomDown{n,i}=subtract(polyShell{1,i},polyShell{1,i+1});
                for j=i+1:zEnd
                    inters=intersect(roomDown{n,i},polyShell{1,j});
                    if size(inters.Vertices)>1
                        roomDown{n,i}=subtract(roomDown{n,i},inters);
                    end
                end
                polyShell{1,i}=subtract(polyShell{1,i},roomDown{n,i});
            end
            polyShell{1,zEnd-(n-1)}=subtract(polyShell{1,zEnd-(n-1)},roomDown{n,zEnd-(n-1)});
        end
    
        shell1=cell(size(shell,1),size(shell,2));
    
        for i=1:size(polyShell,2)
            [xt,yt]=boundary(polyShell{1,i});
            shell1{1,i}=[xt yt z_new(i)*ones(size(boundary(polyShell{1,i}),1),1)];
        end
    
        shell=shell1;
    
    
        %% Crop della mesh dal layer selezionato in su
    
        % per ogni layer non planare faccio un crop della mesh e proietto la
        % silohuette in 2D per ottenere il perimetro non planare fluttuante
    
        polyFloat=cell(1,nNonPlanarLayersDown);
    
        for n=1:nNonPlanarLayersDown
    
            for i=1:size(roomDown,2)
                if ~isempty(roomDown{n, i})
                    if ~isempty(roomDown{n, i}.Vertices)
                        zThres=i+n-1; % primo layer
                        break
                    end
                end
            end
    
            if zThres<=1
                TR_new=TR;
                indexTri=[];
            else
    
                % costruisco una nuova mesh con solo i triangoli sopra il primo layer in
                % modo da calcolarmi il perimetro della surface non planare
    
                index_new=index_tri;
    
                ind1=find(points(:,3)<z_new(zThres));
                temp=[];
    
                for i=1:size(ind1,1)
                    [i1,~]=find(index_new==ind1(i));
                    temp=[temp;i1];
                end
    
                for i=1:size(temp,1)
                    index_new(temp(i),:)= [nan nan nan];
                end
    
                index_new(any(isnan(index_new), 2), :) = [];
    
                TR_new=triangulation(index_new,points); % mesh croppata con i triangoli INTERAMENTE sopra il primo layer
                indexTri=find(tri(:,3)<z_new(zThres) & tri(:,6)>z_new(zThres) | tri(:,3)>z_new(zThres) & tri(:,6)<z_new(zThres) | tri(:,6)<z_new(zThres) & tri(:,9)>z_new(zThres)|tri(:,6)>z_new(zThres) & tri(:,9)< z_new(zThres) |tri(:,3)>z_new(zThres) & tri(:,9)< z_new(zThres)|tri(:,3)<z_new(zThres) & tri(:,9)>z_new(zThres));
    
            end
    
            if ~isempty(indexTri)
    
                trDown=tri(indexTri,:); % triangoli che sono a cavallo del primo layer
    
                U=zeros(size(indexTri,1),6); % punti di intersezione del piano con i triangoli
    
                for i=1:size(indexTri,1)
                    [~,temp]=plane_tri_inter([0 0 z_new(zThres)], [0 0 1], tri(indexTri(i),1:3), tri(indexTri(i),4:6), tri(indexTri(i),7:9));
                    U(i,:)=temp;
                end
    
                terzaColonna=trDown(:,3);
                sestaColonna=trDown(:,6);
                nonaColonna=trDown(:,9);
    
                R=cell(1,size(terzaColonna,1));
    
                for i=1:size(terzaColonna,1)
                    i3=find(terzaColonna>=z_new(zThres));
                    i6=find(sestaColonna>=z_new(zThres));
                    i9=find(nonaColonna>=z_new(zThres));
                end
    
                for i=1:size(i3,1)
                    R{1,i3(i)}=[trDown(i3(i),1) trDown(i3(i),2) trDown(i3(i),3)];
                end
    
                for i=1:size(i6,1)
                    R{2,i6(i)}=[trDown(i6(i),4) trDown(i6(i),5) trDown(i6(i),6)];
                end
    
                for i=1:size(i9,1)
                    R{3,i9(i)}=[trDown(i9(i),7) trDown(i9(i),8) trDown(i9(i),9)];
                end
    
                triMid=nan*ones(size(indexTri,1),9); % triangoli creati dall'intersezione con il layer
                quadMid=nan*ones(size(indexTri,1),12); % quadrilateri creati dall'intersezione con il layer
    
                for i=1:size(R,2)
                    count=0;
                    for j=1:3
                        if ~isempty(R{j,i})
                            count=count+1;
                        end
                    end
    
                    R{4,i}=count;
    
                    if count==1
                        triMid(i,:)=[U(i,1:3) U(i,4:6) cell2mat(R(1:3,i))];
                    elseif count==2
                        quad1=nan*ones(2,3);
                        for j=1:3
                            if ~isempty(R{j,i})
                                quad1(j,:)=R{j,i};
                            end
                        end
                        quad1(any(isnan(quad1), 2), :) = [];
                        quadMid(i,:)=[U(i,1:3) U(i,4:6) quad1(1,1) quad1(1,2) quad1(1,3) quad1(2,1) quad1(2,2) quad1(2,3)];
                    end
                end
    
                quadMid(any(isnan(quadMid), 2), :) = [];
                triMid(any(isnan(triMid), 2), :) = [];
    
                %ordino i punti dei quadrilateri
    
                for i=1:size(quadMid,1)
                    xyz =[quadMid(i,1:3);quadMid(i,4:6);quadMid(i,7:9);quadMid(i,10:12)];
                    xyzc = mean(xyz,1);
                    P = xyz - xyzc;
                    [~,~,V] = svd(P,0);
                    [~,is] = sort(atan2(P*V(:,1),P*V(:,2)));
                    xyz = xyz(is([1:end 1]),:);
                    quadMid(i,:)=[xyz(1,:) xyz(2,:) xyz(3,:) xyz(4,:)];
                end
    
                %sommo tutti i triangoli e quadrilateri originati dall'intersezione con
                %il primo layer
    
                pout1=polyshape([triMid(1,1),triMid(1,4),triMid(1,7)],[triMid(1,2),triMid(1,5),triMid(1,8)]);
    
                for i=2:size(triMid,1)
                    poly=polyshape([triMid(i,1),triMid(i,4),triMid(i,7)],[triMid(i,2),triMid(i,5),triMid(i,8)]);
                    poly=polybuffer(poly,0.001,'JointType','miter');
                    pout1=union(pout1,poly);
                end
    
                for i=1:size(quadMid,1)
                    poly=polyshape([quadMid(i,1),quadMid(i,4),quadMid(i,7),quadMid(i,10)],[quadMid(i,2),quadMid(i,5),quadMid(i,8),quadMid(i,11)]);
                    poly=polybuffer(poly,0.001,'JointType','miter');
                    pout1=union(pout1,poly);
                end
    
                % aggiungo anche tutti i triangoli con tutti i punti superiori al primo
                % layer, ottenendo cosi la proiezione della mesh 2D croppata
    
                pout=polyshape(TR_new.Points(TR_new.ConnectivityList(1,:),1),TR_new.Points(TR_new.ConnectivityList(1,:),2));
    
                for i = 2:size(TR_new.ConnectivityList,1)
                    p1 = polyshape(TR_new.Points(TR_new.ConnectivityList(i,:),1),TR_new.Points(TR_new.ConnectivityList(i,:),2));
                    pout=union(pout,p1);
                end
    
                pout=union(pout,pout1);
            else
                pout=polyshape(TR_new.Points(TR_new.ConnectivityList(1,:),1),TR_new.Points(TR_new.ConnectivityList(1,:),2));
    
                for i = 2:size(TR_new.ConnectivityList,1)
                    p1 = polyshape(TR_new.Points(TR_new.ConnectivityList(i,:),1),TR_new.Points(TR_new.ConnectivityList(i,:),2));
                    pout=union(pout,p1);
                end
            end
    
            pout=polybuffer(pout,0.00001,'JointType','miter');
            %     pout=rmholes(pout); % rimuovo eventuali buchi nel perimetro
            polyFloat{1,n}=pout;
    
            % rimuovo i boundary con area infinitesima generati erronamente dal
            % calcolatore
    
            a2=polyFloat{1,n}.NumHoles+polyFloat{1,n}.NumRegions;
    
            polyArea=[];
            flag=0;
            j=1;
    
            while j<=a2
                polyArea=area(polyFloat{1,n},j);
                if abs(polyArea)<0.1
                    polyFloat{1,n}=rmboundary(polyFloat{1,n},j);
                    a2=polyFloat{1,n}.NumHoles+polyFloat{1,n}.NumRegions;
                else
                    j=j+1;
                end
            end
            polyFloat{1,n}=polybuffer(polyFloat{1,n},-off_nonplanar/1000,'JointType','miter');
        end
    
        %% Offset negativo sul perimetro più esterno per poterlo proiettare correttamente
    
        polyInnerFloatDown=cell(1,size(polyFloat,2));
    
        for i=1:size(polyFloat,2)
            if nNonPlanarPerimetersDown>=2
                for n=2:nNonPlanarPerimetersDown
                    polyFloat{n,i}=polybuffer(polyFloat{1,i},off_nonplanar*(n-1),'JointType','miter');
                end
            end
            polyInnerFloatDown{1,i}=polyFloat{end,i};
        end
    
        perimFloatDown=cell(1,size(polyFloat,2));
        innerFloatDown=cell(1,size(polyFloat,2));
    
        for i=1:size(polyFloat,2)
            w=1;
            polyInnerFloatDown{1,i}=polybuffer(polyInnerFloatDown{1,i},off_nonplanar_small,'JointType','miter');
            for n=1:size(polyFloat,1)
                polyFloat{n,i}=polybuffer(polyFloat{n,i},off_nonplanar_small,'JointType','miter');
                temp=polyFloat{n,i}.NumHoles+polyFloat{n,i}.NumRegions;
                for j=1:temp
                    [perimFloatDown{w,i}(:,1),perimFloatDown{w,i}(:,2)]=boundary(polyFloat{n,i},j);
                    perimFloatDown{w,i}(:,3)=(max(points(:,3))+layer_height+i)*ones(size(perimFloatDown{w,i},1),1);
                    w=w+1;
                end
            end
        end
        %% Perimetri interni (se ci sono)
    
        off_infill=-extrusion_width*(1-infill_overlap);
    
        for i=1:size(polyInnerFloatDown,2)
            w=1;
            polyInnerFloatDown{1,i}=polybuffer(polyInnerFloatDown{1,i},off_infill,'JointType','miter');
            temp=polyInnerFloatDown{1,i}.NumHoles+polyInnerFloatDown{1,i}.NumRegions;
            for j=1:temp
                [innerFloatDown{w,i}(:,1),innerFloatDown{w,i}(:,2)]=boundary(polyInnerFloatDown{1,i},j);
                innerFloatDown{w,i}(:,3)=(max(points(:,3))+layer_height+i)*ones(size(innerFloatDown{w,i},1),1);
                w=w+1;
            end
        end
    
        %% Infill
    
        infillFloatDown=lines1(innerFloatDown,dx,percNonPlanar,-theta1NonPlanar,-theta_relative_nonPlanar,points);
    
    
        %% Remesh
    
        theta_up=90; % non planar surface superiore (angolo soglia 90 gradi)
    
        F = faceNormal(TR); %normali alle face dei triangoli
    
        u=[0 0 1]; %normale all'asse z
        ThetaInDegrees=zeros(length(F),1); % calcolo gli angoli tra le normali e l'asse z
        k=1;
    
        for i=1:length(F)
            CosTheta = max(min(dot(u,F(i,:))/(norm(u)*norm(F(i,:))),1),-1);
            ThetaInDegrees(k) = real(acosd(CosTheta));
            k=k+1;
        end
    
        indexDOWN=find(ThetaInDegrees<theta_up); %scarto tutte le normali che formano un angolo maggiore di 90 con l'asse z
    
        tri_down=tri(indexDOWN,:);
        index_down=index_tri(indexDOWN,:);
    
        %% Resample
    
        resPerimDown=resample_cell(perimFloatDown,spatial_resolution);
        resInfillDown=resample_cell(infillFloatDown,spatial_resolution);
    
    
        %% Projection
        warning('on');
    
        % con findtria trovo gli inidici dei triangoli in cui sono contenuti i
        % punti da proiettare
        % con trianglerayint calcolo la proiezione dei punti sui triangoli trovati
        % prima
    
        projectPerimDown=cell(size(resPerimDown,1),size(resPerimDown,2));
        in1=cell(1,size(resPerimDown,2));
        projectInfillDown=cell(size(resInfillDown,1),size(resInfillDown,2));
        in2=cell(1,size(resInfillDown,2));
        listPerimDown=cell(size(resInfillDown,1),size(resInfillDown,2));
    
        for w=1:size(resPerimDown,2)
            for i=1:size(resPerimDown,1)
                k=1;
                if ~isempty(resPerimDown{i,w})
                    [tp,ti] = findtria(points(:,1:2),index_down,resPerimDown{i,w}(:,1:2));
                    for j=1:size(resPerimDown{i,w},1)
                        IDperim=ti(tp(j,1):tp(j,2));
    
                        if ~isempty(IDperim)
                            [~,~,~,~,inters] = TriangleRayIntersection(resPerimDown{i,w}(j,:), [0 0 -1], tri_down(IDperim,1:3), tri_down(IDperim,4:6), tri_down(IDperim,7:9),'planeType','one sided','border','inclusive');
                            intersIDperim=[inters IDperim];
                            intersIDperim= intersIDperim(~any(isnan(intersIDperim),2),:);
                            intersIDperim=sortrows(intersIDperim,3,'descend');
                            if ~isempty(intersIDperim)
                                projectPerimDown{i,w}(k,:)=intersIDperim(1,:);
                                k=k+1;
                            end
                        end
                    end
                end
    
                if sample_density=="Mesh" % rimuovo i punti ridondanti, ovvero ottengo 2 punti per ogni triangolo tranne per i punti piu esterni
                    if ~isempty(projectPerimDown{i,w})
                        for t=2:size(projectPerimDown{i,w},1)-1
                            if projectPerimDown{i,w}(t,4)==projectPerimDown{i,w}(t-1,4) && projectPerimDown{i,w}(t,4)==projectPerimDown{i,w}(t+1,4)
                                projectPerimDown{i,w}(t,1:3)=nan;
                            end
                        end
                    end
                end
    
                if ~isempty(projectPerimDown{i,w})
                    projectPerimDown{i,w}= projectPerimDown{i,w}(~any(isnan(projectPerimDown{i,w}),2),:);
                end
    
                projectPerimDown{i,w}=projectPerimDown{i,w}(:,1:3)-[0 0 layer_height*(w-1)];
    
                % check if all points are inside the mesh
    
                %if w~=1
                %in1{1,w} = intriangulation(points,index_tri,projectPerimDown{i,w});
                %             projectPerimDown{i,w}=[projectPerimDown{i,w}(in1{1,w},1) projectPerimDown{i,w}(in1{1,w},2) projectPerimDown{i,w}(in1{1,w},3)];
                %outPoints=in1{1,w}(in1{1,w}==0);
                %             check=size(outPoints,1)/size(in1{1,w},1);
                %end
            end
        end
    
        for w=1:size(resInfillDown,2)
            tempp=[];
            for i=1:size(resInfillDown,1)
                k=1;
                if ~isempty(resInfillDown{i,w})
                    [tp,ti] = findtria(points(:,1:2),index_down,resInfillDown{i,w}(:,1:2));
                    for j=1:size(resInfillDown{i,w},1)
                        IDinfill=ti(tp(j,1):tp(j,2));
                        if ~isempty(IDinfill)
                            [~,~,~,~,inters1] = TriangleRayIntersection(resInfillDown{i,w}(j,:), [0 0 -1], tri_down(IDinfill,1:3), tri_down(IDinfill,4:6), tri_down(IDinfill,7:9),'planeType','one sided','border','inclusive');
                            intersIDinfill=[inters1 IDinfill];
                            intersIDinfill= intersIDinfill(~any(isnan(intersIDinfill),2),:);
                            intersIDinfill=sortrows(intersIDinfill,3,'descend');
                            if ~isempty(intersIDinfill)
                                projectInfillDown{i,w}(k,:)=intersIDinfill(1,:);
                                checkk=intersect(infillFloatDown{i,w},resInfillDown{i,w}(j,:),'rows');
                                if ~isempty(checkk)
                                    tempp=[tempp;projectInfillDown{i,w}(k,:)];
                                end
                                k=k+1;
                            end
                        end
                    end
    
                    listPerimDown{i,w}=tempp;
    
                    if sample_density=="Mesh"
                        if ~isempty(projectInfillDown{i,w})
                            for t=2:size(projectInfillDown{i,w},1)-1
                                checkk2=intersect(listPerimDown{i,w},projectInfillDown{i,w}(t,:),'rows');
                                if projectInfillDown{i,w}(t,4)==projectInfillDown{i,w}(t-1,4) && projectInfillDown{i,w}(t,4)==projectInfillDown{i,w}(t+1,4) && isempty(checkk2)
                                    projectInfillDown{i,w}(t,1:3)=nan;
                                end
                            end
                        end
                    end
    
                    if ~isempty(projectInfillDown{i,w})
                        projectInfillDown{i,w}= projectInfillDown{i,w}(~any(isnan(projectInfillDown{i,w}),2),:);
                    end
    
                    projectInfillDown{i,w}=projectInfillDown{i,w}(:,1:3)-[0 0 layer_height*(w-1)];
    
                    %         check if all points are inside the mesh
    
                    %if w~=1
                    %in2{1,w} = intriangulation(points,index_tri,projectInfillDown{i,w});
                    %                     projectInfillDown{i,w}=[projectInfillDown{i,w}(in2{1,w},1) projectInfillDown{i,w}(in2{1,w},2) projectInfillDown{i,w}(in2{1,w},3)];
                    %end
                end
            end
        end
    
        warning('off')
    
        polyShell=cell(size(shell,1),size(shell,2));
        tempShell=cell(size(shell,1),size(shell,2));
    
        for i=1:size(shell,2)
            if ~isempty(shell{1,i})
                polyShell{1,i}=polyshape(shell{1,i}(:,1),shell{1,i}(:,2));
                a2=polyShell{1,i}.NumHoles+polyShell{1,i}.NumRegions;
    
                polyArea=[];
                flag=0;
                j=1;
    
                while flag==0
                    if a2>0
                        polyArea=area(polyShell{1,i},j);
                        if abs(polyArea)<0.1
                            polyShell{1,i}=rmboundary(polyShell{1,i},j);
                            a2=polyShell{1,i}.NumHoles+polyShell{1,i}.NumRegions;
                            j=1;
                        else
                            j=j+1;
                        end
    
                        if j==a2+1
                            flag=1;
                        end
                    else
                        break
                    end
                end
    
                [tempShell{1,i}(:,1),tempShell{1,i}(:,2)]=boundary(polyShell{1,i});
                tempShell{1,i}(:,3)=shell{1,i}(1,3)*ones(size(tempShell{1,i},1),1);
            end
        end
    
        shell=tempShell;
    
        shell=flip(shell,2);
    
        for i=1:size(shell,2)
            for j=1:size(shell,1)
                if ~isempty(shell{j,i})
    
                    shell{j,i}(:,1)=shell{j,i}(:,1)-mediaX;
                    shell{j,i}(:,2)=shell{j,i}(:,2)-mediaY;
                    shell{j,i}(:,3)=shell{j,i}(:,3)-mediaZ;
    
                    shell{j,i}=(roty(angle)*shell{j,i}(:,:)')';
    
                    shell{j,i}(:,1)=shell{j,i}(:,1)+mediaX;
                    shell{j,i}(:,2)=shell{j,i}(:,2)+mediaY;
                    shell{j,i}(:,3)=shell{j,i}(:,3)+mediaZ;
                end
            end
        end
    
        points(:,1)=points(:,1)-mediaX;
        points(:,2)=points(:,2)-mediaY;
        points(:,3)=points(:,3)-mediaZ;
    
        points=roty(angle)*points';
        points=points';
    
        points(:,1)=points(:,1)+mediaX;
        points(:,2)=points(:,2)+mediaY;
        points(:,3)=points(:,3)+mediaZ;
    
        polyShell=cell2poly(shell);
    end
    
    %% Planar
    
    polyInner=cell(1,size(polyShell,2));
    perim=cell(1,size(polyShell,2));
    perimInner=cell(1,size(polyShell,2));
    
    for i=1:size(polyShell,2)
        if ~isempty(polyShell{1,i})
            polyShell{1,i}=polybuffer(polyShell{1,i},off_planar_small,'JointType','miter');
            if nPlanarPerimeters>=2
                for n=2:nPlanarPerimeters
                    polyShell{n,i}=polybuffer(polyShell{1,i},off_planar*(n-1),'JointType','miter');
                end
            end
            polyInner{1,i}=polyShell{end,i};
        end
    end
    
    for i=1:size(polyShell,2)
        w=1;
        for n=1:size(polyShell,1)
            if ~isempty(polyShell{1,i})
                temp=polyShell{n,i}.NumHoles+polyShell{n,i}.NumRegions;
                for j=1:temp
                    [perim{w,i}(:,1),perim{w,i}(:,2)]=boundary(polyShell{n,i},j);
                    perim{w,i}(:,3)=shell{1,i}(1,3)*ones(size(perim{w,i},1),1);
                    w=w+1;
                end
            end
        end
    end
    
    for i=1:size(polyInner,2)
        w=1;
        for n=1:size(polyInner,1)
            if ~isempty(polyInner{1,i})
                temp=polyInner{n,i}.NumHoles+polyInner{n,i}.NumRegions;
                for j=1:temp
                    [perimInner{w,i}(:,1),perimInner{w,i}(:,2)]=boundary(polyInner{n,i},j);
                    perimInner{w,i}(:,3)=shell{1,i}(1,3)*ones(size(perimInner{w,i},1),1);
                    w=w+1;
                end
            end
        end
    end
    
    %% Infill planare
    
    infill=lines1(perimInner,dx,percPlanar,theta1Planar,theta_relative_planar,points);
    
    for i=1:size(infill,2)
        for j=1:size(infill,1)
            if ~isempty(infill{j,i})
                infill{j,i}(any(isnan(infill{j,i}), 2), :) = [];
            end
        end
    end
    
    %% Concateno
    
    if non_planar==1
        for i=1:size(projectPerimDown,2)
            for j=1:size(projectPerimDown,1)
                C1=[];
                if ~isempty(projectPerimDown{j,i})
                    C1(:,1)=projectPerimDown{j,i}(:,1)-mediaX;
                    C1(:,2)=projectPerimDown{j,i}(:,2)-mediaY;
                    C1(:,3)=projectPerimDown{j,i}(:,3)-mediaZ;
    
                    C1=(roty(angle)*C1')';
    
                    C1(:,1)=C1(:,1)+mediaX;
                    C1(:,2)=C1(:,2)+mediaY;
                    C1(:,3)=C1(:,3)+mediaZ;
    
                    projectPerimDown{j,i}=C1;
                end
            end
    
            for j=1:size(projectInfillDown,1)
                C1=[];
                if ~isempty(projectInfillDown{j,i})
                    C1(:,1)=projectInfillDown{j,i}(:,1)-mediaX;
                    C1(:,2)=projectInfillDown{j,i}(:,2)-mediaY;
                    C1(:,3)=projectInfillDown{j,i}(:,3)-mediaZ;
    
                    C1=(roty(angle)*C1')';
    
                    C1(:,1)=C1(:,1)+mediaX;
                    C1(:,2)=C1(:,2)+mediaY;
                    C1(:,3)=C1(:,3)+mediaZ;
    
                    projectInfillDown{j,i}=C1;
                end
            end
        end
    end
    
    s1=[size(perim,1);size(infill,1)];
    smax=max(s1);
    
    A=cell(smax,size(perim,2)*2);
    
    for i=1:size(perim,2)
        A(1:size(perim,1),i+(i-1))=perim(:,i);
        A(1:size(infill,1),i+1+(i-1))=infill(:,i);
    end
    
    keep = any(~cellfun('isempty',A), 1);
    A = A(:,keep);
    
    if non_planar==1
        s1=[size(projectInfill,1);size(projectPerim,1)];
        smax=max(s1);
    
        B=cell(smax,size(projectPerim,2)*2);
    
        for i=1:size(projectPerim,2)
            B(1:size(projectPerim,1),i+(i-1))=projectPerim(:,i);
            B(1:size(projectInfill,1),i+1+(i-1))=projectInfill(:,i);
        end
    
        keep = any(~cellfun('isempty',B), 1);
        B = B(:,keep);
    
        s1=[size(projectInfillDown,1);size(projectPerimDown,1)];
        smax=max(s1);
    
        C=cell(smax,size(projectPerimDown,2)*2);
    
        for i=1:size(projectPerimDown,2)
            C(1:size(projectPerimDown,1),i+(i-1))=projectPerimDown(:,i);
            C(1:size(projectInfillDown,1),i+1+(i-1))=projectInfillDown(:,i);
        end
    
        keep = any(~cellfun('isempty',C), 1);
        C = C(:,keep);
    
        rC=size(C,1);
        cC=size(C,2);
    
        rA=size(A,1);
        cA=size(A,2);
    
        rB=size(B,1);
        cB=size(B,2);
    
        rr=[rA;rB;rC];
        rMax=max(rr);
        cMax=cA+cB+cC;
    
        TOT=cell(rMax,cMax);
    
        TOT(1:rC,1:cC)=C;
        TOT(1:rA,cC+1:cC+cA)=A;
        TOT(1:rB,cC+cA+1:cC+cA+cB)=B;
    end
    
    new_A={};
    A_true=[];
    
    new_B={};
    B_true=[];
    
    new_C={};
    C_true=[];
    
    for i=1:size(A,2)
        for j=1:size(A,1)
            if ~isempty(A{j,i})
                new_A{j,i}=[A{j,i}(1,1:2) zRet; ...
                    A{j,i}(1,:); ...
                    nan nan nan;A{j,i}(2:end,:);...
                    nan nan 0;
                    A{j,i}(end,1:2) zRet];
            end
    
            if ~isempty(A{j,i})
                A_true=[A_true;A{j,i}; nan nan nan];
            end
        end
    end
    
    if non_planar==1
        for i=1:size(B,2)
            for j=1:size(B,1)
                if ~isempty(B{j,i})
                    new_B{j,i}=[B{j,i}(1,1:2) zRet; ...
                        B{j,i}(1,:); ...
                        nan nan nan;B{j,i}(2:end,:);...
                        nan nan 0;
                        B{j,i}(end,1:2) zRet];
                end
    
                if ~isempty(B{j,i})
                    B_true=[B_true;B{j,i}; nan nan nan];
                end
            end
        end
    
        new_C={};
        C_true=[];
    
        for i=1:size(C,2)
            for j=1:size(C,1)
                if ~isempty(C{j,i})
                    new_C{j,i}=[C{j,i}(1,1:2) zRet; ...
                        C{j,i}(1,:); ...
                        nan nan nan;C{j,i}(2:end,:);...
                        nan nan 0;
                        C{j,i}(end,1:2) zRet];
                end
    
                if ~isempty(C{j,i})
                    C_true=[C_true;C{j,i}; nan nan nan];
                end
            end
        end
    
        for i=1:size(TOT,2)
            for j=1:size(TOT,1)
                if ~isempty(TOT{j,i})
                    new_TOT{j,i}=[TOT{j,i}(1,1:2) zRet; ...
                        TOT{j,i}(1,:); ...
                        nan nan nan;TOT{j,i}(2:end,:);...
                        nan nan 0;
                        TOT{j,i}(end,1:2) zRet];
                end
            end
        end
    end
    
    TOT_nan=[];
    TOT_number=[];
    
    BOTTOM_nan=[];
    BOTTOM_number=[];
    
    TOP_nan=[];
    TOP_number=[];
    
    PLANAR_nan=[];
    
    for i=1:size(new_A,2)
        PLANAR_nan=[PLANAR_nan;cell2mat(new_A(:,i))];
    end
    
    PLANAR_nan=round(PLANAR_nan,5);
    
    i=1;
    while i<length(PLANAR_nan)-1
        if PLANAR_nan(i,1)==PLANAR_nan(i+1,1) && PLANAR_nan(i,2)==PLANAR_nan(i+1,2) && PLANAR_nan(i,3)==PLANAR_nan(i+1,3)
            PLANAR_nan(i,:)=[];
        end
        i=i+1;
    end
    
    PLANAR_number=PLANAR_nan;
    PLANAR_number(any(isnan(PLANAR_number), 2), :) = [];
  
    
    if non_planar==1
    
        for i=1:size(new_B,2)
            TOP_nan=[TOP_nan;cell2mat(new_B(:,i))];
        end
    
        TOP_nan=round(TOP_nan,5);
    
        i=1;
        while i<length(TOP_nan)-1
            if TOP_nan(i,1)==TOP_nan(i+1,1) && TOP_nan(i,2)==TOP_nan(i+1,2) && TOP_nan(i,3)==TOP_nan(i+1,3)
                TOP_nan(i,:)=[];
            end
            i=i+1;
        end
    
        TOP_number=TOP_nan;
        TOP_number(any(isnan(TOP_number), 2), :) = [];
    
        for i=1:size(new_C,2)
            BOTTOM_nan=[BOTTOM_nan;cell2mat(new_C(:,i))];
        end
    
    
        BOTTOM_nan=round(BOTTOM_nan,5);
    
        i=1;
        while i<length(BOTTOM_nan)-1
            if BOTTOM_nan(i,1)==BOTTOM_nan(i+1,1) && BOTTOM_nan(i,2)==BOTTOM_nan(i+1,2) && BOTTOM_nan(i,3)==BOTTOM_nan(i+1,3)
                BOTTOM_nan(i,:)=[];
            end
            i=i+1;
        end
    
        BOTTOM_number=BOTTOM_nan;
        BOTTOM_number(any(isnan(BOTTOM_number), 2), :) = [];
        
        for i=1:size(new_TOT,2)
            TOT_nan=[TOT_nan;cell2mat(new_TOT(:,i))];
        end
    
        TOT_nan=round(TOT_nan,5);
    
        i=1;
        while i<length(TOT_nan)-1
            if TOT_nan(i,1)==TOT_nan(i+1,1) && TOT_nan(i,2)==TOT_nan(i+1,2) && TOT_nan(i,3)==TOT_nan(i+1,3)
                TOT_nan(i,:)=[];
            end
            i=i+1;
        end
    
        TOT_number=TOT_nan;
        TOT_number(any(isnan(TOT_number), 2), :) = [];
    end
end