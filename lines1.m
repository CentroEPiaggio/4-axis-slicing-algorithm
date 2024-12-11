function infill=lines1(perimInner,dx,percPlanar,theta1Planar,theta_relative_planar,points)

theta2=theta1Planar+theta_relative_planar;

[center,radius] = minboundcircle(points(:,1),points(:,2)); 
%circonferenza centrata nel valor medio di points e che contiene l'intera loro proiezione
infill=cell(1,size(perimInner,2));

for n=1:size(perimInner,2)
    if ~isempty(perimInner{1,n}) % se la prima cella di ogni piano di infill è vuota passo al piano successivo
        w=1;
        check=(-1)^n; %piano con indice pari o dispari  
        ind=1;
        temp=0;
        x=[];
        
        indexFirst=nan*ones(size(perimInner,1),1); %matrici degli indici dei primi e ultimi elementi di ogni cella
        indexLast=nan*ones(size(perimInner,1),1);
        
        % parametri delle rette che intersecano la geometria
        
        interval=[center(1)-radius:dx:center(1)+radius;center(1)-radius:dx:center(1)+radius];
        
        max1=floor(center(2)-radius)-1;
        max2=ceil(center(2)+radius)+1;
        rette=nan*ones(2*size(interval,2),2);
        k=1;

        if size(interval,2)>2
            for i=1:2:length(interval)-1
                rette(k,:)=[interval(1,i) max1];
                k=k+1;
                rette(k,:)=[interval(2,i) max2];
                k=k+1;
                rette(k,:)=[interval(2,i+1) max2];
                k=k+1;
                rette(k,:)=[interval(1,i+1) max1];
                k=k+1;
            end 
        else
             rette(k,:)=[interval(1,1) max1];
             k=k+1;
             rette(k,:)=[interval(2,1) max2];  
        end

        rette=round(rette,8);

        rotMatrix1= [cosd(theta1Planar) -sind(theta1Planar); sind(theta1Planar) cosd(theta1Planar)];
        rotMatrixInv1=[cosd(-theta1Planar) -sind(-theta1Planar); sind(-theta1Planar) cosd(-theta1Planar)];
        
        rotMatrix2= [cosd(theta2) -sind(theta2); sind(theta2) cosd(theta2)];
        rotMatrixInv2=[cosd(-theta2) -sind(-theta2); sind(-theta2) cosd(-theta2)];
        
        for j=1:size(perimInner,1)
            if ~isempty(perimInner{j,n}) && check<0
                shape=perimInner{j,n}(:,1:2); 
                shape(:,1)=shape(:,1)-center(1);
                shape(:,2)=shape(:,2)-center(2);
                shape=rotMatrix1*shape';
                shape=shape';
                shape(:,1)=shape(:,1)+center(1);
                shape(:,2)=shape(:,2)+center(2);
                shape=round(shape,8);
            elseif ~isempty(perimInner{j,n}) && check>0
                shape=perimInner{j,n}(:,1:2); 
                shape(:,1)=shape(:,1)-center(1);
                shape(:,2)=shape(:,2)-center(2);
                shape=rotMatrix2*shape';
                shape=shape';
                shape(:,1)=shape(:,1)+center(1);
                shape(:,2)=shape(:,2)+center(2);
                shape=round(shape,8);
            else
                break
            end

            rette(any(isnan(rette), 2), :) = [];
            X=poly2poly(shape',rette'); %intersezioni delle rette con ogni cella di ogni piano
            X=X';
            X=round(X,8);
            for i=1:size(X,1)-1
                if X(i,1)==X(i+1,1) &&  X(i,2)==X(i+1,2) %elimino i doppioni
                    X(i,:)=nan;
                end
            end
        
            X(any(isnan(X), 2), :) = [];
            
            if size(X,1)>1 & X(1,:)==X(end,:) %elimino i doppioni agli estremi dell'array
                    X=X(1:end-1,:);
            end
        
            for i=1:size(X,1)
                a1=find(X(:,1)==X(i,1));
                count=size(a1,1);
                if (-1)^count~=1
                    rette(:,1)=rette(:,1)+dx/10;
                    rette=round(rette,8);
                    X=poly2poly(shape',rette'); %intersezioni delle rette con ogni cella di ogni piano
                    X=X';
                    X=round(X,8);
                    for q=1:size(X,1)-1
                        if X(q,1)==X(q+1,1) &&  X(q,2)==X(q+1,2) %elimino i doppioni
                            X(q,:)=nan;
                        end
                    end
                    X(any(isnan(X), 2), :) = [];
                    if size(X,1)>1 & X(1,:)==X(end,:) %elimino i doppioni agli estremi dell'array
                        X=X(1:end-1,:);
                    end
                end
            end
            X(any(isnan(X), 2), :) = [];
            X=round(X,8);

            for i=1:size(X,1)
                if ~isnan(X(i,1))
                    x=[x;X(i,:)];
                end
            end  

            x=round(x,8);
            
            if size(X,1)>1
                indexFirst(ind,1)=temp+1;
                indexLast(ind,1)=size(x,1);
                temp=size(x,1);
                ind=ind+1;
            end
        end

        if size(x,1)>1
            
            allRette=unique(rette(:,1),'rows');
            numTot=size(allRette,1);
            nx=ceil(numTot*percPlanar);
            trueRette=round(numTot/nx);
            xi=[];

            for i=1:trueRette:size(allRette,1)
                xi=[xi;allRette(i)];
            end

            xi=round(xi,8);

            tot=[];
            sharedvalue=intersect(x(:,1),xi);
            sharedvalue=round(sharedvalue,8);
    
            %cerco tra i punti rimasti del perimetro qualcuno da
            %cui ripartire
    
            if ~isempty(sharedvalue)
                idx=find(x(:,1)==sharedvalue(1));
                tot=[tot;x(idx(1),:)];
                x(idx(1),:)=[nan nan];
            else
                continue
                %passo al layer successivo
            end
    
            idx=find(x(:,1)==tot(end,1)); %cerco il nuovo indice su cui spostarmi con la stessa ascissa
            val=size(idx,1);
            flag=0;
  
            for i=1:size(x,1)
                if val==1 %due sole intersezioni
                    tot=[tot; x(idx,:)];
                    x(idx,:)=[nan nan];
                else % se ci sono piu di due intersezioni 
                    t=[idx;nan];
                    t(:,2)=[x(idx(:,1),2);tot(end,2)];
                    val1=size(t,1);
                    t=sortrows(t,2,'descend');
                    f=find(isnan(t(:,1)));
                    if f==1
                        tot=[tot;x(t(2,1),:)];
                        idx=(t(2,1));
                        x(t(2,1),:)=[nan nan];
                    elseif f==val1
                        t=flip(t(:,1));
                        tot=[tot;x(t(2),:)];
                        idx=t(2,1);
                        x(t(2,1),:)=[nan nan];
                    else 
                        %il campione successivo lo aggiungo sempre dal lato dove ho un
                        %numero dispari di campioni con la stessa ascissa
                        up=f-1;
                        down=val1-f;
                        valu=(-1)^up;
                        vald=(-1)^down;
            
                        if vald<0 
                            % ho un numero dispari di campioni sotto
                             tot=[tot;x(t(f+1,1),:)];
                             x(t(f+1,1),:)=[nan nan];
                             idx=t(f+1,1);
                        elseif valu<0 
                            % ho un numero dispari di campioni sopra
                             t=flip(t,1);
                             f=find(isnan(t(:,1)));
                             tot=[tot;x(t(f+1,1),:)];
                             x(t(f+1,1),:)=[nan nan];
                             idx=t(f+1,1);
                        end
                    end    
                end
            
                if ~isempty(idx)

                %controllo se il perimetro va verso dx o sx per
                %capire se scorrere la x verso l'alto o verso il basso
                        
                    if isempty(find(idx==indexLast, 1))
                        a=tot(end,1)-x(idx+1,1);
                    elseif ~isempty(find(idx==indexLast(:,1), 1))
                        ii= find(idx==indexLast(:,1), 1);
                        a=tot(end,1)-x(indexFirst(ii),1);
                    end

                    if isempty(find(idx==indexFirst, 1))
                        b=tot(end,1)-x(idx-1,1);
                    elseif ~isempty(find(idx==indexFirst(:,1), 1))
                        jj=find(idx==indexFirst(:,1), 1);
                        b=tot(end,1)-x(indexLast(jj),1);
                    end

                    if b<0 
                    % scorro verso l'alto

                        if isempty(find(idx==indexFirst, 1)) && ~isnan(x(idx-1,1))
                            tot=[tot; x(idx-1,:)];
                            x(idx-1,:)=[nan nan];
                            idx=idx-1;
                        elseif ~isempty(find(idx==indexFirst(:,1), 1))
                            jj=find(idx==indexFirst(:,1), 1);
                            if ~isnan(x(indexLast(jj),1))
                                tot=[tot; x(indexLast(jj),:)];
                                x(indexLast(jj),:)=[nan nan];
                                idx=indexLast(jj);
                            end
                        else
                            if check>0
                                tot(:,1)=tot(:,1)-center(1);
                                tot(:,2)=tot(:,2)-center(2);
                                tot=rotMatrixInv2*tot'; 
                                tot=tot';
                                tot(:,1)=tot(:,1)+center(1);
                                tot(:,2)=tot(:,2)+center(2);
                            else
                                tot(:,1)=tot(:,1)-center(1);
                                tot(:,2)=tot(:,2)-center(2);
                                tot=rotMatrixInv1*tot'; 
                                tot=tot';
                                tot(:,1)=tot(:,1)+center(1);
                                tot(:,2)=tot(:,2)+center(2);
                            end
                              
                            zVector=perimInner{1,n}(1,3)*ones(size(tot,1),1);
                            infill{w,n}=[tot zVector];
                            tot=[];
                            w=w+1;
                            
                            sharedvalue=intersect(x(:,1),xi);
                            %cerco tra i punti rimasti del perimetro qualcuno da
                            %cui ripartire

                            if ~isempty(sharedvalue)
                                idx=find(x(:,1)==sharedvalue(1));
                                tot=[tot;x(idx(1),:)];
                                x(idx(1),:)=[nan nan];
                                idx=idx(2:end);
                                val=size(idx,1);
                            else
                                flag=1;                        
                                % esco da tutti i cicli for e finisce l'algoritmo
                            end
                                
                            break
                        end
                        for j=1:size(x,1)
                            if ~isempty(find(xi==tot(end,1), 1))
                                %ho trovato un punto in cui spostarmi in verticale
                                idx=find(x(:,1)==tot(end,1));
                                val=size(idx,1);
                                break                   
                            elseif isempty(find(idx==indexFirst, 1)) && ~isnan(x(idx-1,1))
                                tot=[tot;x(idx-1,:)];
                                x(idx-1,:)=[nan nan];
                                idx = idx-1; 
                                val=size(idx,1);
                            elseif  ~isempty(find(idx==indexFirst(:,1), 1)) 
                                jj=find(idx==indexFirst(:,1), 1);
                                if ~isnan(x(indexLast(jj),1))
                                    tot=[tot;x(indexLast(jj),:)];
                                    x(indexLast(jj),:)=[nan nan];
                                    idx=indexLast(jj);
                                    val=size(idx,1);
                                end
                            else
                                if check>0
                                   
                                       tot(:,1)=tot(:,1)-center(1);
                                       tot(:,2)=tot(:,2)-center(2);
                                       tot=rotMatrixInv2*tot'; 
                                       tot=tot';
                                       tot(:,1)=tot(:,1)+center(1);
                                       tot(:,2)=tot(:,2)+center(2);
                                  else
                                       tot(:,1)=tot(:,1)-center(1);
                                       tot(:,2)=tot(:,2)-center(2);
                                       tot=rotMatrixInv1*tot'; 
                                       tot=tot';
                                       tot(:,1)=tot(:,1)+center(1);
                                       tot(:,2)=tot(:,2)+center(2);
                        
                                end

                                zVector=perimInner{1,n}(1,3)*ones(size(tot,1),1);
                                infill{w,n}=[tot zVector];
                                tot=[];
                                w=w+1;
                             
                                sharedvalue=intersect(x(:,1),xi);
                                %cerco tra i punti rimasti del perimetro qualcuno da
                                %cui ripartire

                                if ~isempty(sharedvalue)
                                    idx=find(x(:,1)==sharedvalue(1));
                                    tot=[tot;x(idx(1),:)];
                                    x(idx(1),:)=[nan nan];
                                    idx=idx(2:end);
                                    val=size(idx,1);
                                else
                                    flag=1;                       
                                    % esco da tutti i cicli for e finisce l'algoritmo
                                end
                                break
                            end
                        end
                    elseif a<0 
                            %scorro verso il basso
                                if  isempty(find(idx==indexLast, 1)) && ~isnan(x(idx+1,1))
                                    tot=[tot;x(idx+1,:)];
                                    x(idx+1,:)=[nan nan];
                                    idx=idx+1;
                                elseif ~isempty(find(idx==indexLast, 1))
                                    ii=find(idx==indexLast, 1);
                                    if ~isnan(x(indexFirst(ii),1))
                                        tot=[tot;x(indexFirst(ii),:)];
                                        x(indexFirst(ii),:)=[nan nan];
                                        idx=indexFirst(ii);
                                    end
                                else
                                    if check>0
                                       
                                           tot(:,1)=tot(:,1)-center(1);
                                           tot(:,2)=tot(:,2)-center(2);
                                           tot=rotMatrixInv2*tot'; 
                                           tot=tot';
                                           tot(:,1)=tot(:,1)+center(1);
                                           tot(:,2)=tot(:,2)+center(2);
                                     else
                                           tot(:,1)=tot(:,1)-center(1);
                                           tot(:,2)=tot(:,2)-center(2);
                                           tot=rotMatrixInv1*tot'; 
                                           tot=tot';
                                           tot(:,1)=tot(:,1)+center(1);
                                           tot(:,2)=tot(:,2)+center(2);
                                    
                                    end
                                    
                                    zVector=perimInner{1,n}(1,3)*ones(size(tot,1),1);
                                    infill{w,n}=[tot zVector];
                                    tot=[];
                                    w=w+1;
                                    
                                    sharedvalue=intersect(x(:,1),xi);
                                    %cerco tra i punti rimasti del perimetro qualcuno da
                                    %cui ripartire

                                    if ~isempty(sharedvalue)
                                        idx=find(x(:,1)==sharedvalue(1));
                                        tot=[tot;x(idx(1),:)];
                                        x(idx(1),:)=[nan nan];
                                        idx=idx(2:end);
                                        val=size(idx,1);
                                    else
                                        flag=1; 
                                        % esco da tutti i cicli for e finisce l'algoritmo
                                    end
                                    break
                                end
                    
                            for j=1:size(x,1) 
                                if ~isempty(find(xi==tot(end,1), 1))
                                    idx=find(x(:,1)==tot(end,1));
                                    val=size(idx,1);
                                    break
                                        
                                elseif  isempty(find(idx==indexLast, 1))  && ~isnan(x(idx+1,1))
                                    tot=[tot;x(idx+1,:)];
                                    x(idx+1,:)=[nan nan];
                                    idx=idx+1;
                                
                                elseif  ~isempty(find(idx==indexLast, 1)) 
                                    ii=find(idx==indexLast, 1);
                                    if ~isnan(x(indexFirst(ii),1))
                                        tot=[tot;x(indexFirst(ii),:)];
                                        x(indexFirst(ii),:)=[nan nan];
                                        idx=indexFirst(ii);
                                    end
                                else
                                    if check>0
                                    
                                       tot(:,1)=tot(:,1)-center(1);
                                       tot(:,2)=tot(:,2)-center(2);
                                       tot=rotMatrixInv2*tot'; 
                                       tot=tot';
                                       tot(:,1)=tot(:,1)+center(1);
                                       tot(:,2)=tot(:,2)+center(2);
                                    else
                                       tot(:,1)=tot(:,1)-center(1);
                                       tot(:,2)=tot(:,2)-center(2);
                                       tot=rotMatrixInv1*tot'; 
                                       tot=tot';
                                       tot(:,1)=tot(:,1)+center(1);
                                       tot(:,2)=tot(:,2)+center(2);
                
                                    end
                                    
                                    zVector=perimInner{1,n}(1,3)*ones(size(tot,1),1);
                                    infill{w,n}=[tot zVector];
                                    tot=[];
                                    w=w+1;
                                    
                                    sharedvalue=intersect(x(:,1),xi);
                                    %cerco tra i punti rimasti del perimetro qualcuno da
                                    %cui ripartire
            
                                    if ~isempty(sharedvalue)
                                        idx=find(x(:,1)==sharedvalue(1));
                                        tot=[tot;x(idx(1),:)];
                                        x(idx(1),:)=[nan nan];
                                        idx=idx(2:end);
                                        val=size(idx,1);
                                    else
                                        flag=1; 
                                        % esco da tutti i cicli for e finisce l'algoritmo
                                    end
                                    break
                                end
                            end   
                    else

                        if check>0
                        
                            tot(:,1)=tot(:,1)-center(1);
                            tot(:,2)=tot(:,2)-center(2);
                            tot=rotMatrixInv2*tot'; 
                            tot=tot';
                            tot(:,1)=tot(:,1)+center(1);
                            tot(:,2)=tot(:,2)+center(2);
                        else
                            tot(:,1)=tot(:,1)-center(1);
                            tot(:,2)=tot(:,2)-center(2);
                            tot=rotMatrixInv1*tot'; 
                            tot=tot';
                            tot(:,1)=tot(:,1)+center(1);
                            tot(:,2)=tot(:,2)+center(2);
                        
                        end
                        
                        zVector=perimInner{1,n}(1,3)*ones(size(tot,1),1);
                        infill{w,n}=[tot zVector];
                        tot=[];
                        w=w+1;
                        
                        sharedvalue=intersect(x(:,1),xi);

                        if ~isempty(sharedvalue)
                            idx=find(x(:,1)==sharedvalue(1));
                            tot=[tot;x(idx(1),:)];
                            x(idx(1),:)=[nan nan];
                            idx=idx(2:end);
                            val=size(idx,1);
                        else
                            flag=1;     
                        end
                    end 
                end
                if flag==1
                    break
                end 
            end
        end            
    end
end