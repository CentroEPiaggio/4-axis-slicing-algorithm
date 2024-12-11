function [shell,triTot,triThresUp,triThresDown]=head_tail_true(tri,z,thresUp,thresDown,TR)
shell= {};  
w=1;
triTot=[];

for i=1:length(z)
      punti = [];
      intersections=[];
      k=1;
      T3=tri(:,3);
      T6=tri(:,6);
      T9=tri(:,9);
      indexTri=find(T3<=z(i) & T6>=z(i) | T3>=z(i) & T6<=z(i) | ...
          T6<=z(i) & T9>=z(i)| T6>=z(i) & T9<= z(i) |...
          T3>=z(i) & T9<=z(i)| T3<=z(i) & T9>=z(i));

      indexTriCop=[];
      indexTriCop=find(tri(:,3)==z(i) & tri(:,6)==z(i) & tri(:,9)==z(i));
      zTrue=z(i);

      if ~isempty(indexTriCop)

          V1=tri(indexTriCop(1,1),1:3);
          V2=tri(indexTriCop(1,1),4:6);
          V3=tri(indexTriCop(1,1),7:9);

          normale=cross(V1-V3,V2-V3);
 
          Area=norm(normale)/2;
          normale=normale/(2*Area);

          while ~isempty(indexTriCop)
              if normale(1,3)==-1
                  zTrue=zTrue+10^-8;
              else
                  zTrue=zTrue-10^-8;
              end
              indexTriCop=find(tri(:,3)==zTrue & tri(:,6)==zTrue & tri(:,9)==zTrue);
          end
          indexTri=find(T3<=zTrue & T6>=zTrue | T3>=zTrue & T6<=zTrue | ...
          T6<=zTrue & T9>=zTrue| T6>=zTrue & T9<= zTrue |...
          T3>=zTrue & T9<=zTrue| T3<=zTrue & T9>=zTrue);
      end

      k=1;
      for j=1:size(indexTri,1)
          p3=tri(indexTri(j),3);
          p6=tri(indexTri(j),6);
          p9=tri(indexTri(j),9);

          % primo if - se il piano intereseca un solo punto non lo considerare
          % if successivi - se il piano interseca due vertici aggiungi l'edge tra i due vertici
          % ultimo if - 

            if p3==zTrue &&  p6 >zTrue && p9 >zTrue ||...
                    p3==zTrue &&  p6 <zTrue && p9 <zTrue ||...
                    p3>zTrue &&  p6 ==zTrue && p9>zTrue || ...
                    p3<zTrue &&  p6 ==zTrue && p9<zTrue || ...
                    p3>zTrue &&  p6>zTrue && p9 ==zTrue || ...
                    p3<zTrue &&  p6<zTrue && p9 ==zTrue
                    j=j+1;
            elseif p3==zTrue && p6==zTrue && p9~=zTrue
                intersections(k,:)=[tri(indexTri(j),1:3) tri(indexTri(j),4:6)];
                k=k+1;
            elseif p3~=zTrue && p6==zTrue && p9==zTrue
                intersections(k,:)=[tri(indexTri(j),4:6) tri(indexTri(j),7:9)];
                k=k+1;
            elseif p3==zTrue && p6~=zTrue && p9==zTrue
                intersections(k,:)=[tri(indexTri(j),1:3) tri(indexTri(j),7:9)];
                k=k+1;
            else
               [flag(k), intersections(k,:)] =plane_tri_inter([0 0 zTrue], [0 0 1], tri(indexTri(j),1:3), tri(indexTri(j),4:6), tri(indexTri(j),7:9));
               k=k+1;
            end

            if ~isempty(intersections)
                temp_int=[intersections(end,4:6) intersections(end,1:3)];
                [~,~,iA1]=intersect(intersections(end,:),intersections,"rows");
                [~,~,iA2]=intersect(temp_int,intersections,"rows");
    
                if (length(iA1)>=1 && length(iA2)>=1) || length(iA1)>1
                    intersections(end,:)=[];
                    k=size(intersections,1)+1;
                end
            end
      end
      
      if ~isempty(intersections) 

          triTot(i,1)=size(indexTri,1); % n di triangoli totali su quel layer
          TR_temp=triangulation(TR.ConnectivityList(indexTri,:),TR.Points);
          F = faceNormal(TR_temp); %normali alle face dei triangoli
          u=[0 0 1]; %normale all'asse z
          ThetaInDegrees=zeros(length(F),1); % calcolo gli angoli tra le normali e l'asse z
          k=1;

           for j=1:length(F)
              CosTheta = max(min(dot(u,F(j,:))/(norm(u)*norm(F(j,:))),1),-1);
              ThetaInDegrees(k) = real(acosd(CosTheta));
              k=k+1;
           end
 
          UP=find(ThetaInDegrees(:)<=thresUp);
          DOWN=find(ThetaInDegrees(:)>=thresDown);

          triThresDown(i,1)=size(DOWN,1);
          triThresUp(i,1)=size(UP,1);
    
          intersections(:,3)=round(intersections(:,3),7);
          intersections(:,6)=round(intersections(:,6),7);
          
          temp1=nan*ones(size(intersections,1),1);
          temp2=nan*ones(size(intersections,1),1);

          punti(1,:)=intersections(1,1:3);
          punti(2,:)=intersections(1,4:6);
          intersections(1,4:6)=nan*ones(1,3);
          
          n=2;

          for j=1:size(intersections,1)-1
             for h=1:size(intersections,1)
                temp1(h)=sqrt((punti(n,1)-intersections(h,1))^2+(punti(n,2)-intersections(h,2))^2);
                temp2(h)=sqrt((punti(n,1)-intersections(h,4))^2+(punti(n,2)-intersections(h,5))^2);
             end
             
            id=[];
            [a,ii]=min(temp1);
            [b,jj]=min(temp2);
            c=min(a,b);

            if c==a && ~isnan(intersections(ii,4)) && ~isnan(intersections(ii,5)) && ~isnan(intersections(ii,6))
                 n=n+1;
                 punti(n,:)=intersections(ii,4:6);
                 intersections(ii,:)=nan*ones(1,6);
               
            elseif c==a && isnan(intersections(ii,4)) && isnan(intersections(ii,5)) && isnan(intersections(ii,6))
                 n=n+1;
                 punti(n,:)=[nan nan nan];
                 intersections(ii,1:3)=[nan nan nan];
                 id=find(~isnan(intersections(:,1)),1);        

            elseif c==b && ~isnan(intersections(jj,1)) && ~isnan(intersections(jj,2)) && ~isnan(intersections(jj,3))
                  n=n+1;
                  punti(n,:)=intersections(jj,1:3);
                  intersections(jj,:)=nan*ones(1,6);
               
            else 
                  n=n+1;
                  punti(n,:)=[nan nan nan];
                  intersections(ii,4:6)=[nan nan nan];
                  id=find(~isnan(intersections(:,1)),1);
            end
            
             if ~isempty(id)
                  n=n+1;
                  punti(n,:)=intersections(id,1:3);
                  n=n+1;
                  punti(n,:)=intersections(id,4:6);  
                  intersections(id,4:6)=[nan nan nan];
             end
          end

          shell{1,w} = [punti; nan nan nan];
          w=w+1;
      else 
          shell{1,w}={};
          w=w+1;
      end
end 
