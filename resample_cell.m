function new=resample_cell(old,d)

new=cell(size(old,1),size(old,2)); % nuovi punti infittiti

for i=1:size(old,1)
    for j=1:size(old,2)
        if size(old{i,j})>1
            temp=[];
            for k=1:size(old{i,j},1)-1
                    iniz=old{i,j}(k,1:2);
                    fin=old{i,j}(k+1,1:2);
                    dist=sqrt((fin(1)-iniz(1))^2+(fin(2)-iniz(2))^2);
               if (dist/d)>3
                    x=linspace(old{i,j}(k,1),old{i,j}(k+1,1),dist/d)';
                    y=linspace(old{i,j}(k,2),old{i,j}(k+1,2),dist/d)';
                    z1=linspace(old{i,j}(k,3),old{i,j}(k+1,3),dist/d)';
                    temp=[temp;x(1:end-1) y(1:end-1) z1(1:end-1)]; 
               else
                    temp=[temp;old{i,j}(k,:)]; 
               end
            end
            new{i,j}=[temp;old{i,j}(end,:)];
            new{i,j}=round(new{i,j},5);
        end
    end
end