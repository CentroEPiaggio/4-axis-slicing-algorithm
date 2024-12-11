function poly=cell2poly(shell)

poly=cell(1,size(shell,2)); 

for i=1:size(shell,2)
       if ~isempty(shell{1,i})
         X=shell{1,i}(:,1);
         Y=shell{1,i}(:,2);
         poly{1,i}=polyshape(X,Y);
       end
end