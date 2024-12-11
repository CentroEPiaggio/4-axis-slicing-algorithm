% ==============================================================
% 
%  N : Normal to plane , |N| = 1, N->(nx3)
%  Q : Point into the plane, Q->(nx3)
% 
%  P1,P2,P3 : Triangle vertices, P1,P2,P3->(nx3)
% 
%  return : 0 if no intersection
%  U : Intersection points, U->(3x2)
% 
% =============================================================

function [flag, U] = plane_tri_inter(Q, N, P1, P2, P3)

   
    d1 = (P1(1)-Q(1))*N(1) + (P1(2)-Q(2))*N(2) + (P1(3)-Q(3))*N(3);
    d2 = (P2(1)-Q(1))*N(1) + (P2(2)-Q(2))*N(2) + (P2(3)-Q(3))*N(3); 
    d3 = (P3(1)-Q(1))*N(1) + (P3(2)-Q(2))*N(2) + (P3(3)-Q(3))*N(3); 

    s1 = 0; 
    s2 = 0; 
    s3 = 0; 
    if d1<0
       s1 = 1;
    end 
    if d2<0
       s2 = 1;
    end 
    if d3<0
       s3 = 1;
    end

    if s2~=s1
        U=zeros(1,6);
        U(1) = (d2*P1(1) - d1*P2(1))/(d2-d1) ;
        U(2) = (d2*P1(2) - d1*P2(2))/(d2-d1) ;
        U(3) = (d2*P1(3) - d1*P2(3))/(d2-d1) ;
        if s1~=s3
            U(4) = (d1*P3(1) - d3*P1(1))/(d1-d3) ;
            U(5) = (d1*P3(2) - d3*P1(2))/(d1-d3) ;
            U(6) = (d1*P3(3) - d3*P1(3))/(d1-d3) ;
         else 
            U(4) = (d2*P3(1) - d3*P2(1))/(d2-d3) ;
            U(5) = (d2*P3(2) - d3*P2(2))/(d2-d3) ;
            U(6) = (d2*P3(3) - d3*P2(3))/(d2-d3) ; 
        end 
    elseif s3~=s1
            U=zeros(1,6);
            U(1) = (d3*P1(1) - d1*P3(1))/(d3-d1) ;
            U(2) = (d3*P1(2) - d1*P3(2))/(d3-d1) ;
            U(3) = (d3*P1(3) - d1*P3(3))/(d3-d1) ;

            U(4) = (d2*P3(1) - d3*P2(1))/(d2-d3) ;
            U(5) = (d2*P3(2) - d3*P2(2))/(d2-d3) ;
            U(6) = (d2*P3(3) - d3*P2(3))/(d2-d3) ;  
    end 
    flag = 1; 
    
     if  s2==s1 && s1==s3 && s3==s1
          U = ones(1, 6)*NaN; 
         flag = 0; 
     end
end 
