
% Mesh Creation
w = 2;
l = 2;
Nx = 30;
Ny = 30;
x = linspace(-w/2,w/2, Nx);
y = linspace(-l/2,l/2, Ny);
[X, Y] = meshgrid(x, y) ;
X = X(:) ;
Y = Y(:) ;
TR = delaunay(X, Y) ;
t = TR';
t(4, :) = 1;
p = [X, Y, zeros(length(X),1)].';
figure;
trimesh(TR, X, Y, zeros(1, length (X)),'EdgeColor','b');
xlabel('x');
ylabel('y');
zlabel('z');

%RWG1 Geometry calculations

TrianglesTotal=length(t);

%Find areas of separate triangles
for m=1:TrianglesTotal
   N=t(1:3,m);
   Vec1=p(:,N(1))-p(:,N(2));
   Vec2=p(:,N(3))-p(:,N(2));
   Area(m) =norm(cross(Vec1,Vec2))/2;
   Center(:,m)=1/3*sum(p(:,N),2);
end

%Find all edge elements "Edge_" with at least two adjacent triangles
Edge_=[];
n=0;
for m=1:TrianglesTotal
    N=t(1:3,m);
    for k=m+1:TrianglesTotal
        M=t(1:3,k);      
        a=1-all([N-M(1) N-M(2) N-M(3)]);
        if(sum(a)==2) %triangles m and k have common edge
            n=n+1;
            Edge_=[Edge_ M(find(a))]; 
            TrianglePlus(n)=m;
            TriangleMinus(n)=k; 
        end; 
    end
end
EdgesTotal=length(Edge_);


%All structures of this chapter have EdgeIndicator=2
EdgeIndicator=t(4,TrianglePlus)+t(4,TriangleMinus);

%Find edge length
for m=1:EdgesTotal
   EdgeLength(m)=norm(p(:,Edge_(1,m))-p(:,Edge_(2,m)));
end

%Find nine sub-triangle midpoints 
IMT=[];
for m=1:TrianglesTotal
    n1=t(1,m);
    n2=t(2,m);
    n3=t(3,m); 
    M=Center(:,m);
    r1=    p(:,n1);
    r2=    p(:,n2);
    r3=    p(:,n3);
    r12=r2-r1;
    r23=r3-r2;
    r13=r3-r1;
    C1=r1+(1/3)*r12;
    C2=r1+(2/3)*r12;
    C3=r2+(1/3)*r23;
    C4=r2+(2/3)*r23;
    C5=r1+(1/3)*r13;
    C6=r1+(2/3)*r13;
    a1=1/3*(C1+C5+r1);
    a2=1/3*(C1+C2+M);
    a3=1/3*(C2+C3+r2);
    a4=1/3*(C2+C3+M);
    a5=1/3*(C3+C4+M);
    a6=1/3*(C1+C5+M);
    a7=1/3*(C5+C6+M);
    a8=1/3*(C4+C6+M);
    a9=1/3*(C4+C6+r3);
    Center_(:,:,m)=...
        [a1 a2 a3 a4 a5 a6 a7 a8 a9];
end
%PLUS
for m=1:EdgesTotal
    NoPlus=TrianglePlus(m);
    n1=t(1,NoPlus);
    n2=t(2,NoPlus);
    n3=t(3,NoPlus); 
    if((n1~=Edge_(1,m))&(n1~=Edge_(2,m))) NODE=n1; end;
    if((n2~=Edge_(1,m))&(n2~=Edge_(2,m))) NODE=n2; end;
    if((n3~=Edge_(1,m))&(n3~=Edge_(2,m))) NODE=n3; end;
    FreeVertex=p(:,NODE);
    
    RHO_Plus(:,m)   =+Center(:,NoPlus)-FreeVertex;
    %Nine rho's of the "plus" triangle
    RHO__Plus(:,:,m)  =...
        +Center_(:,:,NoPlus)-repmat(FreeVertex,[1 9]);
end
%MINUS
for m=1:EdgesTotal
    NoMinus=TriangleMinus(m);
    n1=t(1,NoMinus);
    n2=t(2,NoMinus);
    n3=t(3,NoMinus); 
    if((n1~=Edge_(1,m))&(n1~=Edge_(2,m))) NODE=n1; end;
    if((n2~=Edge_(1,m))&(n2~=Edge_(2,m))) NODE=n2; end;
    if((n3~=Edge_(1,m))&(n3~=Edge_(2,m))) NODE=n3; end;
    FreeVertex=p(:,NODE);
    
    RHO_Minus(:,m)   =-Center(:,NoMinus) +FreeVertex;
    %Nine rho's of the "minus" triangle
    RHO__Minus(:,:,m)=...
        -Center_(:,:,NoMinus)+repmat(FreeVertex,[1 9]);
end

%Save result
save RWGM   p ...
            t ...            
            TrianglesTotal ...
            EdgesTotal ...
            Edge_ ...
            TrianglePlus ...
            TriangleMinus ...
            EdgeLength ...
            EdgeIndicator ...
            Area ...
            RHO_Plus ...
            RHO_Minus ...
            RHO__Plus ...
            RHO__Minus ...
            Center ...
            Center_