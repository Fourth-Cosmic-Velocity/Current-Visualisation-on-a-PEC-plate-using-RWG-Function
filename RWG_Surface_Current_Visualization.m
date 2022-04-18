%RWG_Surface_Current_Visualization visualizes the surface current magnitude 
%   Inputs: RWGM.mat from RWG_Mesh_Geometry_Calculation, Current.mat from RWG_Surface_Current_Calculation.
clear all

%load the data
load('RWGM');
load('Current');

Index=find(t(4,:)<=1);
Triangles=length(Index);

%Find the current density for every triangle
for k=1:Triangles
    i=[0 0 0]';
    for m=1:EdgesTotal
        IE=I(m)*EdgeLength(m);
        if(TrianglePlus(m)==k)
            i=i+IE*RHO_Plus(:,m)/(2*Area(TrianglePlus(m)));
        end
        if(TriangleMinus(m)==k)
            i=i+IE*RHO_Minus(:,m)/(2*Area(TriangleMinus(m)));
        end
    end
    CurrentNorm(k)=abs(norm(i));
end

Jmax=max(CurrentNorm);
MaxCurrent=strcat(num2str(Jmax),'[A/m]')
CurrentNorm1=CurrentNorm/max(CurrentNorm);
for m=1:Triangles
    N=t(1:3,m);
    X(1:3,m)=[p(1,N)]';
    Y(1:3,m)=[p(2,N)]';
    Z(1:3,m)=[p(3,N)]';      
end
C=repmat(CurrentNorm1,3,1);


h=fill3(X, Y, Z, C); %linear scale
colormap cool;
axis('equal');
rotate3d

