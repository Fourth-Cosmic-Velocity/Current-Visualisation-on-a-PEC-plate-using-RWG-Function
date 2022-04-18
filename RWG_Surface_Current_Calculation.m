%RWG_Surface_Current_Calculation Solves MoM equations for the scattering problem
%   Inputs: RWGM.mat from RWG_Mesh_Geometry_Calculation, IMP.mat from RWG_Mesh_Impedance_Calculation.

clear all

%load the data
load('RWGM');
load('IMP');
 

%Plate - normal incidence
d       =[0 0 -1];     
Pol     =[1 0 0];      

k=omega/c_;
kv=k*d;

for m=1:EdgesTotal    
   ScalarProduct=sum(kv.*Center(:,TrianglePlus(m))');
   EmPlus =Pol.'*exp(-j*ScalarProduct);      
   ScalarProduct=sum(kv.*Center(:,TriangleMinus(m))');
   EmMinus=Pol.'*exp(-j*ScalarProduct);      
   ScalarPlus =sum(EmPlus.* RHO_Plus(:,m));
   ScalarMinus=sum(EmMinus.*RHO_Minus(:,m));
   V(m)=EdgeLength(m)*(ScalarPlus/2+ScalarMinus/2);   
end

%System solution
I=Z\V.';
save Current   f ...
               omega ...            
               mu_ ...
               epsilon_ ...
               c_ ...
               eta_ ...
               I ...
               V ...
               d ...
               Pol ...      