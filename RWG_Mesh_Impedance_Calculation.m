%RWG_Mesh_Impedance_Calculation Calculates the impedance matrix using
%function impmet
%   Uses the mesh file from RWG_Mesh_Geometry_Calculation, RWGM.mat, as an input.

clear all

%Load the data
load('RWGM');

%EM parameters 
f           =25e9;  
epsilon_    =8.854e-012;
mu_         =1.257e-006;

%Speed of light
c_=1/sqrt(epsilon_*mu_);

%Free-space impedance 
eta_=sqrt(mu_/epsilon_);

%Contemporary variables - introduced to speed up 
%the impedance matrix calculation
omega       =2*pi*f;                                            
k           =omega/c_;
K           =j*k;
Constant1   =mu_/(4*pi);
Constant2   =1/(j*4*pi*omega*epsilon_);
Factor      =1/9;    

FactorA     =Factor*(j*omega*EdgeLength/4)*Constant1;
FactorFi    =Factor*EdgeLength*Constant2;

for m=1:EdgesTotal
    RHO_P(:,:,m)=repmat(RHO_Plus(:,m),[1 9]);   %[3 9 EdgesTotal]
    RHO_M(:,:,m)=repmat(RHO_Minus(:,m),[1 9]);  %[3 9 EdgesTotal]
end
FactorA=FactorA.';
FactorFi=FactorFi.';

%Impedance matrix Z

Z=  impmet( EdgesTotal,TrianglesTotal,...
            EdgeLength,K,...
            Center,Center_,...
            TrianglePlus,TriangleMinus,...
            RHO_P,RHO_M,...
            RHO__Plus,RHO__Minus,...
            FactorA,FactorFi);

%Save result
save IMP   f ...
           omega ...            
           mu_ ...
           epsilon_ ...
           c_ ...
           eta_ ...
           Z ...