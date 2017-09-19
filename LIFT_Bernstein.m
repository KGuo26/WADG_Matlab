N=4;


E54=Elevation2D(BernsteinBasis1D(5),BernsteinBasis1D(4),5);
E43=Elevation2D(BernsteinBasis1D(4),BernsteinBasis1D(3),4);
E32=Elevation2D(BernsteinBasis1D(3),BernsteinBasis1D(2),3);
E21=Elevation2D(BernsteinBasis1D(2),BernsteinBasis1D(1),2);
E10=Elevation2D(BernsteinBasis1D(1),BernsteinBasis1D(0),1);

L0=(N+1)^2/2*E54'*E54;

E1=-nchoosek(N,1)/(2)*E43';
E2=nchoosek(N,2)/(3)*E32'*E43';
E3=-nchoosek(N,3)/(4)*E21'*E32'*E43';
E4=nchoosek(N,4)/(5)*E10'*E21'*E32'*E43';
I=eye(5);
E_1=[I;E1;E2;E3;E4];
L1=E_1*L0;
%face2
E_2([5,9,12,14,15],:)=I;
E_2([4,8,11,13],:)=E1;
E_2([3,7,10],:)=E2;
E_2([2,6],:)=E3;
E_2(1,:)=E4;
L2=E_2*L0;
%face3
E_3([1,6,10,13,15],:)=I;
E_3([2,7,11,14],:)=E1;
E_3([3,8,12],:)=E2;
E_3([4,9],:)=E3;
E_3(5,:)=E4;
L3=E_3*L0;

LIFT1=[L1 L2 L3];


