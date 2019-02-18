clear all
clc


% initial value of twist
uz_0=[0 0 0];
uxy_0=[0 0];
% desired joint positions
q=[0 0 0 0 0 0];

R =[ 1.0000         0         0; 0    0.0404   -0.9992; 0    0.9992    0.0404];
F=R*[0 0 0]';
% initial position of joints
z0=[-0.2858   -0.2025   -0.0945  0 0 0];
z=z0+q;

[r1,r2,r3,Uz] =moving_CTR(z,uz_0,uxy_0,F);



plot3(r1(:,1),r1(:,2),r1(:,3),'linewidth',1)
hold on
plot3(r2(:,1),r2(:,2),r2(:,3),'linewidth',2)
plot3(r3(:,1),r3(:,2),r3(:,3),'linewidth',3)
xlabel('X')
ylabel('Y')
zlabel('Z')
axis equal

%% main ode solver
function [r1,r2,r3,Uz] = moving_CTR(q,uz_0,uxy_0,F)

E=[ 6.4359738368e+10,5.2548578304e+10,4.7163091968e+10];
J= 1.0e-11 *[ 0.0120    0.0653    0.1686];
I= 1.0e-12 *[ 0.0601    0.3267    0.8432 ];
G=[2.5091302912e+10,2.1467424256e+10,2.9788923392e+10]; 
Ux=[21.3 13.108 3.5];
Uy=[0 0 0];
n=3;
 
% length of tubes
l=1e-3.*[328+103  219+113 174];
% length of the curved part of tubes
l_k=1e-3*[103 113 174-40];

% q1 o q3 are robot base movments, q3 to q6 are rbot base rotation angle.
if ~(iscolumn(uz_0))
    uz0=uz_0';
else
    uz0=uz_0;
end
if ~(iscolumn(uxy_0))
    uxy_0=uxy_0';
else
    uxy_0=uxy_0;
end

B=[q(1)  q(2)  q(3)];  % length of tubes before template
%initial angles

alpha=[0 q(5)-q(4) q(6)-q(4)]-B.*uz0';
alpha_1=q(4);

% segmenting tubes  
% check all inputs must have n elements, n is number of tubes
[L,d_tip,EE,UUx,UUy] = segmenting(E,Ux,Uy,l,B,l_k);

SS=L;
for i=1:length(L)
    SS(i)=sum(L(1:i));
%     plot((B(1)+SS(i))*ones(1,10),1:10,'b' ,'LineWidth',2)
end

% S is segmented abssica of tube after template
 S=SS(SS+min(B)>0)+min(B);
 E=zeros(n,length(S)); Ux=E; Uy=E;
 for i=1:n
    E(i,:)=EE(i,SS+min(B)>0); Ux(i,:)=UUx(i,SS+min(B)>0); Uy(i,:)=UUy(i,SS+min(B)>0);
 end
 % each (i,j) element of above matrices correspond to the jth segment of
 % ith tube, 1st tube is the most inner

 %% Solving ode for shape

span=[0 S];       % vector of tube abssica starting at zero
Length=[]; r=[]; U_z=[]; a=[];% solved length, curvatures, and twist angles
%U1_after=[0;0;0];             % 1st tube initial curvature at segment beginning
r0=[ 0 0 0]'; R0=[cos(alpha_1) sin(alpha_1) 0; -sin(alpha_1) cos(alpha_1) 0; 0 0 1];
R0=reshape(R0,[9,1]);
%alpha=alpha-B.*uz_0'; 


% initial curvature of tube 1
EI=E(:,1).*I';
uxy_0(1)= (1/(EI(1)+EI(2)+EI(3)))* (...
      EI(1)*Ux(1)+ ...
      EI(2)*Ux(2)*cos(alpha(1)-alpha(2))+ EI(2)*Uy(2)*sin(alpha(1)-alpha(2))  + ...
      EI(3)*Ux(3)*cos(alpha(1)-alpha(3))+ EI(3)*Uy(3)*sin(alpha(1)-alpha(3)) );
uxy_0(2)= (1/(EI(1)+EI(2)+EI(3)))* (...
      + EI(1)*Uy(1)  + ...
      -EI(2)*Ux(2)*sin(alpha(1)-alpha(2))+ EI(2)*Uy(2)*cos(alpha(1)-alpha(2))  + ...
      -EI(3)*Ux(3)*sin(alpha(1)-alpha(3))+ EI(3)*Uy(3)*cos(alpha(1)-alpha(3)) );

for seg=1:length(S)
    
s_span = [span(seg) span(seg+1)-0.0000001];
y0_1=[r0 ; R0];

y0_2=zeros(2*n,1);
y0_2(n+1:2*n)=alpha;
y0_2(1:n)=uz0;
y0_3=uxy_0;

y_0=[y0_3;y0_2;y0_1];

[s,y] = ode23(@(s,y) ode(s,y,Ux(:,seg),Uy(:,seg),E(:,seg).*I',G.*J,F,n), s_span, y_0);
% first n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% last n elements of y are twist angles, alpha_i
shape=[y(:,2+2*n+1),y(:,2+2*n+2),y(:,2+2*n+3)];
Length=[Length; s];
r=[r; shape];
U_z=[U_z; y(:,2+1:2+n )];
a=[a; y(:,2+n+1:2+2*n )]; % twist angle

%calculating new boundary conditions
r0=shape(end,:)';
R0=y(end,2+2*n+4:2+2*n+12)';
alpha=y(end,2+n+1:2+2*n);
u1=[y(end,1) y(end,2) y(end,3)]';


dtet2=y(end,2+2)-y(end,3);
dtet3=y(end,2+3)-y(end,3);

GG=G';
GG(E(i,seg)==0)=0;

K1=diag( [E(1,seg).*I(1) E(1,seg).*I(1) GG(1).*J(1)] );
K2=diag( [E(2,seg).*I(2) E(2,seg).*I(2) GG(2).*J(2)] );
K3=diag( [E(3,seg).*I(3) E(3,seg).*I(3) GG(3).*J(3)] );
U1=[Ux(1,seg) Uy(1,seg) 0]'; U2=[Ux(2,seg) Uy(2,seg) 0]'; U3=[Ux(3,seg) Uy(3,seg) 0]';
R=[y(end,2+2*n+4) y(end,2+2*n+5) y(end,2+2*n+6);y(end,2+2*n+7) y(end,2+2*n+8) y(end,2+2*n+9);y(end,2+2*n+10) y(end,2+2*n+11) y(end,2+2*n+12)]
R_tet2=[ cos(y(end,2+n+2)) -sin(y(end,2+n+2)) 0; sin(y(end,2+n+2)) cos(y(end,2+n+2)) 0; 0 0 1]';
R_tet3=[ cos(y(end,2+n+3)) -sin(y(end,2+n+3)) 0; sin(y(end,2+n+3)) cos(y(end,2+n+3)) 0; 0 0 1];
e3=[0 0 1]';
u2=R_tet2'*u1+dtet2*e3;
u3=R_tet3'*u1+dtet3*e3;
if seg<length(S)

GG=G';
GG(E(i,seg+1)==0)=0;

K1_new=diag( [E(1,seg+1).*I(1) E(1,seg+1).*I(1) GG(1).*J(1)] );
K2_new=diag( [E(2,seg+1).*I(2) E(2,seg+1).*I(2) GG(2).*J(2)] );
K3_new=diag( [E(3,seg+1).*I(3) E(3,seg+1).*I(3) GG(3).*J(3)] );
U1_new=[Ux(1,seg+1) Uy(1,seg+1) 0]'; U2_new=[Ux(2,seg+1) Uy(2,seg+1) 0]'; U3_new=[Ux(3,seg+1) Uy(3,seg+1) 0]';

%u1_new= (   K1_new + K2_new + K3_new   )\(K1_new*U1_new+R_tet2*K2_new*U2_new+R_tet3*K3_new*U3_new);
 
u1_new= (   K1_new + K2_new + K3_new   )...
     \(  K1*(u1-U1)+R_tet2*K2*(u2-U2)+R_tet3*K3*(u3-U3)+K1_new*U1_new+R_tet2*K2_new*U2_new+R_tet3*K3_new*U3_new ...
     -R_tet2*K2_new*dtet2*e3-R_tet3*K3_new*dtet3*e3  );
 
% u1_new= (   K1_new + K2_new + K3_new   )...
%     \(  K1*(u1-U1)+R_tet2*K2*(R_tet2'*u1-U2)+R_tet3*K3*(R_tet2'*u1-U3)+K1_new*U1_new+R_tet2*K2_new*U2_new+R_tet3*K3_new*U3_new )

u2_new=R_tet2'*u1_new+dtet2*e3;
u3_new=R_tet3'*u1_new+dtet3*e3;
uz0=[u1_new(3) u2_new(3) u3_new(3)];
uz0=U_z(end,:)';
uxy_0=[u1_new(1) u1_new(2)]';

end
end

Uz=zeros(n,1);
for i=1:n
[~,index] =  min( abs(Length-d_tip(i)+0.0001) );
Uz(i)= U_z(index,i);
end

r1=r;
[~, tube2_end] = min(abs(Length-d_tip(2)));
r2=[r(1:tube2_end,1),r(1:tube2_end,2),r(1:tube2_end,3)];
[~, tube3_end] = min(abs(Length-d_tip(3)));
r3=[r(1:tube3_end,1),r(1:tube3_end,2),r(1:tube3_end,3)];
end


%% ODE
function dydt = ode(~,y,Ux,Uy,EI,GJ,F,n)

dydt=zeros(2+2*n+12,1);
% first element of y is curvature along x for first tube,
% 2nd element of y is curvature along y for first tube
% next n elements of y are curvatures along z, e.g., y= [ u1_z  u2_z ... ]
% second n elements of y are twist angles, alpha_i
% last 12 elements are r (position) and R (orientations), respectively

% calculating 2nd and 3rd tube's curvatures
tet2=y(2+n+2); tet3=y(2+n+3);
u2=  [ cos(tet2) -sin(tet2) 0; sin(tet2) cos(tet2) 0; 0 0 1]'*[y(1) y(2) y(3)]'+ dydt(2+n+2)*[0 0 1]';
u3=  [ cos(tet3) -sin(tet3) 0; sin(tet3) cos(tet3) 0; 0 0 1]'*[y(1) y(2) y(3)]'+ dydt(2+n+3)*[0 0 1]';
% u=[u1x u1y u1z u2x u2y u2z u3x u3y u3z] is vector of curvatures
u=[y(1) y(2) y(3) u2(1) u2(2) y(4) u3(1) u3(2) y(5)]';


% odes for twist
GJ=GJ';
GJ(EI==0)=0;
for i=1:n
    if EI(i)==0
    dydt(2+i)=  0;  % ui_z assuming no pre twist in tube
    dydt(2+n+i)= 0;   %alpha_i
    else
    dydt(2+i)=  (  (EI(i))/(GJ(i))  ) * ( u((i-1)*n+1)* Uy(i) -  u((i-1)*n+2)* Ux(i) );  % ui_z assuming no pre twist in tube
    dydt(2+n+i)=  y(2+i)-y(3);   %alpha_i
    end
end

% calculating 1st tube's curvatures in x and y direction
K_inv=diag(  [1/sum(EI) 1/sum(EI) 1/sum(GJ) ]  );
duxy=zeros(n,1);
for i=1:n
   duxy=duxy+ [ cos(y(2+n+i)) -sin(y(2+n+i)) 0; sin(y(2+n+i)) cos(y(2+n+i)) 0; 0 0 1]* ( diag( [EI(i) EI(i) GJ(i)] ) * (...
     dydt(2+n+i) *   [ -sin(y(2+n+i)) -cos(y(2+n+i)) 0; cos(y(2+n+i)) -sin(y(2+n+i)) 0; 0 0 1].' * [y(1) y(2) y(3)].' ) + ...
   [0 -u(3*(i-1)+3) u(3*(i-1)+2) ; u(3*(i-1)+3) 0 -u(3*(i-1)+1); -u(3*(i-1)+2) u(3*(i-1)+1) 0]  *diag( [EI(i) EI(i) GJ(i)] ) * ...
   ([u(3*(i-1)+1) u(3*(i-1)+2) u(3*(i-1)+3)].'-[Ux(i) Uy(i) 0].')  );
end

R1=[y(2+2*n+4) y(2+2*n+5) y(2+2*n+6);y(2+2*n+7) y(2+2*n+8) y(2+2*n+9);y(2+2*n+10) y(2+2*n+11) y(2+2*n+12)];

duxy=-K_inv*duxy-K_inv* ([0 -1 0;1 0 0;0 0 0]*R1'*[F(1) F(2) F(3)].');
dydt(1)=duxy(1);
dydt(2)=duxy(2);

% odes for r and R
e3=[0 0 1]';              
u_hat=[0 -y(3) y(2) ; y(3) 0 -y(1) ; -y(2) y(1) 0 ];
dr1 = R1*e3;
dR1=R1*u_hat;


dydt(2+2*n+1)=dr1(1);dydt(2+2*n+2)=dr1(2);dydt(2+2*n+3)=dr1(3);
dR=dR1'; 
dR=dR(:);
for i=4:12
   dydt(2+2*n+i)=dR(i-3);
end

end


%% code for segmenting tubes
function [L,d1,E,Ux,Uy] = segmenting(E,Ux,Uy,l,B,l_k)

% all vectors must be sorted, starting element belongs to the most inner tube
%E, U, I, G, J   stifness, curvature, inertia, torsion constant, and second moment of inertia vectors for each tube
%l vector of tube length
%B  vector of tube movments with respect to template position, i.e., s=0 (always negative)
%l_k vecot oftube's curved part length


d1= l+B; % position of tip of the tubes
d2=d1-l_k; % position of the point where tube bending starts
points=[0 B d2 d1];
[L, index]=sort(points);
L = 1e-5*floor(1e5*diff(L));  % length of each segment 
%(used floor because diff command doesn't give absolute zero sometimes)
% 
% for i=1:k-1
% if B(i)>B(i+1)
%     sprintf('inner tube is clashing into outer tubes')
%     E=zeros(k,length(L));
%     I=E; G=E; J=E; Ux=E; Uy=E;
%     return
% end
% end

EE=zeros(3,length(L)); UUx=EE; UUy=EE;

for i=1:3

[~, a] = min(abs(index-i+1)); % find where tube begins
[~, b] = min(abs(index-(1*3+i+1))); % find where tube curve starts
[~, c] = min(abs(index-(2*3+i+1))); % find where tube ends

if L(a)==0; a=a+1;  end
if L(b)==0; b=b+1;  end
if c<=length(L)
    if L(c)==0; c=c+1; end
end
    
EE(i,a:c-1)=E(i);
UUx(i,b:c-1)=Ux(i);
UUy(i,b:c-1)=Uy(i);
end

l=L(~(L==0));  % get rid of zero lengthes
E=zeros(3,length(l)); Ux=E; Uy=E;
 for i=1:3
    E(i,:)=EE(i,~(L==0)); Ux(i,:)=UUx(i,~(L==0)); Uy(i,:)=UUy(i,~(L==0));
 end
L=L(~(L==0));

end


