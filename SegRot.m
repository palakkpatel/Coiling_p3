function[t1_new]=SegRot(t1,t2,t3,AngRot)
%% Calculate Rotation Axis
a=t3-t2;
b=t1-t2;
u=cross(b,a);
u=u/norm(u);

%% Calculate Rotation Matrix
AngRot=pi*AngRot/180;
c=cos(AngRot);
s=sin(AngRot);
c1=1-c;

R(:,1)=[c+u(1)*u(1)*c1; u(1)*u(2)*c1+u(3)*s; u(1)*u(3)*c1-u(2)*s];
R(:,2)=[u(1)*u(2)*c1-u(3)*s; c+u(2)*u(2)*c1; u(2)*u(3)*c1+u(1)*s];
R(:,3)=[u(1)*u(3)*c1+u(2)*s; u(2)*u(3)*c1-u(1)*s; c+u(3)*u(3)*c1];

%% Rot segment
b1=b*R';
t1_new= t2 + b1;




end