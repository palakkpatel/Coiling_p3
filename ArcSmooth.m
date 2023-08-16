function [curve]=ArcSmooth(p1,q1,r1,nodes)
% Smoothing the Arcs. That is the bend between two segments
    pr=r1-p1;
    qp=p1-q1;
    qr=r1-q1;
    %% Rotate Normal to z axis
    n=-cross(-qp,pr);
    n=n/norm(n);
    z=[0 0 1];
    theta=acos(dot(n,z)/norm(n)/norm(z));
    axis=cross(n,z);
    axis=axis/norm(axis);
    
    u=axis;
    c=cos(theta);
    s=sin(theta);
    c1=1-c;

    R(:,1)=[c+u(1)*u(1)*c1; u(1)*u(2)*c1+u(3)*s; u(1)*u(3)*c1-u(2)*s];
    R(:,2)=[u(1)*u(2)*c1-u(3)*s; c+u(2)*u(2)*c1; u(2)*u(3)*c1+u(1)*s];
    R(:,3)=[u(1)*u(3)*c1+u(2)*s; u(2)*u(3)*c1-u(1)*s; c+u(3)*u(3)*c1];
    
    %% Rotate PR to x-axis
    pr_r=pr*R';
    pq_r=(-qp)*R';
    
    x=[1 0 0];
    theta1=acos(dot(pr_r,x)/norm(pr_r)/norm(x));
    axis=cross(pr_r,x);
    axis=axis/norm(axis);
    
    u=axis;
    c=cos(theta1);
    s=sin(theta1);
    c1=1-c;
    
    R1(:,1)=[c+u(1)*u(1)*c1; u(1)*u(2)*c1+u(3)*s; u(1)*u(3)*c1-u(2)*s];
    R1(:,2)=[u(1)*u(2)*c1-u(3)*s; c+u(2)*u(2)*c1; u(2)*u(3)*c1+u(1)*s];
    R1(:,3)=[u(1)*u(3)*c1+u(2)*s; u(2)*u(3)*c1-u(1)*s; c+u(3)*u(3)*c1];
    
    pr_r=pr_r*R1';
    pq_r=pq_r*R1';
    
    %% Arc Length
    p=[0 0];
    q=pr_r(1,1:2);
    s0=norm(qp)+norm(qr);
    e=10^-10;
    N=10000;
    g=10;

    for j=1:N
        sg=S(q,g);
        ds=dS(q,g);
        dg=(sg-s0)/ds;
        g=g-dg;
        if abs(dg)<=e
            break;
        end
    end
    a=g;
    b=q(2)/q(1)-a*q(1);
    x=linspace(p(1),q(1),nodes);
    y=a.*x.^2+b.*x;
%     ii=[x(5) y(5)];
%     jj=pq_r(1,1:2);
%     dd=ii-jj;
%     dis=sqrt(sum(dd.^2));
    if y(nodes/2)*pq_r(2)<0
        y=-y;
    end
    curve=zeros(nodes,3);
    curve(:,1:2)=[x' y'];
    curve=curve*R1;
    curve=curve*R;
    curve=[curve(:,1)+p1(1) curve(:,2)+p1(2) curve(:,3)+p1(3)];
    
%     plot3(curve(:,1),curve(:,2),curve(:,3));
    %% UDFs
    function [di]=dI(t)
        di=sqrt(1+t^2);
    end

    function [i]=I(t)
        r=sqrt(1+t^2);
        i=0.5*(t * r + log(t + r));
    end

    function [s]=S(q,a)
        u1=q(2)/q(1) + a*q(1);
        l1=q(2)/q(1) - a*q(1);
        iu=I(u1);
        il=I(l1);
        s=0.5*(iu-il)/a;
    end

    function ds=dS(q,a)
        u1=q(2)/q(1) + a*q(1);
        l1=q(2)/q(1) - a*q(1);
        diu=dI(u1);
        dil=dI(l1);
        iu=I(u1);
        il=I(l1);
        ds=0.5*(a*q(1)*(diu+dil)+il-iu)/(a*a);
    end

end
