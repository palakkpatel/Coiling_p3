function [meshpoint]=smoothing_pp(coil_r,coils)
    % This function smmoths by sub-dividing the segment in nodes segments and smooths the curve
nodes=30;
meshpoint=coil_r;
meshpoint=meshpoint(:,1:3);
preshape=coils;
smoothcoil=meshpoint;
smoothpre=preshape;

[a_d] = sharpAngles(meshpoint, preshape);
idx=length(meshpoint);
j=1;
jmax=idx;
mm=2;
for i=2:idx-1
    j=j+1;
    if a_d(i)<(80*pi/180) || a_d(i)>(100*pi/180)
        continue;
    else
        if mm==i-1
            decoy=smoothcoil(j-1,:);
            decoy_pre=smoothpre(j-1,:);
        else
            decoy=meshpoint(i-1,:);
            decoy_pre=preshape(i-1,:);
        end
        [curve]=ArcSmooth(decoy,meshpoint(i,:),meshpoint(i+1,:),nodes);
        curve_preshape(:,1)=[(linspace(decoy_pre(1),preshape(i+1,1),nodes))'];
        curve_preshape(:,2)=[(linspace(decoy_pre(2),preshape(i+1,2),nodes))'];
        curve_preshape(:,3)=[(linspace(decoy_pre(3),preshape(i+1,3),nodes))'];
        k=j;
        jmax=jmax+(nodes-3);
        j=j+(nodes-3);        
        smoothcoil(j+1:jmax,:)=meshpoint(i+1:idx,:);
        smoothcoil(k-1:j+1,:)=curve;
        smoothpre(j+1:jmax,:)=preshape(i+1:idx,:);
        smoothpre(k-1:j+1,:)=curve_preshape;
        mm=i;
    end
end

meshpoint=smoothcoil;
%plot3(meshpoint(:,1),meshpoint(:,2),meshpoint(:,3));
end

        
    
    