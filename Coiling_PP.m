clear all;
close all;
clc;

% load coil preshape
No =1;
coils = cell(No,1);
idx_coil=0;
idx_coil = 1;
coordinates=load('CJA_coilpreshape.txt');
%load('./deployment_rob/new_preshape/SML/6x15.mat');
coils{idx_coil} = coordinates;
% idx_coil = idx_coil + 1;
% coordinates=load('./deployment_rob/new_preshape/LKC/4x10_final.txt');
% %load('./deployment_rob/new_preshape/SML/3x8.mat');
% coils{idx_coil} = coordinates;
% idx_coil = idx_coil + 1;
% coordinates=load('./deployment_rob/new_preshape/SML/3x8_final.txt');
% %load('./deployment_rob/new_preshape/SML/3x8.mat');
% coils{idx_coil} = coordinates;
% idx_coil = idx_coil + 1;
% coordinates=load('./deployment_rob/new_preshape/SML/2x4_final.txt');
% %load('./deployment_rob/new_preshape/SML/2x4.mat');
% coils{idx_coil} = coordinates;
% idx_coil = idx_coil + 1;
% coordinates=load('./deployment_rob/new_preshape/SML/2x2_final.txt');
% %load('./deployment_rob/new_preshape/SML/2x2.mat');
% coils{idx_coil} = coordinates;

% load IA Sac geometry
% load('WSN_Sac.mat');
%load('./deployment_rob/new_preshape/LKC_Sac.mat');
[fout,vout]=stlread('CJA_Final_sac.stl');
mid_vout=mean(vout);
% user input: starting point and starting direction
starting_point = [mid_vout;mid_vout;mid_vout;mid_vout];% moslty centriod of IA Sac
centroid=[211.3,337.2,399.5];% centriod of neckplane
% centroid=[0 0 0];
sp=centroid + (starting_point(1,:)-centroid)*1;
coordinates = coils{1};
% Starting Directions
R = eye(3); % Rotation matrix
b = centroid-starting_point(1,:);
b=-b;
a = coordinates(2,:)-coordinates(1,:);
v = cross(a,b)/norm(a)/norm(b);
s = norm(v);
c = a*b'/norm(a)/norm(b);
vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = (eye(3) + vx + vx^2/(1+c))*R;
R_pp=R;

angg=90;
x2=[cosd(angg) sind(angg) 0];
x1=[1 0 0];
a=x1;b=x2;
v = cross(a,b)/norm(a)/norm(b);
s = norm(v);
c = a*b'/norm(a)/norm(b);
vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
R = (eye(3) + vx + vx^2/(1+c))*R;

% set some constants
eps = 10^-6;
stepsize = 1;
radius = 0.17;%HEC 0.1524; WSN 0.1778; CJA: 0.15875
n_cnt = 300;
grid_n = 20;
v_min = min(vout);
v_max = max(vout);
coil_r = cell(No,1);
coil_h = cell(No,1); %helping table for previous coils
coil_c = cell(No,1); %count for helping table
bad_ratio = 0.05; % smooth the 0.05 points with large bending energy
bad_max = 0.4; % if more then 40% of points have large bending energy, converge
avg_len = 10; % average the bending energy with 20 neighbours
tic;
collision_cnt = 0;
for coil_idx = 1:No % Coils
    % coordinates = coordinates(1:1000,:);
    coordinates = coils{coil_idx};
    [node_cnt, dim] = size(coordinates);
    offsets = coordinates(2:stepsize:node_cnt,:) - coordinates(1:stepsize:node_cnt-1,:);
    meshpoint = zeros(size(offsets,1)+1, dim);
    grid_location = zeros(size(meshpoint));
    [help_t,cnt_t] = neighbour_table(fout);
    %% insert the coil inside the vessel step by step
%     meshpoint(1,:) = starting_point(coil_idx,:);
    meshpoint(1,:) = sp;
    grid_location(1,:) = update_grid(meshpoint(1,:),v_min,v_max,grid_n);
    idx = 1;
%     R = eye(3);
    rot=ones(node_cnt,1);
    cnt_pp=ones(node_cnt,1);
    ii=1;
    ip_size=1;
    flag=1;
%     R=R_pp;
 
    while ii < node_cnt % Segment
        ii=ii+1;
        idx = idx + 1;
        cc_cnt=1;
        
        meshpoint(idx,:) = meshpoint(idx-1,:) + offsets(idx-1,:)*(R'^rot(idx));    
        seg_pp=[meshpoint(idx,:);meshpoint(idx-1,:)];
        R_pp=R;
        
         flag=1;
         cnt_pp(idx)=cnt_pp(idx)+1;
         fprintf('idx = %i \t cnt= %i \n',idx,cnt_pp(idx));
      
        while flag < 300 % Collision Detection 
%             flag(idx) = flag(idx)+1;
            seg = [meshpoint(idx,:);meshpoint(idx-1,:)];
            [collision, tri_idx] = CollisionDetection(seg,vout,fout,radius);
            collision_cnt = collision_cnt + collision;
            
            if (collision > 0) % collision response       
            
                if idx>5
                    p = meshpoint(idx-2,:);
                else
                    p = starting_point(coil_idx,:);
                end
                cnt = 0;
                new_seg = seg;
                
                while(collision > 0 && cnt < 1000) % Rotate Segment until there is no collision with IA Sac
                    new_seg = SegmentRotation(new_seg,p,vout,fout,tri_idx,radius, eps);
                    tris = tri_neighbour(tri_idx,fout,n_cnt, help_t, cnt_t);
                    [collision, tri_idx] = CollisionDetection(new_seg,vout,fout(tris,:),radius);
%                     [collision, tri_idx] = CollisionDetection(new_seg,vout,fout(tris,:),radius);
                    if collision > 0
                        tri_idx = tris(tri_idx);
                    end
                    cnt = cnt + 1;
                end % End While-loop for Segment Rotation
                
                if cnt >= 1000
                    fprintf('\tRotation error in Segment idx= %i... \n\tMarching 2 points backward\n',idx);
%                     rot(idx-1)=rot(idx-1)+1;
%                     ii=ii-2;
%                     idx=idx-2;
%                     fprintf('\tidx = %i and ii=%i \n',idx,ii);
                    
                    %%%%%%%%
                    ii=ii-3;
                    idx=idx-3;
                    ip(ip_size,:)=meshpoint(idx,:);
                    ip_size=ip_size+1;    
                    seg11 = meshpoint(idx-1:idx,:); 
                    AngRot=5;
                    [meshpoint(idx,:)]=SegRot(meshpoint(idx,:),meshpoint(idx-1,:),meshpoint(idx-2,:),AngRot);
                    seg = [meshpoint(idx,:);meshpoint(idx-1,:)];
                    [collision, tri_idx] = CollisionDetection(seg,vout,fout,radius);
%                     [collision, tri_idx] = CollisionDetection(seg,vout,fout,radius);
                    if collision>=1
                        AngRot=-2*AngRot;
                        [meshpoint(idx,:)]=SegRot(meshpoint(idx,:),meshpoint(idx-1,:),meshpoint(idx-2,:),AngRot);
                    end
                    a = seg11(2,:)-seg11(1,:);
                    b = meshpoint(idx,:)-meshpoint(idx-1,:);

                    v = cross(a,b)/norm(a)/norm(b);
                    s = norm(v);
                    c = a*b'/norm(a)/norm(b);

                    vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                    R = (eye(3) + vx + vx^2/(1+c))*R;
                    %%%%%%%%%%
                    
                    break;                    
                else
                    if ~isreal(new_seg)
%                         rot=1;
                        meshpoint(idx,:) = seg(1,:);
                    else
%                         rot=1;
                        meshpoint(idx,:) = new_seg(1,:);
                    end
                    a = seg(1,:)-seg(2,:);
                    b = meshpoint(idx,:)-meshpoint(idx-1,:);

                    v = cross(a,b)/norm(a)/norm(b);
                    s = norm(v);
                    c = a*b'/norm(a)/norm(b);

                    vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                    R = (eye(3) + vx + vx^2/(1+c))*R;
                end
                    
                % compute the rotation matrix
%                 a = seg(1,:)-seg(2,:);
%                 b = meshpoint(idx,:)-meshpoint(idx-1,:);
% 
%                 v = cross(a,b)/norm(a)/norm(b);
%                 s = norm(v);
%                 c = a*b'/norm(a)/norm(b);
% 
%                 vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
%                 R = (eye(3) + vx + vx^2/(1+c))*R;
            end
            
            if cnt_pp(idx)>=3
                break;
            end

            %check collision with coil itself
            seg = meshpoint(idx-1:idx,:);  
            grid_location(idx,:) = update_grid(meshpoint(idx,:),v_min,v_max,grid_n);
            [collision,seg_idx] = collisionSeg(seg,meshpoint(1:idx-2,:),radius,...
                grid_location(idx-1:idx,:),grid_location(1:idx-2,:));
            if collision > 0 && idx>4
%                 meshpoint(idx,:) = meshpoint(idx,:)+0.1*rand(1,3);
%                 meshpoint(idx,:) = meshpoint(idx-1,:)+...
%                     (meshpoint(idx,:)-meshpoint(idx-1,:)).*norm(offsets(idx-1,:));
%                 meshpoint(idx,:) = meshpoint(idx-1,:) + offsets(idx-1,:)*(R'^rot(idx));
%                 rot(idx)=rot(idx)+1;
                  AngRot=-0.15;
%                  AngRot=0.15*flag*((-1)^(flag+1));
                [meshpoint(idx,:)]=SegRot(meshpoint(idx,:),meshpoint(idx-1,:),meshpoint(idx-2,:),AngRot);

                % compute the rotation matrix
                a = seg(2,:)-seg(1,:);
                b = meshpoint(idx,:)-meshpoint(idx-1,:);

                v = cross(a,b)/norm(a)/norm(b);
                s = norm(v);
                c = a*b'/norm(a)/norm(b);

                vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                R = (eye(3) + vx + vx^2/(1+c))*R;

            end
            c1=collision;
            
%             % check with itself
%             seg = meshpoint(idx-1:idx,:);
%             for jj = 1:coil_idx
%                 [collision,seg_idx] = collisionSeg(seg,coil_r{coil_idx},radius,...
%                     grid_location(idx-1:idx,:),grid_location(1:idx-2,:));
%                 if collision > 0
% %                     meshpoint(idx,:) = meshpoint(idx,:)+rand(1,3);
% %                     meshpoint(idx,:) = meshpoint(idx-1,:)+...
% %                         (meshpoint(idx,:)-meshpoint(idx-1,:)).*norm(offsets(idx,:));
%                       collision
%                       meshpoint(idx,:) = meshpoint(idx-1,:) + offsets(idx-1,:)*(R'^rot);
%                 end
%             end
            c2=collision;
            
            % check with previous coils
            if coil_idx > 1
                seg = meshpoint(idx-1:idx,:);
                for jj = 1:coil_idx-1
                    temp_p = coil_r{jj};%previous coil
                    temp_h = coil_h{jj};%helping table for previous coil
                    temp_c = coil_c{jj};%count for helping table
                    [collision, seg_idx] = collisionPreSeg(seg,temp_p,temp_h,temp_c,radius,v_min,v_max,grid_n);
                    if collision > 0 && idx>4
%                         meshpoint(idx,:) = meshpoint(idx,:)+rand(1,3);
%                         meshpoint(idx,:) = meshpoint(idx-1,:)+...
%                             (meshpoint(idx,:)-meshpoint(idx-1,:)).*norm(offsets(idx-1,:));
%                     meshpoint(idx,:) = meshpoint(idx-1,:) + offsets(idx-1,:)*(R'^(rot(idx)+1));
%                     rot(idx)=rot(idx)+1;
                         AngRot=-0.15;
%                          AngRot=0.15*flag*((-1)^(flag+1));
                        [meshpoint(idx,:)]=SegRot(meshpoint(idx,:),meshpoint(idx-1,:),meshpoint(idx-2,:),AngRot);

                        % compute the rotation matrix
                        a = seg(2,:)-seg(1,:);
                        b = meshpoint(idx,:)-meshpoint(idx-1,:);

                        v = cross(a,b)/norm(a)/norm(b);
                        s = norm(v);
                        c = a*b'/norm(a)/norm(b);

                        vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
                        R = (eye(3) + vx + vx^2/(1+c))*R;

                    end
                end
            end
            c3=collision;
            
            collision=c1+c2+c3;
            flag = flag+1;
            
            if collision==0                
                break;
            end
        end % End While-loop for Collision Detection
        
        if flag >= 300 && idx>5
            fprintf('\tCollision within Coil itself.... \n')
%             rot(idx)=rot(idx)+1;
            ii=ii-3;
            idx=idx-3;
            ip(ip_size,:)=meshpoint(idx,:);
            ip_size=ip_size+1;   
            seg11 = meshpoint(idx-1:idx,:); 
            AngRot=3;
            [meshpoint(idx,:)]=SegRot(meshpoint(idx,:),meshpoint(idx-1,:),meshpoint(idx-2,:),AngRot);
            seg = [meshpoint(idx,:);meshpoint(idx-1,:)];
            [collision, tri_idx] = CollisionDetection(seg,vout,fout,radius);
            if collision>=1
                AngRot=-2*AngRot;
                [meshpoint(idx,:)]=SegRot(meshpoint(idx,:),meshpoint(idx-1,:),meshpoint(idx-2,:),AngRot);
            end
            a = seg11(2,:)-seg11(1,:);
            b = meshpoint(idx,:)-meshpoint(idx-1,:);

            v = cross(a,b)/norm(a)/norm(b);
            s = norm(v);
            c = a*b'/norm(a)/norm(b);

            vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
            R = (eye(3) + vx + vx^2/(1+c))*R;
%             ii=ii-1;
%             idx=idx-1;
            %break;
        end
            grid_location(idx,:) = update_grid(meshpoint(idx,:),v_min,v_max,grid_n);
%             a = seg_pp(1,:)-seg_pp(2,:);
%             b = meshpoint(idx,:)-meshpoint(idx-1,:);
% 
%             v = cross(a,b)/norm(a)/norm(b);
%             s = norm(v);
%             c = a*b'/norm(a)/norm(b);
% 
%             vx = [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
%             R = (eye(3) + vx + vx^2/(1+c))*R_pp;

    end% for ii = 1+stepsize : stepsize : node_cnt

    
    coil_r{coil_idx} = [meshpoint,grid_location];
    [coil_h{coil_idx}, coil_c{coil_idx}]= previousCoils(meshpoint,grid_location);
    
   runningTime(coil_idx) = toc;
end % for coil_idx

%smoothing 
for coil_idx=1:No
    coil_r{coil_idx}=smoothing_pp(coil_r{coil_idx},coils{coil_idx});
end
runningTime(coil_idx+1) = toc;

% Plotting Figures
figure(1);
node_xyz=vout;
face_num=length(fout);
IV.vertices=vout;
IV.faces=fout;
patch(IV,'edgecolor','none','facecolor','w');
hold on
alpha 0.4;
axis off;
camlight;
hold on
view(10,0);
%nnn = idx;
color = ['b','r','g','k'];
color = repmat(color,1,ceil(No/4));
for ii = 1:No
    meshpoint = coil_r{ii};
    meshpoint = meshpoint(:,1:3);
%     meshpoint = meshpoint(1:idx,:);
    t1 = cscvn(meshpoint');
    pp = fnplt(t1);
    q = linspace(0,2*pi,30);
    base = 0.175*[cos(q);sin(q)];%HEC; WSN 0.1778; CJA: 0.15875
    traj = pp;
    [X,Y,Z] = extrude(base,traj);
    surf(X,Y,Z,'FaceColor',color(ii),'EdgeColor',color(ii),'FaceAlpha',1);
    % surf2stl('coil_test.stl',X,Y,Z,'ascii')
end
title(['Starting Point= ',num2str(starting_point(1,:))]);
% a=char(datetime);
% b='_';
% str1=[a(1:2),a(4:6),a(8:11),b,a(13:14),a(16:17)];
% fig1name=sprintf('./Results/%s_fig1',str1);
% savefig(fig1name);
% 
% figure(2)
% plot3(meshpoint(:,1),meshpoint(:,2),meshpoint(:,3));
% hold on
% plot3(ip(:,1),ip(:,2),ip(:,3),'.r');
% hold on
% fig2name=sprintf('./Results/%s_fig2',str1);
% savefig(fig2name);

% Save Output in .mat file
%coildata=sprintf('./Results/%s_CoilData',str1);
%save(coildata,'starting_point','fout','vout','coil_r','coils','runningTime');
