clc;clear;close all;
% Defining the number of grids in rows and columns
N=7;                    % Number N defines the size of the rectangular grid here it is 7*7
NP=((N-1)*5)^2;         % Gives the number of pixels in the grid for pixel width of 0.2
% generates the 16 number of nodes
x=linspace(1,N,N);
y=linspace(1,N,N);
l1=zeros(N,2);
l2=zeros(N-2,2);
l3=zeros(N,2);
l4=zeros(N-2,2);
% L stores the coorinates of nodes
L=zeros(4*N,2);
% v_d the distance from focus to vertex that controls the width of ellipse
v_d=0.01;
beta=0.75; 
% alpha=6.75; 
mu=2.35;
for i=1:N
l1(i,:)=[x(1),y(i)];
l3(i,:)=[x(N),y(N-i+1)];
end
for j=1:N-2
    l2(j,:)=[x(N-j),y(1)];
    l4(j,:)=[x(j+1),y(N)];
end
L=[l1;l4;l3;l2];
plot(L(:,1),L(:,2),'or','MarkerSize',12,'MarkerFaceColor',[0.5,0.5,0.5]);
rectangle('Position',[1 1 N-1 N-1],'LineWidth',2)
xlabel('Distance X(meters)');
            ylabel('Distance Y(meters)');
            title('True SLF with sensor node represented on boarder of rectangular grid')
hold on;
% xlim([0, N+1]);
% ylim([0, N+1]);
xlim([1, N]);
ylim([1, N]);
grid on;
% -------- creating the object -------- %
% generalized square object creation %
ov=[2 3 1 1];
o1=[ov(1) ov(1) ov(1)+ov(3)  ov(1)+ov(3)];o2=[ov(2)  ov(2)+ov(4)  ov(2)+ov(4) ov(2)];
i=1;
for r=1:ov(3)
    for c=1:ov(4)
    pp(i,:)=[((ov(1)+r-1)+(ov(1)+r))/2,((ov(2)+c-1)+(ov(2)+c))/2];
    i=i+1;
    end
end
% creating the number of links eliminating the same point as tx and rx  %
p=0;
for i=1:length(L)
    for k=1:length(L)
    if L(k,1)== L(i,1)&&L(k,2)==L(i,2)
        disp('No line is possible to draw');
    else
        xx=[L(i,1),L(k,1)];yy=[L(i,2),L(k,2)];
        line(xx,yy);    
        p=p+1;
        d=L(i,:)-L(k,:);
        lin(p).tx=L(i,:);
        lin(p).rx=L(k,:);
        lin(p).dist=sqrt(d*d'); 
%         if lin(p).att==1
        lin(p).c=lin(p).dist/2;
        lin(p).Elips_ma=v_d+lin(p).c;              % Gives ellipse majoraxis
        lin(p).Elips_mi=sqrt((lin(p).Elips_ma)^2-lin(p).c^2); % Gives ellipse minoraxis
        lin(p).cn=[(xx(1)+xx(2))/2,(yy(1)+yy(2))/2];    % Gives ellipse center
        lin(p).ang=180-((atan((yy(2)-yy(1))/(xx(2)-xx(1))))*180/pi);  % Gives ellipse angle
        [XX,YY] = calculateEllipse(lin(p).cn(1),lin(p).cn(2),lin(p).Elips_ma,lin(p).Elips_mi,lin(p).ang);
       
        plot(XX, YY, 'g','Linewidth',.5);
plot(L(:,1),L(:,2),'or','MarkerSize',12,'MarkerFaceColor',[0.5,0.5,0.5]);
        hold on;   
     end
    end
end

%  Assigning Pixel coordinates for the square grid  %
cn=1;
for i=1:length(lin)
    lin(i).pix=[];
    for m=1:.2:N-0.2
        for n=1:.2:N-0.2
        pix=[n+.1,m+.1];
        lin(i).pix=[lin(i).pix;pix]; 
        end
    end
end 

% steps to obtain the error covariance matrix
% finding number of pixels in SLF
m=1;
k=1;
pix=zeros(NP,2);
for i=1:sqrt(NP)
    for j=1:sqrt(NP)
        pix(m,:)=[k+0.2,j+0.2];
        m=m+1;
    end
    k=k+1;
end
dist_pix=[];
pt=1;
for k=1:NP
    for j=1:NP   
         dist_pix(pt,:)=norm(pix(k,:)-pix(j,:));
%        dist_pix(pp,:)=norm(lin(i).pix(k,:)-lin(i).pix(j,:));
        pt=pt+1;      
    end    
end
% Error covariance that relates the SLF with transreceiver pairs

sigma=0.4;                               % pixel variance , reducing this will smooth more and increasing creates sharp peak
delta=2.05;                               %Pixel correlation constant, increasing this will blur the image more 
fact=sigma/delta;
pixin=fact*exp(-(dist_pix)/delta);
pixmat=vec2mat(pixin,NP);
pixm_inv=pinv(sqrtm(pixmat));

% attenuation calculation for Tx, Rx for synthetic loss feild values  %

TP=30;
for i=1:length(lin)    
    t=lin(i).tx;
    r=lin(i).rx;
    distance=polyintersect(t,r,o1,o2);
    lin(i).distance=distance;
    lin(i).TP=30;
    RP=TP-lin(i).distance;
    lin(i).rss=RP;
    if  distance>0
    lin(i).att= 1;
    else
        lin(i).att=0;
    end
end

% Calculate ellipse for the loss feilds only  %  
for i=1:length(lin)
%     for i=24
      if lin(i).att==1
            [XX1,YY1]=calculateEllipse(lin(i).cn(1),lin(i).cn(2),lin(i).Elips_ma,lin(i).Elips_mi,lin(i).ang);
%             lin(i).atenElips=
            figure(2)
%             plot(L(:,1),L(:,2),'or',X1, Y1, 'g','Linewidth',0.5);
            plot(XX1, YY1, 'r','Linewidth',0.5);
            hold on;
            rectangle('Position',ov,'LineWidth',1,'FaceColor','b')
            grid on;
            xlabel('Distance X(meters)');
            ylabel('Distance Y(meters)');
            title('ellipse fit for transreceivers that are affected by an object')
            
      end
end

% weighing the attenuated pixels inside the ellipse for the loss fields %
lin(i).gap=[];
for i=1:length(lin)
     T=lin(i).tx;
     R=lin(i).rx;
     lin(i).pixin=[];
     lin(i).pixid=[];
%       lin(i).pixid1=[];
     if lin(i).att==1         
      for j=1:length(lin(i).pix)    
        gap= norm(lin(i).pix(j,:)-T)+norm(lin(i).pix(j,:)-R);
        lin(i).gap=[lin(i).gap,gap];
         if gap <=2*lin(i).Elips_ma            
            pixin=lin(i).pix(j,:);
            lin(i).pixin=[lin(i).pixin;pixin];
        end
      end
     end
end

% Assigning zeros to make all weight vectors of length 1*16  %
nul=zeros(NP,2);
% w=zeros(1,16);
        for i=1:length(lin)
            w=zeros(1,NP);
            ad=zeros(1,NP);
            [row,col]=size(lin(i).pix);
            for r=1:row
                [row1,col1]=size(lin(i).pixin);
                for c=1:row1
                    if isempty(lin(i).pixin)
                        lin(i).pixin=nul;
                    else
                     tmp=find((lin(i).pix(r,1)==lin(i).pixin(c,1)) && (lin(i).pix(r,2)==lin(i).pixin(c,2)));
                      if ~isempty(tmp)                         
                      lin(i).pixid=[lin(i).pixid,r];
                      lin(i).w(lin(i).pixid)=1/lin(i).dist;
                      end
                    end
                end
            end
        end       
ad=zeros(1,NP);
 W_M=zeros(length(lin),NP);
 w1=zeros(1,NP);
 
 % Obtaing the weight matrix for all 552 links %
 for i=1:length(lin)
     if isempty(lin(i).pixid)
         W_M(i,:)=zeros(1,NP);
     else
         lin(i).w1=zeros(1,NP);
         lin(i).w1(1,1:length(lin(i).w))=lin(i).w;
         W_M(i,:)=lin(i).w1;       
     end
 end 
% Generating the tikhonov matrix using identity matrix %
% Q=eye(16);
Q=eye(NP);
QT=transpose(Q);
del_y=zeros(552,1);
% cn=1;
for i=1:length(lin)
    del_y(i,1)=lin(i).distance;
end
% Adding noise which is combination of fading and measurement noise %
% Ts=0.01; %Given sampling period
Ts=0.1;
SNRdb = 10; %Given SNRdb
% variance = 1/(Ts*(10^(SNRdb/10)));  
% variance = 1/(Ts*(10^(SNRdb/10))); 
no= sqrt(1).*randn(1,size(del_y,1));
% no= sqrt(variance).*randn(1,size(del_y,1));
noise=0.1*(no');
del_y1=del_y+noise; % loss feild with noise used for estimating the imaging vector %
% RSS with quantization
% Use quantize to quantize data to a fixed-point type with a wordlength of 
% 3 bits, a fraction length of 2 bits, convergent rounding, and wrap on overflow.
q = quantizer('fixed','convergent','wrap',[8 2]);
% x = (-2:eps(q)/4:2)';
dely = del_y';
% dely= (del_y(1):eps(q)/4:del_y(552))';
% r=range(q);
yq = quantize(q,dely);
yq=abs(yq)';

W_MT= transpose(W_M); 
% alpha=0.75;         % Lower alpha fits the data
%  alpha=2;          % Higher alpha matches prior information
% W_Minv=pinv(W_MT*W_M);
W_Minv=(W_MT*W_M);
T_reg=beta*(QT*Q);
W_reg=pinv(W_Minv+T_reg);
X_LS=(W_reg*W_MT)*yq;
I_TM=vec2mat(X_LS,30);
I_TM=I_TM/max(max(I_TM));
% figure(3);
% % x1=[1.5 4.5];
% % y1=[1.5 4.5];
% x1=[1 7];
% y1=[1 7];
% imagesc(x1,y1,I_TM);
% set(gca,'YDir','normal');
%             xlabel('Distance X(meters)');
%             ylabel('Distance Y(meters)');
% title('Image reconstruction using tikhonov identity matrix');
% 
% figure(33);
% colormap(gray);
% x1=[1 7];
% y1=[1 7];
% imagesc(x1,y1,I_TM);
% set(gca,'YDir','normal');
%             xlabel('Distance X(meters)');
%             ylabel('Distance Y(meters)');
% title('Tikhonov regularization using Identity matrix in gray scale');
% 
% figure(34);
% surf(I_TM);
% title('surface');
% 
% % % Reconstruction using noisy data %
% W_Minv1=(W_MT*W_M);
% T_reg1=beta*(QT*Q);
% W_reg1=pinv(W_Minv1+T_reg1);
% X_LSN=(W_reg1*W_MT)*yq;
% I_TMN=vec2mat(X_LSN,30);
% figure(35);
% % x2=[1.5 4.5];
% % y2=[1.5 4.5];
% x2=[1 7];
% y2=[1 7];
% imagesc(x2,y2,I_TMN);
% set(gca,'YDir','normal');
% title('Image reconstruction using tikhonov identity matrix along with 0.1% noise');

% generation of error covariance matrix
% pixm_inv=pinv(pixmat);
W_MT= transpose(W_M);
W_MF=(W_MT*W_M); 
W_regco=pinv(W_MF+0.75*pixm_inv);
W_regco1=0.075*pinv(pixm_inv'*pixm_inv);
% W_regco1=pinv(pixm_inv'*pixm_inv)/0.75;

% Reconstruction using error covariance matrix
% pixm_inv=pinv(pixmat);
W_MT= transpose(W_M);
W_MF=(W_MT*W_M); 
W_regco=pinv(W_MF+10.75*pixm_inv);
% W_regco=pinv(W_MF+5.75*pixm_inv'*pixm_inv);
I_tc=(W_regco*W_MT)*yq; 
I_TMc=vec2mat(I_tc,30);
I_TMc=I_TMc/max(max(I_TMc));
figure(36);
x1=[1 7];
y1=[1 7];
imagesc(x1,y1,I_TMc);
imagesc(x1,y1,I_TMc);
set(gca,'YDir','normal');
            xlabel('Distance X(meters)');
            ylabel('Distance Y(meters)');
title('Image reconstruction with quantization effects');

%  Linear eps-insensitive loss, nu-version with nu=0.5 and C=1.0, primal solution
X=W_M;
y=yq;
ell=size(X,1);
 n=size(X,2);
 lambda = 10;
%  C = 1/500;      % Initial assumption
 C = 1;
 epsilon = 0.5;
 H = diag([ones(1,n) 0 C*zeros(1,2*ell)]);
 f = [zeros(n+1,1);C*ones(2*ell,1)];
 A = [X ones(ell,1) zeros(ell) -eye(ell);
-X -ones(ell,1) -eye(ell) zeros(ell)];
 b = [epsilon + y; epsilon - y];
 Aeq = zeros(1,n + 1 + 2*ell);
 beq = 0;
 LB = [-inf*ones(n+1,1);zeros(2*ell,1)];
 wbxipxim = quadprog(H,f,A,b,Aeq,beq,LB);
w = wbxipxim(1:n);
 b = wbxipxim(n+1);
 yest = X*w + b;
% Linear epsilon-insensitive loss with epsilon=0.5 and C=0.1, dual solution
G=X*X';
 H = [G -G; -G G];
 f = [epsilon - y; epsilon + y];
 A = zeros(1,2*ell);
 b = 0;
 Aeq = [ones(1,ell) -ones(1,ell)]; beq = 0;
 LB = zeros(2*ell,1);
 UB = C*ones(2*ell,1);
 alphapalpham = quadprog(H,f,A,b,Aeq,beq,LB,UB);
alpha = alphapalpham(1:ell)-alphapalpham(ell+1:2*ell);
i = find(alpha>0.0000001 & alpha<C-0.0000001);
b = y(i)-G(i,:)*alpha-epsilon;
b_f = mean(b);
yest = G*alpha + b_f;
wt = X'*alpha;


% W_T=W_regco*wt;
W_T=W_regco1*wt;
% W_T=wt;
W_Tf=vec2mat(W_T,30);
W_Tf=W_Tf/max(max(W_Tf));
figure(80);
x1=[1 7];
y1=[1 7];
imagesc(x1,y1,W_Tf);
imagesc(x1,y1,W_Tf);
set(gca,'YDir','normal');
xlabel('Distance X(meters)');
ylabel('Distance Y(meters)');
title('Image reconstruction using epsilon-SVR');

