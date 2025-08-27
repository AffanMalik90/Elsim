%%{vv
clear all;clc;

px = 186+0; py = 1460+0; pz = 320+0;


ite = 2;

fb = 111 ;
fr = 131 ;
%fr = fb  ;
fe = 135;

prx = 'BQ' ;
ser = '112' ;

it2 = 10;

px*py*pz

'Start!!'

pathA = '/mnt/scratch/hcy/SBMES_2022/2022_0720_A/DataCHA/';
pathB = '';
pathC = '/mnt/scratch/hcy/SBMES_2022/2022_1026_A/FIGs/';


TpNm = 'CP';
preG = [TpNm 'ASMB'];

%TpNm = 'CC';
%preC = [TpNm 'ASMB'];


ext = '.dat'; 
% ext = '.txt'; 

crss = 1;
isovL = 0.45;


%% grid system
%dh = 0.1;
dh = 0.5*0.325;

Xpp(1) = 0;
for i = 2:px
    Xpp(i) = Xpp(i-1)+dh;
end

Ypp(1) = 0;
for j = 2:py
    Ypp(j) = Ypp(j-1)+dh;
end


Zpp(1) = 0.0d0;
for k = 2:pz
    Zpp(k)= Zpp(k-1)+dh;
end

[Ypp(1) Ypp(py) Xpp(1) Xpp(px) Zpp(1) Zpp(pz)]

% plot(Ypp,'o')

[X,Y,Z] = meshgrid([Ypp],[Xpp],[Zpp]);
size(X)
%% domain parameters

fname = [pathA 'PSASMB' num2str(1000+1) ext]

fid = fopen(fname);
skip1 = fread(fid,1,'int32');
BD = fread(fid,inf,'single');
%BD = fread(fid,inf,'double');
fclose(fid);


[r1 c1] = size(BD)
dp = reshape(BD(1:r1-1),px,py,pz);
%dp = reshape(BD(1:r1-0),px,py,pz);

%NewX =  X(px/2:ite:px,1:ite:py,1:ite:pz);
%NewY =  Y(px/2:ite:px,1:ite:py,1:ite:pz);
%NewZ =  Z(px/2:ite:px,1:ite:py,1:ite:pz);
%NewP = dp(px/2:ite:px,1:ite:py,1:ite:pz);

NewX =  X(1:ite:px,1:ite:py,1:ite:pz);
NewY =  Y(1:ite:px,1:ite:py,1:ite:pz);
NewZ =  Z(1:ite:px,1:ite:py,1:ite:pz);
NewP = dp(1:ite:px,1:ite:py,1:ite:pz);

clear X Y Z dp

dpg = NewP;

psi = dpg;

%% fraction
fname = [pathB 'G' prx ser 'TmFx.txt']
fid = fopen(fname);
BD = fscanf(fid,'%g');
fclose(fid); 

SZ = size(BD);

cl = SZ(1)/7

Rslt = reshape(BD,7,cl);
    

X = Rslt(4,:);
CV = -0.1 - Rslt(7,:);
% % CR = Rslt(5,:);

FF = [[fb:it2:fr] fe];
SZ = size(FF);

for k = 1: SZ(2)
i = FF(k);

%% concentration  
%for i = fb: it2: fr

    fname = [pathA preG num2str(1000+i) ext]
%     
    fid = fopen(fname);
    skip1 = fread(fid,1,'int32');
   BD = fread(fid,inf,'single');
%   BD = fread(fid,inf,'double');
    fclose(fid);
    
    [r1 c1] = size(BD)

	CnG = reshape(BD(1:r1-1),px,py,pz);
%   CnG = reshape(BD(1:r1-0),px,py,pz);

    
   
    h1 = figure(1)

    nCn = CnG(1:ite:px,1:ite:py,1:ite:pz);
%	nCn = CnG(px/2:ite:px,1:ite:py,1:ite:pz);
    size(nCn)
   
    upb = max(max(max(CnG)))
    lwb = min(min(min(CnG)))
    
    p1 = patch(isosurface(NewX,NewY,NewZ,psi,isovL));
    p2 = patch(isocaps(NewX,NewY,NewZ,psi,isovL));

    isocolors(NewX,NewY,NewZ,nCn,p1)
    set(p1,'FaceColor','interp','EdgeColor','none')
    
    isonormals(NewX,NewY,NewZ,dpg,p1)

    isocolors(NewX,NewY,NewZ,nCn,p2)
    set(p2,'FaceColor','interp','EdgeColor','none')
    
  
    

    set(gca,'Projection','perspective')
    view(108,27); daspect([1,1,1]);
    axis([Ypp(1) Ypp(py) Xpp(1) Xpp(px) Zpp(1) Zpp(pz)]);
    box on;
    lightangle(108,27)
    lighting phong   

% %     
    whitebg('w')
    set(gca,'FontSize',48,'FontWeight','bold')
    
%     set(gca,'XTick',[0:4:8])
%     set(gca,'YTick',[0:4:16])
%     set(gca,'ZTick',[0:4:16])
    
%     set(gca,'XTick',[])
%     set(gca,'YTick',[])
%     set(gca,'ZTick',[])

%     
  caxis([0 1])

    colormap(parula)
    colorbar('Position',[0.827 0.11 0.011 0.815])
%     colorbar('YLim',[lwb upb],'YTick',[0.0 0.2 0.4 0.6 0.8 1.0], ...
%         'YTickLabel',{'0.0' '0.2' '0.4' '0.6' '0.8' '1.0'}, ...
%         'FontSize',28,'FontWeight','bold','Position',[0.82 0.1 0.015 0.8])
    
%     colorbar('YLim',[-0.067 1.00],'YTick',[0.0 0.2 0.4 0.6 0.8 1.0], ...
%         'YTickLabel',{'0.0' '0.2' '0.4' '0.6' '0.8' '1.0'}, ...
%         'FontSize',40,'FontWeight','bold','Position',[0.82 0.1 0.015 0.8])    

    set(gcf,'color','w')  
    set(gcf,'Position',[740 85 1460 1260])
    set(gca,'BoxStyle','full')
    
    title(['X = ' num2str(X(i)) ',   ' 'CV = ' num2str(CV(i)) 'V'])
    

    pause(0.5)
    hgsave(h1,[pathC TpNm '_' prx ser '_3D_' num2str(1000+i) '.fig'],'-v7.3')

    pause(1)

    close(h1)

end

clear all;






