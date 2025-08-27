%%{
clear all;clc;


px = 136+0;
pz = 236+0;

py = 170*3+50;

ite = 1;

'Start!!'

pathA = '/mnt/scratch/hcy/SBMES_2025/2025_0818_GraMS/DATs/'
%ext = '.dat'; 
ext = '.txt'; 

for N = 0: 0 % : 6

	frt = num2str(1000+N)
	fln = [pathA 'LSS_3D_OR_136x560x236_Z_II_' frt '.dat']

	fid = fopen(fln,'r');
	skip1 = fread(fid,1,'int32');
	% BD = fread(fid,inf,'double');
	BD = fread(fid,inf,'single');
	fclose(fid);
	
	[r1 c1] = size(BD)	
	dp = reshape(BD(1:r1-1),px,py,pz);
	
	dp2 = imresize3(dp,2);
	dp2 = smooth3(dp2,'box',3);	
	size(dp2)
	px*py*pz
	


    frm = num2str(100+N);
    fnm = ['T' frm]
	
	fid = fopen([pathA 'GrMic_II_LO_272x1120x472_' fnm '.txt'],'w');
	fprintf(fid,'%e\t',dp2);
	fclose(fid);


    %% domain parameters
	fname = [pathA 'GRP_LS_136x560x236_Z_II_' frm ext]   
	
	%   fid = fopen(fname,'r');
	%   skip1 = fread(fid,1,'int32');
	%   % BD = fread(fid,inf,'double');
	%   BD = fread(fid,inf,'single');
	%   fclose(fid);
	   
	fid = fopen(fname,'r');
	BD = fscanf(fid,'%g');
	fclose(fid);  

   [r1 c1] = size(BD)
   
   % px*py*pz
	dp = reshape(BD(1:r1-0),px,py,pz);
	%  dp = reshape(BD(1:r1-1),px,py,pz);
   
   % dp = smooth3(dp,'box',3);
   % dp = 0.5*(1+tanh(dp/1));
   
%   dp = dp(1:px,1:py,1:pz);
	dp2 = imresize3(dp,2);
	dp2 = smooth3(dp2,'box',3);
	
	size(dp2)
	clear dp 
	dp = dp2;
	
	fid = fopen([pathA 'GrMic_II_DS_272x1120x472_' fnm '.txt'],'w');
	fprintf(fid,'%e\t',dp);
	fclose(fid);


   %%
%   fname = [pathA 'POR_LS_114x390x198_Z_II_' frm ext] 
  
%    fid = fopen(fname,'r');
%    BD = fscanf(fid,'%g');
%    skip1 = fread(fid,1,'int32');
%    % % BD = fread(fid,inf,'double');
%    BD = fread(fid,inf,'single');
%    fclose(fid);
   

	fname = [pathA 'POR_LS_136x560x236_Z_II_' frm '.txt']     
	fid = fopen(fname,'r');
	BD = fscanf(fid,'%g');
	fclose(fid);  
	
	[r1 c1] = size(BD)
	
	px*py*pz
	dp = reshape(BD(1:r1-0),px,py,pz);
%  dp = reshape(BD(1:r1-1),px,py,pz);
  
  % dp = smooth3(dp,'box',3);
  % dp = 0.5*(1+tanh(dp/1));
  
%   dp = dp(1:px,1:py,1:pz);
	dp2 = imresize3(dp,2);
	dp2 = smooth3(dp2,'box',3);
	
	
	size(dp2)
	clear dp 
	dp = dp2;
	
	fid = fopen([pathA 'GrMic_II_DL_272x1120x472_' fnm '.txt'],'w');
	fprintf(fid,'%e\t',dp);
	fclose(fid);

    %%
    
    % px = 360;
    % py = 440;
    % pz = 320;
    
    %% grid system
    % dh = 1 *0.5;
    % dh = 1;
    % 
    % Xpp(1) = 0;
    % for i = 2:px
    %     Xpp(i) = Xpp(i-1)+dh;
    % end
    % 
    % Ypp(1) = 0;
    % for j = 2:py
    %     Ypp(j) = Ypp(j-1)+dh;
    % end
    % 
    % Zpp(1) = 0.0d0;
    % for k = 2:pz
    %     Zpp(k)= Zpp(k-1)+dh;
    % end
    % 
    % [Ypp(1) Ypp(py) Xpp(1) Xpp(px) Zpp(1) Zpp(pz)]
    % 
    % % plot(Ypp,'o')
    % 
    % [X,Y,Z] = meshgrid([Ypp],[Xpp],[Zpp]);
    % size(X)
    % 
    % 
    % NewX =  X(1:ite:px,1:ite:py,1:ite:pz);
    % NewY =  Y(1:ite:px,1:ite:py,1:ite:pz);
    % NewZ =  Z(1:ite:px,1:ite:py,1:ite:pz);
    % NewP = dp(1:ite:px,1:ite:py,1:ite:pz);
    % 
    % clear X Y Z dp
    % 
    % dp = NewP;


    size(dp)
    
    % figure()
    % 
    % 
    % isovL = 0.0;
    % 
    % p1 = patch(isosurface(NewX,NewY,NewZ,dp,isovL));
    % p2 = patch(isocaps(NewX,NewY,NewZ,dp,isovL));
    % 
    % isocolors(NewX,NewY,NewZ,dp,p1)
    % set(p1,'FaceColor','interp','EdgeColor','none')
    % 
    % isonormals(NewX,NewY,NewZ,dp,p1)
    % 
    % isocolors(NewX,NewY,NewZ,dp,p2)
    % set(p2,'FaceColor','interp','EdgeColor','none')
    % 
    % 
    % set(gca,'Projection','perspective')
    % view(3); daspect([1,1,1]);
    % % axis([Ypp(1) Ypp(py) Xpp(1) Xpp(px) Zpp(1) Zpp(pz)]);
    % box on;
    % lightangle(45,60)
    % lighting phong   
    
    
%    pln(:,:) = dp(:,:,270);
%     pln(:,:) = dp(:,:,472);
    
%    figure()
 %   surf(pln)
    
%    shading interp
end

clear all;




