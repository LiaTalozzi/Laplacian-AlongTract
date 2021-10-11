%%


%created by Dr. Lia Talozzi
%please cite Talozzi L et al. Along-tract analysis of the arcuate fasciculus using the Laplacian operator to evaluate different tractography methods. Magn Reson Imaging. 2018.
%%

function [final_stats]=Laplacian_segmentation (tract, ID, side, flag_smooth,flag_iper_smooth )
if exist('libraries', 'dir')
    addpath(genpath('libraries')) %add manually from the matlab interface if authomatically doesn't work
end

if ~exist('flag_smooth','var')
    flag_smooth=0;
end

if ~exist('flag_iper_smooth','var')

    flag_iper_smooth=0;
end

dim_carattere=40;  %dim carattere nei plot

home_dir_orig=pwd;

if contains(tract,'GCC') || contains(tract,'SCC') || contains(tract,'Fornix')
cd(tract)
else
cd([tract '_' side])
end

home_dir=pwd;
disp(home_dir)

 cd(home_dir)

  dir_fig=[home_dir '/figure'];
  dir_data=[home_dir '/' ID '/data'];

%%
%compute laplacian

%%
%caratteristiche singoli tratti:

if contains(tract,'GCC') || contains(tract,'SCC') || contains(tract,'Fornix')
  name_tract=[tract '_' th '_ontoMNI'];
else
name_tract=[tract '_' side '_' th '_ontoMNI'];
end

if contains(tract,'CST')
z_min_MNI=20;
z_max_MNI=70;
num_percentili=19;  %inizio a contare da 0
x_label_text='Inferior                                                     Superior';
end

if contains(tract,'AF')
z_min_MNI=40;
y_max_MNI=65;
num_percentili=14; %inizio a contare da 0
x_label_text='Frontal                                                      Temporal';
end

if contains(tract,'ATR')
y_min_MNI=57;
y_max_MNI=75;
x_min_MNI=30;
x_max_MNI=60;
num_percentili=3; %inizio a contare da 0
x_label_text='Frontal                                                      Talamus';
end

if contains(tract,'GCC')
y_max_MNI=87;
num_percentili=5; %inizio a contare da 0
x_label_text='Right                                                      Left';
end

if contains(tract,'SCC')
y_min_MNI=24;
num_percentili=4; %inizio a contare da 0
x_label_text='Right                                                      Left';
end

if contains(tract,'FAT')
y_max_MNI=85;
num_percentili=14; %inizio a contare da 0
x_label_text='Anterior                                                    Posterior';
end

if contains(tract,'OR')
y_max_MNI=40;
%y_max_MNI=75;
num_percentili=14;  %inizio a contare da 0
x_label_text='Anterior                                                    Posterior';
end

if contains(tract,'IFOF')
y_min_MNI=27;
y_max_MNI=78;
num_percentili=4;  %inizio a contare da 0
x_label_text='Frontal                                                     Occipital';
end

if contains(tract,'ILF')
y_min_MNI=35;
y_max_MNI=66;
num_percentili=5;  %inizio a contare da 0
x_label_text='Frontal                                                     Occipital';
end

if contains(tract,'UF')
%y_min_MNI=35;
y_max_MNI=80;
num_percentili=2;  %inizio a contare da 0
x_label_text='Frontal                                                     Temporal';
end

if contains(tract,'CGC')
%y_max_MNI=75;
z_min_MNI=39;
num_percentili=2;  %inizio a contare da 0
x_label_text='Frontal                                                     Occipital';
end


    subject=ID;

cd(home_dir)
cd(subject)
pwd

vol_tract_orig=niftiread(name_tract);
info_nifti=niftiinfo(name_tract);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%tract voxel coordinates extractions

vol_tract=vol_tract_orig>0;

region_stat=regionprops3(vol_tract, 'all');
[n_voxel position]=max(region_stat.Volume);
clear coord
coord_=region_stat.VoxelList{position,1};
coord(:,1)=coord_(:,2); coord(:,2)=coord_(:,1); coord(:,3)=coord_(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (x) coord(i,1) (y) coord(i,2) (z) coord(i,3)

tract_cut=zeros(size(vol_tract));
coord_cut=[];

if contains(tract,'CST')
for i=1:size(coord,1)
    if (coord(i,3) > z_min_MNI && coord(i,3) < z_max_MNI )
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end

if contains(tract,'CGC')
for i=1:size(coord,1)
    if (coord(i,3) > z_min_MNI)
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end

if contains(tract,'GCC')
for i=1:size(coord,1)
    if (coord(i,2) < y_max_MNI)
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end

if contains(tract,'SCC')
for i=1:size(coord,1)
    if (coord(i,2) > y_min_MNI)
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end


if contains(tract,'AF')
for i=1:size(coord,1)
    if (coord(i,2) < y_max_MNI && coord(i,3) > z_min_MNI )
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end

if contains(tract,'ATR')
for i=1:size(coord,1)
    if (coord(i,2) > y_min_MNI && coord(i,2) < y_max_MNI && coord(i,1) > x_min_MNI && coord(i,1) < x_max_MNI )
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end

if contains(tract,'FAT')
for i=1:size(coord,1)
    if (coord(i,2) < y_max_MNI )
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end

if contains(tract,'OR')
for i=1:size(coord,1)
    if (coord(i,2) < y_max_MNI)   %y=60 z=40
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end


if contains(tract,'IFOF')
for i=1:size(coord,1)
    if (coord(i,2) < y_max_MNI && coord(i,2) > y_min_MNI )   %y=60 z=40
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end

if contains(tract,'ILF')
for i=1:size(coord,1)
    if (coord(i,2) < y_max_MNI && coord(i,2) > y_min_MNI )
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end


if contains(tract,'UF')
for i=1:size(coord,1)
    if coord(i,2) < y_max_MNI
       coord_cut=[coord_cut; coord(i,:)];
        tract_cut(coord(i,1),coord(i,2),coord(i,3)) =1;
    end
end
end

coord= coord_cut;

vol=tract_cut;

%niftiwrite(single(tract_cut), ['core_' name_tract], info_nifti);

%%
%modeling mesh

%vol = bwareaopen(vol,50);  %chiedo che i volumi di pixel connessi minori di 50 vengano tolti
%vol = smooth3(vol,'gaussian',3);%faccio uno smoothing del rumore per togliere il rumore e meglio parametrizzare la superficie
 %vol=vol_tract;

region_stat=regionprops3(vol, 'VoxelList');

if size(region_stat.VoxelList,1) > 1
coord_big_delete=region_stat.VoxelList{2,1};
for i=1:size(coord_big_delete,1)
vol(coord_big_delete(i,2),coord_big_delete(i,1),coord_big_delete(i,3))=0;
end
end

vol = bwareaopen(vol,50,6);
%vol = bwareaopen(vol,50);  %chiedo che i volumi di pixel connessi minori di 30 vengano tolti
vol = smooth3(vol,'gaussian',3);

if flag_smooth == '1'
vol = smooth3(vol,'gaussian',3);%faccio uno smoothing del rumore per togliere il rumore e meglio parametrizzare la superficie
end

if flag_iper_smooth == '1'
vol = smooth3(vol,'gaussian',3);%faccio uno smoothing del rumore per togliere il rumore e meglio parametrizzare la superficie
end

s=isosurface(vol,0);
 %serve per fare un rendering di miglior qualit� e non far vedere i vertici

faces_orig=s.faces';
vertices_orig=s.vertices';
faces=faces_orig;
vertices=vertices_orig;
faces(1,:)=faces_orig(2,:); faces(2,:)=faces_orig(1,:); faces(3,:)=faces_orig(3,:);
vertices(1,:)=vertices_orig(2,:); vertices(2,:)=vertices_orig(1,:); vertices(3,:)=vertices_orig(3,:);


options.symmetrize = 1;
options.normalize = 0;
L0 = compute_mesh_laplacian(vertices,faces,'combinatorial',options);
%Performing eigendecomposition
%if size(vertices,2)<1500
    [U,D] = eig(full(L0));
%else
%    opts.disp = 0;
%    [U,D] = eigs(L0,50,'SM',opts);
 %   U = real(U(:,end:-1:1));
%end

L2 = U(:,2);

d=L2+abs(min(L2));  %in L2 ci sono sia valori negativi che positivi, cos� facendo li slitto  tutti sopra lo zero, assegnado il valore zero a negativo pi� basso

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     number of segments

n_seg=num_percentili-1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eigenvalues shifting to positive only

intervallo=(max(d)-min(d));  %definisco l'intervallo della nuova scala di valori positivi
 for i=1:size(d,1)
 d(i,1)= round((d(i,1)*(num_percentili))/intervallo) ;    %ragruppo in 20 valori
 end
 d=d';

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 if contains(tract,'CST') || contains(tract,'FAT')
column_min_inf=find(vertices(3,:) == min(vertices(3,:))  );
 if (round(mean(d(1, column_min_inf))) > (num_percentili/3) )
     d=-(d-num_percentili*ones(size(d,1),size(d,2)));
    %impongo che i percentili inizino dalla parte inferiore con valore 0 fino
 %alla parte superiore con valore 1
 end
 end


  if contains(tract,'AF') || contains(tract,'IFOF') || contains(tract,'OR') || contains(tract,'ILF') ||  contains(tract,'CGC') || contains(tract,'UF') || contains(tract,'ATR')
%column_max_ant=find(vertices(3,:) >  55 & vertices(1,:) == 60  );
%if side == 'L'
 %impongo che i percentili inizino dalla parte frontale con valore 0 fino
 %alla parte posteriore con valore 1
column_max_ant=find(vertices(2,:) == max(vertices(2,:))  );
 if (round(mean(d(1, column_max_ant))) > (num_percentili/3) )
     d=-(d-num_percentili*ones(size(d,1),size(d,2)));
 %end
end
  end

  if contains(tract,'GCC')
%column_max_ant=find(vertices(3,:) >  55 & vertices(1,:) == 60  );
%if side == 'L'
 %impongo che i percentili inizino dalla parte frontale con valore 0 fino
 %alla parte posteriore con valore 1
column_max_side=find(vertices(1,:) == max(vertices(1,:))  );
 if (round(mean(d(1, column_max_side))) < (num_percentili/3) )
     d=-(d-num_percentili*ones(size(d,1),size(d,2)));
 %end
end
  end

  if contains(tract,'SCC')
%column_max_ant=find(vertices(3,:) >  55 & vertices(1,:) == 60  );
%if side == 'L'
 %impongo che i percentili inizino dalla parte frontale con valore 0 fino
 %alla parte posteriore con valore 1
column_max_side=find(vertices(1,:) == max(vertices(1,:))  );
 if (round(mean(d(1, column_max_side))) > (num_percentili/3) )
     d=-(d-num_percentili*ones(size(d,1),size(d,2)));
 %end
end
 end


 d=d+1;
 percentile_mesh=d' ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%from mesh parameterization to volume segmentation


NN=knnsearch(vertices',coord); %evaluate the closest

%associo l'autovalore
percentile=zeros(1,size(coord,1));
for i=1:size(coord,1)
   percentile(1,i)=d(1,NN(i,1));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(pwd)

figure,
%plotmesh(vertices,faces,'facecolor', 'w')
title(['mesh surface for subject ' subject 'tract ' tract ' ' 'brain side ' side] );
p=patch(s);
isonormals(vol,p); %good 3D rendering
daspect([1,1,1])
view(3); axis tight
camlight right; camlight left; lighting gouraud


%if side == 'R'
%axis([10 45 20 80 20 60]) %R
%else axis([50 85 20 80 20 60]) %L
%end
x_l=xlabel('y (mm)');
y_l=ylabel('x (mm)');
z_l=zlabel('z (mm)');

%if side == 'L'
%set(gca, 'YDir', 'reverse')
%end
if side == 'L'
view(145, -42);
end
if side == 'R'
view(-16, -49);
end
%set(gca, 'YDir', 'reverse')


saveas(gcf,[dir_data '/mesh' subject '_' tract '_' side '.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%volume parameterization image
figure,
scatter3(coord(:,1),coord(:,2),coord(:,3),40,percentile,'filled');
title(['scatter plot percentili volume for subject ' subject ' ' 'brain side ' side] );
%legend('percentile all AF volume points','percentile mesh','Location','northwest')
%if side == 'R'
%axis([10 45 20 80 20 60]) %R
%else axis([50 85 20 80 20 60]) %L
%end
x_l=xlabel('x (mm)');
y_l=ylabel('y (mm)');
z_l=zlabel('z (mm)');

%if side == 'L'
%set(gca, 'YDir', 'reverse')
%end
%set(gca, 'YDir', 'reverse')

daspect([1,1,1])
view(3); axis tight
camlight right; camlight left; lighting gouraud
if side == 'R'
view(125, 46);
end
colormap(jet);
colorbar;
saveas(gcf,[dir_data '/scatter_' tract 'vol' subject side '.fig'])


vol_color = zeros(size(tract_cut,1),size(tract_cut,2),size(tract_cut,3)) ;
%coord=cat(2,coord(:,2),coord(:,1),coord(:,3));


for i=1:size(coord,1)

    %coord(i,1)=size(vol,1)+1-coord(i,1);  %dovrei rivertire la destra e la sinistra perch� poi nella conversione a nifti me le reinverte
    vol_color(coord(i,1),coord(i,2),coord(i,3))= percentile(1,i);

end

 cd(pwd)
niftiwrite(single(vol_color), [dir_data '/Laplacian_seg_' tract '_' side ], info_nifti);


%%
%%
   %stats whole tract core RD AD

   if contains(tract,'GCC') || contains(tract,'SCC') || contains(tract,'Fornix')
       [FA_median, FA_iqr, FA_mask]=map_stat([tract '_' th '_FA_ontoMNI'], vol_color);
       [MD_median, MD_iqr, MD_mask]=map_stat([tract '_' th '_MD_ontoMNI'], vol_color);
       [AD_median, AD_iqr, AD_mask]=map_stat([tract '_' th '_AD_ontoMNI'], vol_color);
       [RD_median, RD_iqr, RD_mask]=map_stat([tract '_' th '_FA_ontoMNI'], vol_color);
   else
       [FA_median, FA_iqr, FA_mask]=map_stat([tract '_' side '_' th '_FA_ontoMNI'], vol_color);
       [MD_median, MD_iqr, MD_mask]=map_stat([tract '_' side '_' th '_MD_ontoMNI'], vol_color);
       [AD_median, AD_iqr, AD_mask]=map_stat([tract '_' side '_' th '_AD_ontoMNI'], vol_color);
       [RD_median, RD_iqr, RD_mask]=map_stat([tract '_' side '_' th '_FA_ontoMNI'], vol_color);
   end
volume=nnz(vol_color(:));
%salva risultato
Statistica_core_median=cat(2, volume.*0.008, MD_median.*1000, FA_median, AD_median.*1000, RD_median.*1000);
Statistica_core_iqr=cat(2, 0 , MD_iqr.*1000 , FA_iqr, AD_iqr.*1000, RD_iqr.*1000);
Statistica_core=cat(1, Statistica_core_median, Statistica_core_iqr);

    save([dir_data '/' tract '_' side '_' 'vol_MD_FA_AD_RD_core_lin.mat'],'Statistica_core');


%%
%stats along tract
'stats along tract'

seg_median_FA=zeros(num_percentili,1);
seg_median_MD=zeros(num_percentili,1);
seg_median_AD=zeros(num_percentili,1);
seg_median_RD=zeros(num_percentili,1);

seg_median_vol=zeros(num_percentili,1);
seg_cog=zeros(num_percentili,3);
for i=1:(num_percentili+1)

   mask=(vol_color==i);

    MD_seg=MD_mask .* mask;
    FA_seg=FA_mask .* mask;
    AD_seg=AD_mask .* mask;
    RD_seg=RD_mask .* mask;

    FA_seg=FA_seg(:); FA_seg=FA_seg(FA_seg~=0); seg_median_FA(i,1)= median(FA_seg);
    MD_seg=MD_seg(:); MD_seg=MD_seg(MD_seg~=0); seg_median_MD(i,1)= median(MD_seg);
    AD_seg=AD_seg(:); AD_seg=AD_seg(AD_seg~=0); seg_median_AD(i,1)= median(AD_seg);
    RD_seg=RD_seg(:); RD_seg=RD_seg(RD_seg~=0); seg_median_RD(i,1)= median(RD_seg);


    seg_median_FA(i,1)= median(FA_seg);
    seg_median_MD(i,1)= median(MD_seg);
    seg_median_AD(i,1)= median(AD_seg);
    seg_median_RD(i,1)= median(RD_seg);

    seg_median_vol(i,1)= nnz(mask(:));
    S=regionprops3(mask);
    index=find(max(S.Volume));
    seg_cog(i,:)=S.Centroid(index,:);
end

%save risult

along_cog_y=seg_cog(:,1);
along_cog_x=seg_cog(:,2);
along_cog_z=seg_cog(:,3);

along_cog_x=abs(along_cog_x-46).*2;
along_cog_y=(along_cog_y-(126/2+1)).*2;
along_cog_z=(along_cog_z-(72/2+1)).*2;

Statistica_along_tract=cat(2,seg_median_vol.*0.008, seg_median_MD.*1000, seg_median_FA,along_cog_x,along_cog_y,along_cog_z,seg_median_AD.*1000, seg_median_RD.*1000 );

    save([dir_data '/' tract '_' side '_' 'vol_MD_FA_cog_along_lin_AD_RD.mat'],'Statistic_along_tract');


    %{
%plot results along-tract
data=load(['data/' tract '_' side '_vol_MD_FA_cog_along_lin_AD_RD']);

pz_along_vol=data.Statistica_along_tract(:,1);
pz_along_MD=data.Statistica_along_tract(:,2);
pz_along_FA=data.Statistica_along_tract(:,3);
pz_along_cog_x=data.Statistica_along_tract(:,4);
pz_along_cog_y=data.Statistica_along_tract(:,5);
pz_along_cog_z=data.Statistica_along_tract(:,6);
pz_along_AD=data.Statistica_along_tract(:,7);
pz_along_RD=data.Statistica_along_tract(:,8);

'Figura FA'
 %%
 %inzio figure
%dati_ID=['risultati_' num_par];
n_seg=num_percentili+1;


figure('units','normalized','outerposition',[0 0 1 1])

if side == 'L'
color_intq=[0 0.8 1];
sig='-';
else
color_intq='r';
sig='--';
end
 x=1:n_seg;

h_paz=plot(x,pz_along_FA);
if side == 'L'
    set(h_paz,'Marker','.','color','b','MarkerSize',15);
else
set(h_paz,'Marker','.','color','r','MarkerSize',15);
end
 y_l=ylabel('FA');
 legend(ID);
%legend('controls','interquartile','patient','Location', 'EastOutside');
x_l=xlabel(x_label_text);
xlim([0 n_seg+1])
%ylim([0.2 0.6])
%ti=title(['FA ' side]);
%set(gca,'XGrid','on')
%set(gca,'FontSize',dim_carattere);set(x_l,'FontSize',dim_carattere);set(y_l,'FontSize',dim_carattere);
set(gca,'FontSize',dim_carattere); set(x_l,'FontSize',dim_carattere);set(y_l,'FontSize',dim_carattere);
title([tract ' ' side ' ' ID])


saveas(gcf,[dir_fig '/FA_' ID '.png'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot MD

'Plot MD'


%dati_ID=['risultati_' num_par];
figure('units','normalized','outerposition',[0 0 1 1])
 x=1:n_seg;
h_paz=plot(x,pz_along_MD);
if side == 'L'
    set(h_paz,'Marker','.','color','b','MarkerSize',15);
else
set(h_paz,'Marker','.','color','r','MarkerSize',15);
end
 y_l=ylabel('MD (mm^{2}/s) x 10^{-3}');
%legend('controls','interquartile','patient','Location', 'EastOutside');
x_l=xlabel(x_label_text);

xlim([0 n_seg+1])
%ylim([0.3 0.6])
%ti=title(['MD ' side]);
%set(gca,'XGrid','on')
set(gca,'FontSize',dim_carattere);set(x_l,'FontSize',dim_carattere);set(y_l,'FontSize',dim_carattere);
ax = gca;
%ax.YAxis.Exponent = -3;
set(gca,'FontSize',dim_carattere); set(x_l,'FontSize',dim_carattere);set(y_l,'FontSize',dim_carattere);
legend(ID);
title([tract ' ' side ' ' ID])

cd(dir_fig)
saveas(gcf,['MD_' ID '.png'])

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot vol

%dati_ID=['risultati_' num_par];
figure('units','normalized','outerposition',[0 0 1 1])
 x=1:n_seg;
h_paz=plot(x,pz_along_vol);
if side == 'L'
    set(h_paz,'Marker','.','color','b','MarkerSize',15);
else
set(h_paz,'Marker','.','color','r','MarkerSize',15);
end
 y_l=ylabel('volume (cm^{3})');
%legend('controls','interquartile','patient','Location', 'EastOutside');
x_l=xlabel(x_label_text);

xlim([0 n_seg+1])
%ylim([0 1200])


%ti=title(['voxel ' side]);
%set(gca,'XGrid','on')
set(gca,'FontSize',dim_carattere);set(x_l,'FontSize',dim_carattere);set(y_l,'FontSize',dim_carattere);
legend(ID);
title([tract ' ' side ' ' ID])
cd(dir_fig)
saveas(gcf,['voxel_' ID '.png'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cog x

%dati_ID=['risultati_' num_par];
figure('units','normalized','outerposition',[0 0 1 1])
 x=1:n_seg;
h_paz=plot(x,pz_along_cog_x);
if side == 'L'
    set(h_paz,'Marker','.','color','b','MarkerSize',15);
else
set(h_paz,'Marker','.','color','r','MarkerSize',15);
end
 y_l=ylabel('|x| mm');
%legend('controls','interquartile','patient','Location', 'EastOutside');
x_l=xlabel(x_label_text);

xlim([0 n_seg+1])
%ylim([25 45])
%ylim([0.3 0.6])


%ti=title(['x coord ' side]);

%set(gca,'XGrid','on')
set(gca,'FontSize',dim_carattere);set(x_l,'FontSize',dim_carattere);set(y_l,'FontSize',dim_carattere);
%ax.YAxis.Exponent = -3;
legend(ID);
title([tract ' ' side ' ' ID])
cd(dir_fig)
saveas(gcf,['x_' ID  '.png'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cog y

%dati_ID=['risultati_' num_par];
figure('units','normalized','outerposition',[0 0 1 1])
n_seg=size(along_cog_y,1);
 x=1:n_seg;
h_paz=plot(x,pz_along_cog_y);
if side == 'L'
    set(h_paz,'Marker','.','color','b','MarkerSize',15);
else
set(h_paz,'Marker','.','color','r','MarkerSize',15);
end
 y_l=ylabel('y mm');
%legend('controls','interquartile','patient','Location', 'EastOutside');
x_l=xlabel(x_label_text);

xlim([0 n_seg+1])
%ylim([0.3 0.6])

title([tract '' side ' ' ID])
%set(gca,'XGrid','on')
set(gca,'FontSize',dim_carattere);set(x_l,'FontSize',dim_carattere);set(y_l,'FontSize',dim_carattere);
legend(ID);
%ax.YAxis.Exponent = -3;
cd(dir_fig)
saveas(gcf,['y_' ID  '.png'])


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot cog z
%dati_ID=['risultati_' num_par];
figure('units','normalized','outerposition',[0 0 1 1])
n_seg=size(along_cog_z,1);
 x=1:n_seg;

h_paz=plot(x,pz_along_cog_z);
if side == 'L'
    set(h_paz,'Marker','.','color','b','MarkerSize',15);
else
set(h_paz,'Marker','.','color','r','MarkerSize',15);
end

%legend('controls','interquartile','patient','Location', 'EastOutside');

%ylim([0.3 0.6])
%ti=title(['z coord ' side]);
y_l=ylabel('z (mm)');
x_l=xlabel(x_label_text);
xlim([0 n_seg+1])
%ylim([5 40])
%set(gca,'XGrid','on')
set(gca,'FontSize',dim_carattere);
%set(ti,'FontSize',dim_carattere);
set(gca,'FontSize',dim_carattere); set(x_l,'FontSize',dim_carattere);set(y_l,'FontSize',dim_carattere);
title([tract ' ' side ' ' ID])
legend(ID);
cd(dir_fig)
saveas(gcf,['z_' ID  '.png'])
%}

close all
cd(home_dir_orig)
end

function [map_median, map_iqr, map_mask]=map_stat(map_name, vol_laplacian)

         map=niftiread(map_name);

   mask_=vol_laplacian;
   mask_= mask_>0;

  map_mask=map.* mask_;

map_tract=map_mask(:); map_val_tract=map_tract(map_tract~=0);
map_median=median(map_val_tract(:));
map_iqr=iqr(map_val_tract(:));
end
