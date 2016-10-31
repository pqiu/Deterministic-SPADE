%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this part of the code generates one scatter plot for the dendritic cells
% much of this part is copied from the code that generates Figure 2(b) in
% the main text

load('mouseBM.mat')
close all
h10 = figure(10);
pp = get(h10,'position'); set(h10,'position',[0 pp(2) pp(3)*1.7   pp(4)])

h = subplot(2,4,1);
FlowJo_contour2D(data(ismember(marker_names,'SSC-A'),:),data(ismember(marker_names,'CD11b'),:),5); 
title('All cells', 'FontSize', 14); xlabel('SSC-A', 'FontSize', 14); ylabel('CD11b', 'FontSize', 14); axis([-2 13 -3 9])
set(h,'xticklabel',[],'yticklabel',[])
Not_Myeloids = (data(ismember(marker_names,'CD11b'),:)>=-2.4 & data(ismember(marker_names,'CD11b'),:)<=2.49 & data(ismember(marker_names,'SSC-A'),:)>=1.8 & data(ismember(marker_names,'SSC-A'),:)<=6);
Myeloids = (data(ismember(marker_names,'CD11b'),:)>=2.63 & data(ismember(marker_names,'CD11b'),:)<=7 & data(ismember(marker_names,'SSC-A'),:)>=2 & data(ismember(marker_names,'SSC-A'),:)<=6);
axis([-1.5,13,-3,9])
text(0.7,7.9,'Myeloid cells','FontSize',14,'color',[110 110 110]/255,'Fontweight','bold')
line([ 1.8     6      6    1.8    1.8], ...
     [-2.4  -2.4   2.49   2.49   -2.4], ...
     'color','k','linewidth',2)
line([   2      6    6   2      2], ...
     [2.63   2.63    7   7   2.63], ...
     'color',[110 110 110]/255,'linewidth',2)
 
figure(10);h = subplot(2,4,2);
FlowJo_contour2D(data(ismember(marker_names,'B220'),Not_Myeloids),data(ismember(marker_names,'TCR-B'),Not_Myeloids),5); 
% title('Gated Lymphoid cells', 'FontSize',11);
xlabel('B220', 'FontSize',14); ylabel('TCR-\beta', 'FontSize',14); axis([-2.5 8.5 -1 6.5])
set(h,'xticklabel',[],'yticklabel',[])
Tcells = (Not_Myeloids & data(ismember(marker_names,'TCR-B'),:)>=2 & data(ismember(marker_names,'TCR-B'),:)<=5.5 & data(ismember(marker_names,'B220'),:)>=-2 & data(ismember(marker_names,'B220'),:)<=2.5);
Bcells = (Not_Myeloids & data(ismember(marker_names,'TCR-B'),:)>=-0.6 & data(ismember(marker_names,'TCR-B'),:)<=2 & data(ismember(marker_names,'B220'),:)>=3 & data(ismember(marker_names,'B220'),:)<=6);
Not_T_or_B = (Not_Myeloids & data(ismember(marker_names,'TCR-B'),:)>=-0.6 & data(ismember(marker_names,'TCR-B'),:)<=1.7 & data(ismember(marker_names,'B220'),:)>=-2 & data(ismember(marker_names,'B220'),:)<=1.8);
text(-1.2,6,'T cells','FontSize',14,'color',[15 144 69]/255,'fontweight','bold')
text(6.15,1,{'B', 'cells'},'FontSize',14,'color',[244 124 32]/255,'fontweight','bold')
line([  -2   2.5    2.5     -2     -2], ...
     [ 2.0   2.0    5.5    5.5    2.0], ...
     'color',[15 144 69]/255,'linewidth',2)
line([   3     6      6      3      3], ...
     [-0.6  -0.6      2      2   -0.6], ...
     'color',[244 124 32]/255,'linewidth',2)
line([  -2   1.8    1.8     -2     -2], ...
     [-0.6  -0.6    1.7    1.7   -0.6], ...
     'color','k','linewidth',2)

 
figure(10);h = subplot(2,4,3);
FlowJo_contour2D(data(ismember(marker_names,'CD4'),Tcells),data(ismember(marker_names,'CD8'),Tcells),5); 
% title('T cells', 'FontSize',11); 
xlabel('CD4', 'FontSize',14); ylabel('CD8', 'FontSize',14); axis([-2.5 6.8 -2 6])
set(h,'xticklabel',[],'yticklabel',[])
CD4cells = (Tcells & data(ismember(marker_names,'CD8'),:)>=-1.5 & data(ismember(marker_names,'CD8'),:)<=1.5 & data(ismember(marker_names,'CD4'),:)>=3.6 & data(ismember(marker_names,'CD4'),:)<=6);
CD8cells = (Tcells & data(ismember(marker_names,'CD8'),:)>=2.25 & data(ismember(marker_names,'CD8'),:)<=5 & data(ismember(marker_names,'CD4'),:)>=-2 & data(ismember(marker_names,'CD4'),:)<=2.8);
text(-2.2,5.6,'CD8^+ T cells','FontSize',14,'color',[0 174 239]/255,'fontweight','bold')
text(3.4,3.1,{'CD4^+'},'FontSize',14,'color',[122 39 120]/255,'fontweight','bold')
text(3.4,2.0,{'T cells'},'FontSize',14,'color',[122 39 120]/255,'fontweight','bold')
line([3.6    6      6    3.6     3.6], ...
     [-1.5  -1.5  1.5    1.5    -1.5], ...
     'color',[122 39 120]/255,'linewidth',2)
line([  -2   2.8    2.8     -2     -2], ...
     [2.25  2.25      5      5   2.25], ...
     'color',[0 174 239]/255,'linewidth',2)

figure(10);h = subplot(2,4,4);
FlowJo_contour2D(data(ismember(marker_names,'c-kit'),Not_T_or_B),data(ismember(marker_names,'Sca-1'),Not_T_or_B),2);
xlabel('c-kit', 'FontSize',14); ylabel('Sca-1', 'FontSize',14);axis([-1.5 7.5 -1.5 8])
set(h,'xticklabel',[],'yticklabel',[])
Stemcells = (Not_T_or_B & data(ismember(marker_names,'Sca-1'),:)>=2.8 & data(ismember(marker_names,'Sca-1'),:)<=6.5 & data(ismember(marker_names,'c-kit'),:)>=3.5 & data(ismember(marker_names,'c-kit'),:)<=7);
text(3.8,7.2,'HSPC','FontSize',14,'color',[231 55 37]/255,'fontweight','bold')
line([3.5     7       7    3.5     3.5], ...
     [2.8   2.8     6.5    6.5     2.8], ...
     'color',[231 55 37]/255,'linewidth',2)

annotation(figure(10),'arrow',[0.215 0.33],[0.65 0.65]);
annotation(figure(10),'arrow',[0.41 0.535],[0.825 0.825]);
annotation(figure(10),'line',[0.36,0.36],[0.597,0.51]);
annotation(figure(10),'line',[0.36,0.78],[0.51,0.51]);
annotation(figure(10),'arrow',[0.78 0.78],[0.51 0.57]);





load('SPADE_cluster_mst_upsample_result.mat','idx','mst_tree','node_positions','node_size','tree_bubble_contour','tree_annotations')
downsampled_assign = idx;

load('SPADE_cluster_mst_upsample_result.mat','all_assign')
all_assign = abs(all_assign{1});


subplot(2,4,[5:8])

[cluster_populate] = get_module_mean(double(Myeloids)', all_assign);  
cluster_populate = cluster_populate - (max(cluster_populate)+min(cluster_populate))/2;
subplot(5,6,[19 25]); 
draw_SPADE_tree_annotation(mst_tree, node_positions, node_size, cluster_populate, [-0.5 0.6], 1, 0, zeros(1,size(mst_tree,1)), 'gray', [],[])
title('Myeloid gate','Fontsize',14,'color',[110 110 110]/255,'fontweight','bold') 
draw_colorbar([25,-55,30,5],lbmapv2(100,'gray'),{'0%',' ','100%'})

[cluster_populate] = get_module_mean(double(Bcells)', all_assign); 
cluster_populate = cluster_populate - (max(cluster_populate)+min(cluster_populate))/2;
subplot(5,6,[20 26]); 
draw_SPADE_tree_annotation(mst_tree, node_positions, node_size, cluster_populate, [-0.5 0.6], 1, 0, zeros(1,size(mst_tree,1)), 'gray', [],[])
title('B cells gate','Fontsize',14,'color',[244 124 32]/255,'fontweight','bold') 
draw_colorbar([25,-55,30,5],lbmapv2(100,'gray'),{'0%',' ','100%'})

[cluster_populate] = get_module_mean(double(Tcells)', all_assign); 
cluster_populate = cluster_populate - (max(cluster_populate)+min(cluster_populate))/2;
subplot(5,6,[21 27]); 
draw_SPADE_tree_annotation(mst_tree, node_positions, node_size, cluster_populate, [-0.5 0.6], 1, 0, zeros(1,size(mst_tree,1)), 'gray', [],[])
title('T cells gate','Fontsize',14,'color',[15 144 69]/255,'fontweight','bold') 
draw_colorbar([25,-55,30,5],lbmapv2(100,'gray'),{'0%',' ','100%'})

[cluster_populate] = get_module_mean(double(CD8cells)', all_assign); 
cluster_populate = cluster_populate - (max(cluster_populate)+min(cluster_populate))/2;
subplot(5,6,[22 28]); 
draw_SPADE_tree_annotation(mst_tree, node_positions, node_size, cluster_populate, [-0.5 0.6], 1, 0, zeros(1,size(mst_tree,1)), 'gray', [],[])
title('CD8^+ gate','Fontsize',14,'color',[0 174 239]/255,'fontweight','bold') 
draw_colorbar([25,-55,30,5],lbmapv2(100,'gray'),{'0%',' ','100%'})

[cluster_populate] = get_module_mean(double(CD4cells)', all_assign); 
cluster_populate = cluster_populate - (max(cluster_populate)+min(cluster_populate))/2;
subplot(5,6,[23 29]); 
draw_SPADE_tree_annotation(mst_tree, node_positions, node_size, cluster_populate, [-0.5 0.6], 1, 0, zeros(1,size(mst_tree,1)), 'gray', [],[])
title('CD4^+ gate','Fontsize',14,'color',[122 39 120]/255,'fontweight','bold') 
draw_colorbar([25,-55,30,5],lbmapv2(100,'gray'),{'0%',' ','100%'})
 
[cluster_populate] = get_module_mean(double(Stemcells)', all_assign); 
cluster_populate = cluster_populate - (max(cluster_populate)+min(cluster_populate))/2;
subplot(5,6,[24 30]); 
draw_SPADE_tree_annotation(mst_tree, node_positions, node_size, cluster_populate, [-0.1 0.1], 1, 0, zeros(1,size(mst_tree,1)), 'gray', [],[])
title('HSPC gate','Fontsize',14,'color',[231 55 37]/255,'fontweight','bold') 
draw_colorbar([25,-55,30,5],lbmapv2(100,'gray'),{'0%',' ','10%'})






Bcells_branch = ismember(all_assign, tree_annotations{6});
DC_branch = ismember(all_assign, tree_annotations{7});
Tcell_branch = ismember(all_assign, unique([tree_annotations{2},tree_annotations{8},tree_annotations{9}]));
CD4cells_branch = ismember(all_assign, tree_annotations{8});
CD8cells_branch = ismember(all_assign, tree_annotations{2});
Myeloids_branch = ismember(all_assign, [tree_annotations{4},tree_annotations{5}]);
ckit_branch = ismember(all_assign, tree_annotations{1});


column_sum = [sum(Bcells), sum(Tcells), sum(CD4cells), sum(CD8cells), sum(Myeloids), sum(Stemcells)]

row_sum = [sum(Bcells_branch); sum(DC_branch); sum(Tcell_branch); sum(CD4cells_branch); sum(CD8cells_branch);  sum(Myeloids_branch);  sum(ckit_branch)]

counts = [sum(Bcells_branch & Bcells), sum(Bcells_branch & Tcells), sum(Bcells_branch & CD4cells), sum(Bcells_branch & CD8cells), sum(Bcells_branch & Myeloids), sum(Bcells_branch & Stemcells); ...
          sum(DC_branch & Bcells), sum(DC_branch & Tcells), sum(DC_branch & CD4cells), sum(DC_branch & CD8cells), sum(DC_branch & Myeloids), sum(DC_branch & Stemcells); ...
          sum(Tcell_branch & Bcells), sum(Tcell_branch & Tcells), sum(Tcell_branch & CD4cells), sum(Tcell_branch & CD8cells), sum(Tcell_branch & Myeloids), sum(Tcell_branch & Stemcells); ...
          sum(CD4cells_branch & Bcells), sum(CD4cells_branch & Tcells), sum(CD4cells_branch & CD4cells), sum(CD4cells_branch & CD8cells), sum(CD4cells_branch & Myeloids), sum(CD4cells_branch & Stemcells); ...
          sum(CD8cells_branch & Bcells), sum(CD8cells_branch & Tcells), sum(CD8cells_branch & CD4cells), sum(CD8cells_branch & CD8cells), sum(CD8cells_branch & Myeloids), sum(CD8cells_branch & Stemcells); ...
          sum(Myeloids_branch & Bcells), sum(Myeloids_branch & Tcells), sum(Myeloids_branch & CD4cells), sum(Myeloids_branch & CD8cells), sum(Myeloids_branch & Myeloids), sum(Myeloids_branch & Stemcells); ...
          sum(ckit_branch & Bcells), sum(ckit_branch & Tcells), sum(ckit_branch & CD4cells), sum(ckit_branch & CD8cells), sum(ckit_branch & Myeloids), sum(ckit_branch & Stemcells); ...
          ]
      
freq = counts./repmat(column_sum, 7, 1)

original_column_sum = column_sum;
original_row_sum = [152685; 3996; 17538; 2931; 6301; 202180; 15681];
original_counts = [146017   88      0       0       2246        16;
                   3562     79      77      0       83          0;
                   364      12377   2729    6037    3033        0;
                   0        2858    2713    0       27          0;
                   0        6174    0       5843    32          0;
                   5        75      0       0       199048      0;
                   8        815     2       9       1159        401;
                   ]
original_freq = original_counts./repmat(original_column_sum,7,1)



figure(20)
subplot(1,2,1);
plot(counts(:),original_counts(:),'o'); 
axis([1, 210000,1,210000]); line([1,210000],[1,210000],'color',[1 1 1]*0.7)
xlabel('Deterministic SPADE')
ylabel('Original SPADE')
title('Cell counts of overlap with manual gating')
subplot(1,2,2);
plot(freq(:),original_freq(:),'o'); axis([-0.5,1.5,-0.5,1.5]); line([-0.5,1.5],[-0.5,1.5],'color',[1 1 1]*0.7)
xlabel('Deterministic SPADE')
ylabel('Original SPADE')
title('Cell frequency/percentage of overlap with manual gating')

