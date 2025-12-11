clc
close all
clear all

% wsi size
wsi_size = [53923 91737];

% MPP value of WSI
pix_dims = 0.2125;

% load in bladder tissue image from Xenium
tissue_img = imread("Full_Tissue_Export.png");

% derive tissue whole-slide mask and define non-tissue region
tissue_mask = bwpropfilt(bwfill(imgaussfilt(rgb2gray(tissue_img),2)>3,'holes'),'Area',1);
non_tissue = ~tissue_mask;

% load in bladder lumen mask
lumen_mask = imread("Bladder_Lumen_Mask.png");

% check mask dims to be whole-slide scale
lumen_mask = imresize(lumen_mask,wsi_size);
tissue_mask = imresize(tissue_mask,wsi_size);

% compute distance transform and apply jet color
dist_transform = bwdist(lumen_mask);
dist_transform(tissue_mask==0)=0;
figure,imagesc(dist_transform),colorbar();
hold on;
colormap(flipud(jet));

% compute lumen boundary points, distance, and convert to microns
[lumen_boundary,lumen_label] = bwboundaries(lumen_mask);
lumen_boundary_unstack = [lumen_boundary{1};lumen_boundary{2};lumen_boundary{3}];
lumen_boundary_microns = pix_dims.*lumen_boundary_unstack;

% load in centroid data
centroid_data = readtable("path_to_table.csv");
trans_centroids = centroid_data{:,2:3};

% compute cell type and/or transcript distances relative to bladder lumen
lumen_dists = [];
for i = 1:length(trans_centroids);

    test_pt = trans_centroids(i,:);
    test_pt_x_shift = size(dist_transform,1)-(test_pt(1)./pix_dims);
    test_pt_y_shift = (test_pt(2)./pix_dims);
    test_pt = [test_pt_x_shift,test_pt_y_shift];
    
    dist_calc = pdist2(test_pt,lumen_boundary_unstack);
    dist_min = min(dist_calc);
    dist_min_microns = dist_min*pix_dims;
    
    lumen_dists = [lumen_dists;dist_min_microns];
    
    if mod(i,1000) == 0;
        disp(i);
    end
end 

centroid_data.Lumen_Dists_Microns = lumen_dists;

% save the results appending lumen distances
writetable(centroid_data,"out_path_for_table.csv","WriteVariableNames",true);

%% END