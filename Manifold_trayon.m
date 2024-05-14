
clear
close all
% %% Example data
% DATA_PATH='E:\GPS\cal\BP_FIT\island-master\island-master\SampleData\HeadDirectionData'
% load([DATA_PATH '\full_neuron_firing_per_bin']) %num of spikes cells-by-time
% load([DATA_PATH '\angle_per_temporal_bin'])
% % load([DATA_PATH '\position_per_temporal_bin'])
% load([DATA_PATH '\spike_rate_mat_neuron_by_angle'])
% constants


% 00008_14072201_T4C1_egoN
% 00011_31122101_T6C1_ego

%% BP data
% DataName='BPData_2_106_Total_.mat';
DataName_Event={'Nosepoke_Rcount_1','Nosepoke_Rcount_2','Nosepoke_Rcount_3','Trayon_pump','Trayon_nonpump'};
DataName='HighDimData_Nosepoke_Rcount_1.mat';

load('HighDimData_Trayon_nonpump.mat')
HighDimData1=HighDimData;
load('HighDimData_Trayon_pump.mat')
% HighDimData1=HighDimData;
% load('HighDimData_Nosepoke_Rcount_3.mat')

HighDimData=[HighDimData1;HighDimData];


cmap21=colormap(othercolor('PuBu6',200));
cmap21 = [cmap21; 0 0 0];
cmap22=colormap(othercolor('PiYG3',200));
cmap22 = [cmap22; 0 0 0];
cmap23=colormap(othercolor('YlOrBr7',200));
cmap23 = [cmap23; 0 0 0];
color_index = round(1 + ((Timelag_per_bin - min(Timelag_per_bin)) / (max(Timelag_per_bin) - min(Timelag_per_bin))) * (200 - 1));
l=length(HighDimData)/2;

Timelag_per_bin=[Timelag_per_bin;Timelag_per_bin];

[pathstr, session, ext] = fileparts(DataName);
% load(DataName)
constantsBP
% Timelag_per_bin(201:end)=Timelag_per_bin(201:end)+7;
full_neuron_firing_per_bin=HighDimData;
%%
number_of_neurons = size(full_neuron_firing_per_bin, 2);

% Filter the neurons firing
% filtered_full_neuron_firing = filter_neuron_firing(full_neuron_firing_per_bin);
filtered_full_neuron_firing = full_neuron_firing_per_bin; % SR: No filter

%% Filter according to frames which contain a minimum number of active neurons
% (filter out quiet frames).

% Count number of active neurons per sample
number_of_active_neurons = sum(filtered_full_neuron_firing, 2);

% SR: Set NUMBER_OF_ACTIVE_NEURONS_THRESHOLD = 0. No mask
filter_mask = number_of_active_neurons > NUMBER_OF_ACTIVE_NEURONS_THRESHOLD;



%% Truncate data according to actual frames being used

filtered_full_neuron_firing_per_bin = full_neuron_firing_per_bin(filter_mask, :);
filtered_neuron_firing = filtered_full_neuron_firing(filter_mask, :);

% filtered_angle_per_temporal_bin = angle_per_temporal_bin(filter_mask)';
% filtered_position_per_temporal_bin = position_per_temporal_bin(filter_mask)';

%% Reduce data
% Original P values
P_NEIGHBORS_VEC = [0.075 / 20 0.075];
NUMBER_OF_REDUCED_DIMENSIONS_VEC = [10 10];
    
full_reduced_data = create_reduced_data(filtered_neuron_firing, P_NEIGHBORS_VEC, NUMBER_OF_REDUCED_DIMENSIONS_VEC);
% Take the final results and continue processing on them
reduced_data = full_reduced_data{length(P_NEIGHBORS_VEC) + 1};

% Use decoder to estimate head direction
% estimated_head_direction_angle_per_sample_index = estimate_head_direction(spike_rate_mat_neuron_by_angle, full_neuron_firing_per_bin);

% Filter the estimated head direction
% estimated_head_direction_angle_per_sample_index = estimated_head_direction_angle_per_sample_index(filter_mask);

% Handle missing behavioral entries
if INCLUDE_UNIDENTIFIED_ANGLES == false
%     reduced_data = reduced_data(~isnan(filtered_angle_per_temporal_bin), :);
%     filtered_angle_per_temporal_bin = filtered_angle_per_temporal_bin(~isnan(filtered_angle_per_temporal_bin));
%     estimated_head_direction_angle_per_sample_index = estimated_head_direction_angle_per_sample_index(~isnan(filtered_angle_per_temporal_bin));
end

% head_direction_neurons_indices = find_head_direction_neurons(spike_rate_mat_neuron_by_angle);

%%
% Plot the unlabeled reduced data
f10=figure;
plot3(reduced_data(:,2),reduced_data(:,3),reduced_data(:,4),'.'); %% 忽略掉第一个向量

% Plot the angle on the reduced data
cmap2 = hsv(NUMBER_OF_ANGLE_BINS);
% Add black color for 'nan' values
cmap2 = [cmap2; 0 0 0];

% For sleeping behavioral mode (rem and sws) we use the decoded head
% direction rather than actual head direction (which would probably be
% constant).
% visualization_angle_per_temporal_bin = filtered_angle_per_temporal_bin;

% visualization_angle_per_temporal_bin(find(visualization_angle_per_temporal_bin<330/180*pi&visualization_angle_per_temporal_bin>30/180*pi))=max(visualization_angle_per_temporal_bin)/2;
% visualization_angle_per_temporal_bin(find(visualization_angle_per_temporal_bin<60/180*pi|visualization_angle_per_temporal_bin>120/180*pi))=0;
% visualization_angle_per_temporal_bin(find(visualization_angle_per_temporal_bin<150/180*pi|visualization_angle_per_temporal_bin>210/180*pi))=0;
% visualization_angle_per_temporal_bin(find(visualization_angle_per_temporal_bin<240/180*pi|visualization_angle_per_temporal_bin>300/180*pi))=0;
% index_of_visualization_angle_per_temporal_bin = round(NUMBER_OF_ANGLE_BINS * visualization_angle_per_temporal_bin / ( 2 * pi));
% index_of_visualization_angle_per_temporal_bin(index_of_visualization_angle_per_temporal_bin == 0) = NUMBER_OF_ANGLE_BINS;
% % Color the missing values in black
% index_of_visualization_angle_per_temporal_bin(isnan(index_of_visualization_angle_per_temporal_bin)) = NUMBER_OF_ANGLE_BINS + 1; % 0-120
% The fourth argument is the dot size
f1=figure;
% scatter3(reduced_data(:, 2), reduced_data(:, 3), reduced_data(:, 4), 20, Timelag_per_bin, 'fill');


scatter3(reduced_data(1:l, 2), reduced_data(1:l, 3), reduced_data(1:l, 4), 20, cmap21(color_index, :), 'fill');
hold on 
scatter3(reduced_data(l+1:2*l, 2), reduced_data(l+1:2*l, 3), reduced_data(l+1:2*l, 4), 20, cmap22(color_index, :), 'fill');
% scatter3(reduced_data(2*l+1:3*l, 2), reduced_data(2*l+1:3*l, 3), reduced_data(2*l+1:3*l, 4), 20, cmap23(color_index, :), 'fill');

% colorm0ap(cmap2)

%% K-means
Threshold=10;
k=100;
point_cloud = [reduced_data(:, 2), reduced_data(:, 3), reduced_data(:, 4)];
opts = statset('Display','final');
[idx,C] = kmeans(point_cloud,k,'Distance','cityblock',...
    'Replicates',5,'Options',opts);
HistctCend=histcounts(idx,1:200);
ThresholdIndex=HistctCend>Threshold;
figure
scatter3(C(ThresholdIndex,1),C(ThresholdIndex,2),C(ThresholdIndex,3))
CThre=C(ThresholdIndex,:);
% save('E:\GPS\cal\BP_FIT\island-master\matlab_examples\tutorial_examples\Reduced_data.mat','CThre')

f2=figure;
scatter(reduced_data(:, 2), reduced_data(:, 3), 5, Timelag_per_bin, 'fill');
axis equal;


x=filtered_full_neuron_firing_per_bin;
[V,E,r,R]=PCA(x);
f3=figure







% colormap(othercolor('PuBu6'));
scatter3(x(1:l,:)*V(:,1),x(1:l,:)*V(:,2),x(1:l,:)*V(:,3), 20, cmap21(color_index, :), 'fill');


hold on
% colormap(othercolor('Set27'));
scatter3(x(l+1:2*l,:)*V(:,1),x(l+1:2*l,:)*V(:,2),x(l+1:2*l,:)*V(:,3), 20, cmap22(color_index, :), 'fill');

% colormap(othercolor('YlOrBr7'));
% scatter3(x(2*l+1:3*l,:)*V(:,1),x(2*l+1:3*l,:)*V(:,2),x(2*l+1:3*l,:)*V(:,3), 20, cmap23(color_index, :), 'fill');

% scatter3(x*V(:,1),x*V(:,2),x*V(:,3), 20, Timelag_per_bin, 'fill');
% scatter3(x*V(:,1),x*V(:,2),x*V(:,3), 20, Timelag_per_bin, 'fill');

f4=figure
scatter(x*V(:,1),x*V(:,2), 5, Timelag_per_bin, 'fill');


f5=figure
[pc, score, latent, tsquare]=pca(x);
scatter3(score(:,1),score(:,2),score(:,3), 20, Timelag_per_bin, 'fill');
f6=figure
scatter(score(:,1),score(:,2), 5, Timelag_per_bin, 'fill');


%% K-Means clustering
clustering_labels = k_means_clustering(reduced_data(:, CLUSTERING_DIMENSIONS), NUMBER_OF_CLUSTERS, 1);

%% Clustering visualization
cmap_clustering = jet(NUMBER_OF_CLUSTERS);
f7=figure;
scatter(reduced_data(:, 2), reduced_data(:, 3), 5, cmap_clustering(clustering_labels, :), 'fill');

axis equal;
box;

xlabel('Comp. 1');
ylabel('Comp. 2');


%% Create transition matrix and order it
transition_index_vec = clustering_labels(1:end - 1) + (clustering_labels(2:end) - 1) * NUMBER_OF_CLUSTERS;
[transition_index_count, ~] = histcounts(transition_index_vec, [0.5:1:NUMBER_OF_CLUSTERS^2 + 0.5]);
transition_index_mat = reshape(transition_index_count, [NUMBER_OF_CLUSTERS NUMBER_OF_CLUSTERS])';

transition_mat = transition_index_mat ./ repmat(sum(transition_index_mat, 2), [1 NUMBER_OF_CLUSTERS]);

% Run ordering code to get correct shuffling of the matrix
Ordering_cyclical_Ver000

ordered_clustering_labels = zeros(size(clustering_labels));

for cluster_index = 1:NUMBER_OF_CLUSTERS
    cluster_indices = find(clustering_labels == chosen_shuffle(cluster_index));
    
    if MIRROR_ORDERING == true
        ordered_clustering_results(cluster_indices) = NUMBER_OF_CLUSTERS + 1 - cluster_index;
    else
        ordered_clustering_results(cluster_indices) = cluster_index;
    end
end

estimated_angle_by_clustering = mod(CENTER_OF_CLUSTERING_ANGLE_BINS(ordered_clustering_results) + ACTUAL_VERSUS_CLUSTERING_SHIFT, 2 * pi);

%% Truncate data according to actual frames being used
full_estimated_angle_by_clustering = zeros(size(full_neuron_firing_per_bin, 1), 1);
full_estimated_angle_by_clustering(filter_mask) = estimated_angle_by_clustering;

smoothed_estimated_angle_by_clustering = smooth_estimated_angle(full_estimated_angle_by_clustering, filter_mask)';

%% Create tuning curves by clusters
neuron_by_cluster_spike_count = zeros(NUMBER_OF_CLUSTERS, number_of_neurons);

% Count the number of frames of each cluster
frames_per_cluster_count = histcounts(clustering_labels,  0.5:1:NUMBER_OF_CLUSTERS + 0.5);

for cluster_index = 1:NUMBER_OF_CLUSTERS
    cluster_frames_indices = find(clustering_labels == cluster_index);
    
    neuron_by_cluster_spike_count(cluster_index, :) = sum(filtered_neuron_firing(cluster_frames_indices, :), 1);
end

neuron_firing_rate = (neuron_by_cluster_spike_count ./ repmat(frames_per_cluster_count', [1 number_of_neurons])) * (BEHAVIORAL_SAMPLE_RATE / BEHAVIORAL_SAMPLES_PER_TEMPORAL_BIN);

ordered_neuron_firing_rate = neuron_firing_rate(chosen_shuffle, :);

% %% Plot the estimated polar plots
% firing_rate = create_firing_rate_matrix(filtered_full_neuron_firing_per_bin, smoothed_estimated_angle_by_clustering);
% plot_polar_tuning_curve(firing_rate, head_direction_neurons_indices);
% plot_polar_tuning_curve(spike_rate_mat_neuron_by_angle, head_direction_neurons_indices);
% 
% % Calculate correlation over actual and inferred tuning curve
% number_of_angle_bins = size(spike_rate_mat_neuron_by_angle, 2);
% 
% correlations = ones(1, number_of_neurons);
% 
% % TODO: This should be separated in case the number of bins of the actual
% % and the estimated differs
% CENTER_OF_ANGLE_BINS = [0.5 * (2 * pi) / number_of_angle_bins:...
%                         (2 * pi) / number_of_angle_bins:...
%                         2 * pi - 0.5 * (2 * pi) / number_of_angle_bins];
%                     
% CENTER_OF_HISTOGRAM_BINS = -0.9:0.2:0.9;
% 
% for neuron_index = 1:number_of_neurons
%     current_neuron_actual_firing_rate = spike_rate_mat_neuron_by_angle(neuron_index, :);
%     current_neuron_estimated_firing_rate = firing_rate(neuron_index, :);
%     
%     correlations(neuron_index) = corr(current_neuron_actual_firing_rate', current_neuron_estimated_firing_rate');
% end
% 
% f8=figure;
% hist(correlations, 8);
% 
% hist_mat = [hist(correlations(head_direction_neurons_indices), CENTER_OF_HISTOGRAM_BINS); ...
%             hist(correlations(~ismember(1:number_of_neurons, head_direction_neurons_indices)), CENTER_OF_HISTOGRAM_BINS)];
% 
% f9=figure;
% 
% colormap jet;
% bar(CENTER_OF_HISTOGRAM_BINS, hist_mat', 'stacked');
% 
% xlabel('Correlation');
% ylabel('Count');

exportgraphics(f1,[session '_1' '.jpg'],'Resolution',300)
exportgraphics(f2,[session '_2' '.jpg'],'Resolution',300)
% exportgraphics(f3,[session '_3' '.jpg'],'Resolution',300)
% exportgraphics(f4,[session '_4' '.jpg'],'Resolution',300)
exportgraphics(f5,[session '_5' '.jpg'],'Resolution',300)
exportgraphics(f6,[session '_6' '.jpg'],'Resolution',300)
exportgraphics(f7,[session '_7' '.jpg'],'Resolution',300)
% exportgraphics(f8,[session '_8' '.jpg'],'Resolution',300)
% exportgraphics(f9,[session '_9' '.jpg'],'Resolution',300)
