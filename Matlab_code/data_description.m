%%%%%%%%%%%%%%%%%%%%% Description for the '.mat' datasets %%%%%%%%%%%%%%%%%%%%%

%% data example 1
dataA = load('data_example1.mat');
X = dataA.X;
Y = dataA.Y;
lab = dataA.lab;
group = dataA.group;
true_pattern = dataA.true_pattern;
% In this simulated example, the confounder is a factor with three levels, 
% corresponding to three groups. The confounder contributes globally to all 
% variables(genes). Here are details for the data:
% X: the N by p data matrix, number of samples=30, number of variables=400
% Y: the N by q confounder matrix, q=3, Y is chosen such that the penalty 
% term equals the between groups sum of squares.
% lab: the labels (used in the plots in the user's guide)
% group: the group labels
% true_pattern: the true underlying latent pattern
%%

%% data example 2
dataA = load('data_example2.mat');
X = dataA.X;
Y = dataA.Y;
lab = dataA.lab;
true_pattern = dataA.true_pattern;
confound_pattern = dataA.confound_pattern;
% In this simulated example, the confounder is assumed to be continuous and 
% it contributes a global trend to the data. Here are details for the data:
% X: the N by p data matrix, number of samples=10, number of variables=400
% Y the N by q confounder matrix, q=1, 
% representing a continuous confounder. Y is centered.
% lab: the labels (used in the plots in the user's guide)
% true_pattern: the true underlying latent pattern
% confound_pattern: the pattern for the confounder
%%

%% data example 3
dataA = load('data_example3.mat');
X = dataA.X;
Y = dataA.Y;
lab = dataA.lab;
true_pattern = dataA.true_pattern;
% The confounder is unobserved and we only know the primary variable of 
% interest (the biological conditions). There are 10 biological conditions, 
% each with 3 replicates. The variation is shared among replicates for half 
% of the genes and not shared for the other genes. 
% Here are details for the data:
% X: the N by p data matrix, number of samples=30, number of variables=400
% Y: the N by q confounder matrix, q=30. For a biological condition, 
% treating the 3 replicates as 3 groups, it can be shown that the penalty term 
% equals the summation of the between groups sum of squares over the biological conditions.
% Details are provided in the user's guide.
% lab: labels for the biological conditions.
% true_pattern: the true underlying latent pattern
%%

%% The human brain exon array dataset, time window 2
dataA = load('data_brain_w2.mat');
X = dataA.X;
Y = dataA.Y;
Yid = dataA.Yid;
regions = dataA.regions; 
hemispheres = dataA.hemispheres;
donor_labs = dataA.donor_labs;
% A subset of 1,000 genes are included for demonstration purpose.
% Here are details for the data:
% X: the N by p data matrix, N is the number of samples and p is the number of genes.
% Y: the N by q confounder matrix (See implementation details in the user's guide).
% Yid: labels for the individuals, left and right hemispheres from the same 
% donors are treated as different individuals (See implementation details 
% in the user's guide).
% regions: labels for the brain regions
% hemispheres: '1' represents left hemisphere, '3' represents right hemisphere
% donor_labs: labels for the donors
%%

%% The human brain exon array dataset, time window 5
dataA = load('data_brain_w5.mat');
X = dataA.X;
Y = dataA.Y;
regions = dataA.regions; 
hemispheres = dataA.hemispheres;
donor_labs = dataA.donor_labs;
% A subset of 1,000 genes are included for demonstration purpose.
% Here are details for the data:
% X: the N by p data matrix, N is the number of samples and p is the number of genes.
% Y: the N by q confounder matrix (See implementation details in the user's guide).
% regions: labels for the brain regions
% hemispheres: '1' represents left hemisphere, '3' represents right hemisphere
% donor_labs: labels for the donors
%%

%% The modENCODE RNA-Seq data, fly and worm, embryonic stage
dataA = load('data_fly_worm.mat');
data_fly = dataA.data_fly;
data_worm = dataA.data_worm;
fly_time = dataA.fly_time;
worm_time = dataA.worm_time;
X = dataA.X;
X_time = dataA.X_time;
X_species = dataA.X_species;
Y = dataA.Y;
fly_gene = dataA.fly_gene;
% The modENCODE RNA-Seq data, fly and worm, embryonic stage. For the 
% orthologs in fly that map to multiple orthologs in worm, we took median 
% to get a one to one match, resulting in 4831 ortholog paris. Data downloaded 
% from https://www.encodeproject.org/comparative/transcriptome/.
% data_fly: expression levels for fly.
% data_worm: expression levels for worm.
% fly_time: time windows for fly. In the unit of hours.
% worm_time: time windows for worm In the unit of hours.
% X: data matrix for fly and worm combined. To gain robustness, 
% we used the rank across samples within the same species. 
% The rank matrix was then scaled to have unit variance.
% X_time: labels for the time windows/points
% X_species: labels for the species
% Y: confounder matrix (see user's guide for details).
% fly_gene: the fly ortholog
%%