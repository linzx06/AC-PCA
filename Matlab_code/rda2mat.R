library(acPCA)
library(R.matlab)

data("data_example1")
names(data_example1)
writeMat("data_example1.mat", X=data_example1$X, Y=data_example1$Y, lab=data_example1$lab, colors=data_example1$colors, true_pattern=data_example1$true_pattern)

data("data_example2")
names(data_example2)
writeMat("data_example2.mat", X=data_example2$X, Y=data_example2$Y, lab=data_example2$lab, colors=data_example2$colors, true_pattern=data_example2$true_pattern)

data("data_example3")
names(data_example3)
writeMat("data_example3.mat", X=data_example3$X, Y=data_example3$Y, lab=data_example3$lab, true_pattern=data_example3$true_pattern, confound_pattern=data_example3$confound_pattern)

data("data_example4")
names(data_example4)
writeMat("data_example4.mat", X=data_example4$X, Y=data_example4$Y, lab=data_example4$lab, true_pattern=data_example4$true_pattern)

data("data_example5")
names(data_example5)
writeMat("data_example5.mat", X=data_example5$X, Y=data_example5$Y)

data("data_brain_w1")
names(data_brain_w1)
writeMat("data_brain_w1.mat", X=data_brain_w1$X, Y=data_brain_w1$Y,
         Yid=data_brain_w1$Yid, regions=data_brain_w1$regions, hemispheres=data_brain_w1$hemispheres,
         donor_labs=data_brain_w1$donor_labs)

data("data_brain_w2")
names(data_brain_w2)
writeMat("data_brain_w2.mat", X=data_brain_w2$X, Y=data_brain_w2$Y,
         Yid=data_brain_w2$Yid, regions=data_brain_w2$regions, hemispheres=data_brain_w2$hemispheres,
         donor_labs=data_brain_w2$donor_labs)

data("data_fly_worm")
names(data_fly_worm)
writeMat("data_fly_worm.mat", 
         data_fly=data_fly_worm$data_fly, 
         data_worm=data_fly_worm$data_worm,
         fly_time=data_fly_worm$fly_time, 
         worm_time=data_fly_worm$worm_time, 
         X=data_fly_worm$X,
         Y=data_fly_worm$Y,
         X_time=data_fly_worm$X_time,
         X_species=data_fly_worm$X_species,
         fly_gene=data_fly_worm$fly_gene
         )