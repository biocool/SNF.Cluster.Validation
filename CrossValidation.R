
#Variables that Yaqoing used for doing SNF
#"yaqoing.input.data." is based on final clinical scores and "yaqoing.input.data.intake" is based on intake scores
ROI_var <- 
  colnames(yaqoing.input.data.intake)[c(5:8)]
clinic_var <- 
  colnames(yaqoing.input.data.intake)[c(9:14)]


# normalization variables for each yaqoing.input.data.intakea type

ROI_yaqoing.input.data.intake <- standardNormalization(yaqoing.input.data.intake[,ROI_var])
clinc_yaqoing.input.data.intakea <- standardNormalization(yaqoing.input.data.intake[,clinic_var])

# pair-wise distance
DistROIs = (dist2(as.matrix(ROI_yaqoing.input.data.intake),as.matrix(ROI_yaqoing.input.data.intake)))^(1/2)
DistClinic = (dist2(as.matrix(clinc_yaqoing.input.data.intakea),as.matrix(clinc_yaqoing.input.data.intakea)))^(1/2)

# add row and column names
rownames(DistROIs) <- yaqoing.input.data.intake$subj
rownames(DistClinic) <- yaqoing.input.data.intake$subj

colnames(DistROIs) <- yaqoing.input.data.intake$subj
colnames(DistClinic) <- yaqoing.input.data.intake$subj

# similarity graphs
WROIs <- affinityMatrix(DistROIs, K = 20, sigma = 0.5)
WClinic <- affinityMatrix(DistClinic, K = 20, sigma = 0.5)
yaqoing.data.TwoMainLayers.list <- list(WROIs, WClinic)




layers.list <- yaqoing.data.TwoMainLayers.list 
#cluster.index.vctr <- yaqoing.input.data$index
cluster.index.vctr <- yaqoing.input.data.intake$Cluster



accuracy.vctr.repeated.k.fold.CV.intake <- 
  SNF.repeated.cross.validation(no.of.repeat = 10, 
                              fold.No = 5, 
                              layers.list = layers.list, 
                              cluster.index.vctr = cluster.index.vctr)

save(accuracy.vctr.repeated.k.fold.CV.intake,
     file = "./Result/accuracy.vctr.repeated.k.fold.CV.intake")

# > mean(accuracy.vctr.repeated.k.fold.CV)
# [1] 0.839416
# > sd(accuracy.vctr.repeated.k.fold.CV)
# [1] 0.06239052
# > save(accuracy.vctr.repeated.k.fold.CV,
#        +      file = "./Result/accuracy.vctr.repeated.k.fold.CV")

SNF.repeated.cross.validation <- 
  function(no.of.repeat, fold.No, layers.list, cluster.index.vctr)
{
  accuracy.vctr.repeated.k.fold.CV <- c()
  for (i in c(1:no.of.repeat)) 
  {
    print("================================================")
    print("iteration:")
    print(i)
    accuracy.vctr.k.fold.CV <- 
      SNF.cross.validation(fold.No = 5, 
        layers.list = layers.list, 
        cluster.index.vctr = cluster.index.vctr)
    accuracy.vctr.repeated.k.fold.CV <- 
      c(accuracy.vctr.repeated.k.fold.CV, accuracy.vctr.k.fold.CV)
    
  }
  return(accuracy.vctr.repeated.k.fold.CV)
  
}



SNF.cross.validation <- function(fold.No, layers.list , cluster.index.vctr)
{
  

  flds <- createFolds( c(1:length(yaqoing.data.cluster.index.vctr)), 
                       k = fold.No, list = TRUE, returnTrain = FALSE)
  accuracy.vctr <- c()
  i<-1
  for (fold.index in flds) 
  {
    print("Fold")
    print(i)
    testSampleIndexVctr <- 
      fold.index
    accuracy.for.this.fold <- 
      SNF.Clustering.Validation.Train.Test(layers.list = yaqoing.data.TwoMainLayers.list, 
                                           cluster.index.vctr.from.SNF.on.WholeData = yaqoing.data.cluster.index.vctr, 
                                           testSampleIndexVctr = testSampleIndexVctr)
    accuracy.vctr <- 
      c(accuracy.vctr, accuracy.for.this.fold)
    i <- i + 1
  }
  return(accuracy.vctr)
}

SNF.Clustering.Validation.Train.Test <- 
  function(layers.list, cluster.index.vctr.from.SNF.on.WholeData , testSampleIndexVctr)
{
  # Create train and test data
  no.of.train.samples <- length(cluster.index.vctr.from.SNF.on.WholeData) - length(testSampleIndexVctr)# number of training cases
  no.of.test.samples <- length(cluster.index.vctr.from.SNF.on.WholeData) - no.of.train.samples
  train = lapply(layers.list, function(x) x[-testSampleIndexVctr, ]) # Use a list of subset of the layers (views) for training
  test = lapply(layers.list, function(x) x[testSampleIndexVctr, ]) # Test the rest of the data set
  groups = cluster.index.vctr.from.SNF.on.WholeData[-testSampleIndexVctr] #labels of the training data
  # Set the other parameters
  K = 20
  alpha = 0.5
  t = 20
  # method	
  # A indicator of which method to use to predict the label. method = 0 means to use local and global consistency; method = 1 means to use label propagation.
  method = TRUE
  # Apply the prediction function to the data
  #length of the "newcluster.index.vctr.from.SNF.done.on.WholeData" is equal the number of training data
  #the first "no.of.train.samples" elements of this vector and the rest shows the predicted clusters for the test data
  newcluster.index.vctr.from.SNF.done.on.WholeData = groupPredict(train,test,groups,K,alpha,t,method)

  # Compute the prediction accuracy====
  #"cluster.index.vctr.from.SNF.on.WholeData[testSampleIndexVctr]" is the cluster labels that were assigned to the 
  no.of.correctly.classified.test.instances <- sum(cluster.index.vctr.from.SNF.on.WholeData[testSampleIndexVctr] == 
        newcluster.index.vctr.from.SNF.done.on.WholeData[-c(1:no.of.train.samples)])
  
  accuracy <- (no.of.correctly.classified.test.instances / no.of.test.samples )
  
  return(accuracy)
}


################################################################################################################################################################################################################################################################################################
#SNFtool's Example====
################################################################################################################################################################################################################################################################################################
# Provide an example of predicting the new labels with label propagation
# Load views into list "dataL" and the cluster assignment into vector "label"
data(dataL)
dim(dataL)
View(dataL)#it contains two views (layers)
data(label)
# Create the training and test data
n = floor(0.8*length(label)) # number of training cases
trainSample = sample.int(length(label), n)
train = lapply(dataL, function(x) x[trainSample, ]) # Use the first 150 samples for training
test = lapply(dataL, function(x) x[-trainSample, ]) # Test the rest of the data set
groups = label[trainSample]
# Set the other
K = 20
alpha = 0.5
t = 20
method = TRUE
# Apply the prediction function to the data
newLabel = groupPredict(train,test,groups,K,alpha,t,method)
# Compare the prediction accuracy
accuracy = sum(label[-trainSample] == newLabel[-c(1:n)])/(length(label) - n)
################################################################################################################################################################################################################################################################################################
################################################################################################################################################################################################################################################################################################
