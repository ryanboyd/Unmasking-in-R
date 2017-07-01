


Unmask = function(data, criterion="Author", people_to_unmask = 'X', obs.id='Title', num_folds = 10,
                  features_to_drop = 3, seed=123, n.iter=5, verbose=T, graphs=T, unk.in.results=F,
                  unmask.unk.only=F){
  
  set.seed(seed)
  
  require(e1071)
  require(gtools)
  
  print('Starting the unmask function...')
  
  #this is a copy of the createFolds() function from the 'caret' package
  #I have put it in here to make sure that caret does not need to be installed
  #in order to run this script. great function, thanks to Max Kuhn and contributors!
  rand_group_assign = function (y, k = 10, list = TRUE, returnTrain = FALSE) 
  {
    if (class(y)[1] == "Surv") 
      y <- y[, "time"]
    if (is.numeric(y)) {
      cuts <- floor(length(y)/k)
      if (cuts < 2) 
        cuts <- 2
      if (cuts > 5) 
        cuts <- 5
      breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
      y <- cut(y, breaks, include.lowest = TRUE)
    }
    if (k < length(y)) {
      y <- factor(as.character(y))
      numInClass <- table(y)
      foldVector <- vector(mode = "integer", length(y))
      for (i in 1:length(numInClass)) {
        min_reps <- numInClass[i]%/%k
        if (min_reps > 0) {
          spares <- numInClass[i]%%k
          seqVector <- rep(1:k, min_reps)
          if (spares > 0) 
            seqVector <- c(seqVector, sample(1:k, spares))
          foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
        }
        else {
          foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                                                                 size = numInClass[i])
        }
      }
    }
    else foldVector <- seq(along = y)
    if (list) {
      out <- split(seq(along = y), foldVector)
      names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
                          sep = "")
      if (returnTrain) 
        out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
    }
    else out <- foldVector
    out
  }
  
  
  
  remove_zero_variance <- function(dat) {
    out <- lapply(dat, function(x) length(unique(x)))
    want <- which(!out > 1)
    unlist(want)
  }
  
  
  
  if (!people_to_unmask %in% data[, c(criterion)]) stop('The observations that you are trying to unmask\n do not exist in your author list.')
  if (num_folds < 2) stop('You must include at least 2 folds to run the analysis.')
  if (n.iter < 1) stop('You must include at least 1 iteration.')
  if (n.iter < 5) cat('\nRunning analysis. You probably want to increase n.iter to at least 5...\n\n')
  
  
  
  #this is where we do some basic sanity checking of the input to make sure that it's going to work
  if (features_to_drop == 'auto'){
    features_to_drop = floor(((ncol(data) - 2) / num_folds) / 2)
  } else {
    if (features_to_drop > floor(((ncol(data) - 2) / num_folds) / 2)){
      print('You are trying to drop too many features\nper iteration. Automatically resizing...')
      features_to_drop = floor(((ncol(data) - 2) / num_folds) / 2)
    }
  }
  
  
  #make sure that it's a factor
  data[, c(criterion)] = as.factor(data[, c(criterion)])
  data[, c(obs.id)] = as.factor(data[, c(obs.id)])
  
  cat('\nDropping all items that do not have enough observations\nto have at least one observation in each fold...\n\n')
  #first thing that we want to do is eliminate anything that doesn't have enough obvservations
  #for the number of folds that we're running
  obs.counts = as.data.frame(table(data[,c(obs.id)]))
  obs.to.keep = obs.counts[obs.counts$Freq >= num_folds,][[1]]
  data=subset(data, data[,c(obs.id)] %in% obs.to.keep)
  
  
  #this gets us two dataframes, one with the observations that we want to unmask,
  #and one with all of the remaining data
  #data = subset(data, data[, c(criterion)]!=people_to_unmask)
  
  data[, c(criterion)] = factor(data[, c(criterion)])
  

  #now, what we want to do is go through *each* author that we know,
  #and then go through *each* item that we want to unmask for each
  #author that we know, then do out cross-validation with removing items
  #this will get us a curve for each unknown item for *each* author that we know
  #who it actually is
  
  Known_Authors = levels(as.factor(data[, c(criterion)]))
  
  #do the same thing for the unknowns
  Items = levels(factor(data[,c(obs.id)]))
  
  
  #we'll go ahead and set up a dataframe that contains the results here
  Results_DF = data.frame(matrix(0, nrow=(length(Items)*length(Known_Authors)), ncol=num_folds+3))


  
  #here, we name the columns for the comparison authors and folds
    colnames(Results_DF)[1:3] = c(obs.id, criterion, 'same.different')
    for (i in 4:(num_folds+3)){
      colnames(Results_DF)[i] = paste('Fold', i-3, 'ACC', sep='.')
    }


  #create a temporary list to make sure that we're matching same/different author data
  temp_author_item_list = unique(data[, c(obs.id, criterion)])
  
  
  #and then we do the same thing for the rows
  
    for (Author in 1:(length(Known_Authors))){   
      
      for (Item in 1:(length(Items))){  
    
      
      Results_DF[Item + ((Author-1) * length(Items)), 1] = Items[Item]
      Results_DF[Item + ((Author-1) * length(Items)), 2] = Known_Authors[Author]
      
      #marks the observation as being by the same of a different author
      same_diff = (temp_author_item_list[, c(criterion)][[match(Items[Item], temp_author_item_list[, c(obs.id)])]] == Known_Authors[Author])
      
      if ((Known_Authors[Author] != people_to_unmask) & (temp_author_item_list[, c(criterion)][[match(Items[Item], temp_author_item_list[, c(obs.id)])]] != people_to_unmask)) {
        ifelse(same_diff==T, Results_DF[Item + ((Author-1) * length(Items)), 3] <- 'same', Results_DF[Item + ((Author-1) * length(Items)), 3] <- 'diff')
      } else {
        Results_DF[Item + ((Author-1) * length(Items)), 3] <- NA
      }
      
      }
    }
  


  for(iter in 1:n.iter){
  
    cat(paste('\n\nIteration ', iter, '...', '\n\n', sep=''))
    
    #Now we start looping through each author
    for (Author in 1:(length(Known_Authors))){
      
      Author_Name = Known_Authors[[Author]]
      
      #this subsets the full dataset down to the known author being tested against for this round
      Author_Data = subset(data, data[,c(criterion)]==Author_Name)
      Author_Data[ ,c(criterion)] = factor(Author_Data[ ,c(criterion)])
      
      
      #now that we've subset the data to a single author, we want to go through
      #and do each of the "unknown" items and compare against the author
      #remember that each unknown item can span several observations, so we're going to do
      #it by "title" (or whatever gets passed to "obs.id")
      
      for (Item in 1:(length(Items))){
        
        Item_Name = Items[Item]
        Unknown_Data = subset(data, data[,c(obs.id)] == Items[Item])
        
        #here, we figure out if we want to even compare this item to anything else
        Unknown_Data_Author_Name = Unknown_Data[1 ,c(criterion)]
        if ((unmask.unk.only==T) & (Unknown_Data_Author_Name != people_to_unmask)) next
        
        Author_Data_temp = subset(Author_Data, Author_Data[,c(obs.id)] != Items[Item])
        #assign a random string to the author for this, since we're treating it as distinct from any other author
        Unknown_Data[,c(criterion)] = paste(sample(c(0:9, letters, LETTERS), 20, replace=TRUE), collapse="")
        
        #skip this iteration if 
        if (nrow(Author_Data_temp) == 0){
          Results_DF[Item + ((Author-1) * length(Items)), 4:(num_folds+3)] = NA
          cat(paste(Author_Name, ' only has 1 item, so you cannot do any "same person" comparisons... Skipping...', sep=''))
          next
        }
        
        #we want to randomly remove rows so that there's an equal number of observations for both 
        #the known and unknown authors
        
        while(nrow(Author_Data_temp) > nrow(Unknown_Data)){
          Author_Data_temp = Author_Data_temp[-c(sample(1:nrow(Author_Data_temp), 1)), ]
        }
        
        while(nrow(Author_Data_temp) < nrow(Unknown_Data)){
          Unknown_Data = Unknown_Data[-c(sample(1:nrow(Unknown_Data), 1)), ]
        }
        
        
        #So, what we're left with here is this:
        #Author_Name = the name of the current author that we're looking at
        #Author_Data = data frame for the current author that we're looking at 
        #Item_Name = name of the unknown item
        #Unknown_Data = data for just the item that we're currently looking at
        #data = full known dataset
        #unmask_DF = full unknown dataset
        
        analysis_DF = smartbind(Author_Data_temp, Unknown_Data)
        
        #now we have analysis_DF -- the dataframe that we're going to do all of our magic with
        
        #set up each fold of our analysis
        analysis_DF$unmask_fold_group = rand_group_assign(analysis_DF[,c(criterion)], k=num_folds, list=F)
        
        if (nrow(analysis_DF) > nrow(Unknown_Data)){
        
          if (verbose==F){
            cat(paste(Author_Name, ' / ', Item_Name, '... \n', sep=''))
          }
        
          for (fold in 1:num_folds){
            
            if (verbose==T){
              cat(paste(Author_Name, ' / ', Item_Name, ': Fold ', fold, '... ', sep=''))
            }
            
            train_data = subset(analysis_DF, unmask_fold_group != fold)
            test_data = subset(analysis_DF, unmask_fold_group == fold)
            
            #remove the fold number variable from the data so that it's not included in the analysis
            #we also want to remove the observation Ids, since we don't want to use that for modeling
            train_data = train_data[ , -which(names(train_data) %in% c("unmask_fold_group", obs.id))]
            test_data = test_data[ , -which(names(test_data) %in% c("unmask_fold_group", obs.id))]
            
            train_data_y = train_data[ , which(names(train_data) %in% c(criterion))]
            train_data = train_data[ , -which(names(train_data) %in% c(criterion))]
            
            test_data_y = test_data[ , which(names(test_data) %in% c(criterion))]
            test_data = test_data[ , -which(names(test_data) %in% c(criterion))]
            
            #this leaves us with train_data and test_data, both of which should have a "criterion" column
            #(the name of which is passed to this function), and everything else should be numeric data
            #that came with it
            
            #so, what we need to do now is build a model on the training data, then test it on the
            #test data
            
            n <- names(train_data)
            f <- as.formula(paste("train_data_y ~", paste(n[!n %in% "y"], collapse = " + ")))
            
            model <- svm(f, data=train_data, scale=T, kernal='linear', degree=3, na.action=na.omit, type="C-classification")
            
            #once the model is built, we want to predict on the test data
            pred = predict(model, test_data)
            
            #calculate the accuracy of this model
            Acc_Evaluation = mean(as.numeric(test_data_y == pred)) * 100
            
            if (verbose==T){
              cat(paste('Accuracy = ', Acc_Evaluation, '\n', sep=''))
            }
            
            
            
            #append the accuracy to the Results_DF data
            Results_DF[Item + ((Author-1) * length(Items)), fold+3] = Results_DF[Item + ((Author-1) * length(Items)), fold+3] + Acc_Evaluation
            
            #get the parameter weights of each input variable
            param_weights = as.data.frame(t(model$coefs) %*% model$SV)
            
            
           
            #so, what we're left with now is param_weights, which tells us the parameter weights
            #of ONLY the variables that went into the analysis
            
            #let's sort them in descending order and start cutting variables
            vars_to_drop = character()
            param_weights = sort(param_weights ,decreasing=T)
            for (drop in 1:features_to_drop){
              vars_to_drop[length(vars_to_drop)+1] = names(param_weights[drop])
            }
            param_weights = sort(param_weights,decreasing=F)
            for (drop in 1:features_to_drop){
              vars_to_drop[length(vars_to_drop)+1] = names(param_weights[drop])
            }
            
            
            if (verbose==T){
              cat('\t\tDropping variables:\n\t')
              cat(vars_to_drop, sep=', ')
              cat('\n')
              
              #now we have a "vars_to_drop" vector that tells us the top N variables in either extreme. We can then drop these from
              #both the training and testing dataframes
              if (length(which(colnames(analysis_DF) %in% vars_to_drop)) > 0){
                analysis_DF = analysis_DF[ , -which(colnames(analysis_DF) %in% vars_to_drop)]
              }
              
              cat('\t\tVars remaining: ')
              cat(ncol(analysis_DF))
              cat('\n')
            } else {
              if (length(which(colnames(analysis_DF) %in% vars_to_drop)) > 0){
                analysis_DF = analysis_DF[ , -which(colnames(analysis_DF) %in% vars_to_drop)]
              }
            }
            
          }
        
        }
        
        
      }
      
      
      
    }
    
  }


  #average the accuracy across all iterations
  Results_DF[, c(4:length(Results_DF))] = Results_DF[, c(4:length(Results_DF))] / n.iter

  #prepare to return by cleaning up, adding some additional info into the frame itself
  Results_DF = Results_DF[ order(Results_DF[,2], Results_DF[,3]),]
  

  #crank out some additional data
   
 #accuracy difference between round i and i+1
 if (num_folds>1){
   temp=ncol(Results_DF)
   for(i in 1:(num_folds-1)){
     Results_DF[, c(paste('ACC_Delta', i, i+1, sep='-'))] = Results_DF[, c(i+3)] - Results_DF[, c(i+4)] 
   }
   
   Results_DF[, c('Max_ACC_Drop_1-iter')] = apply(as.data.frame(Results_DF[, c((temp+1):length(Results_DF))]), 1, max)
 }


 #accuracy difference between round i and i+2
  if (num_folds>2){
    temp=ncol(Results_DF)
    for(i in 1:(num_folds-2)){
      Results_DF[, c(paste('ACC_Delta', i, i+2, sep='-'))] = Results_DF[, c(i+3)] - Results_DF[, c(i+4)] 
    }
    Results_DF[, c('Max_ACC_Drop_2-iter')] = apply(as.data.frame(Results_DF[, c((temp+1):length(Results_DF))]), 1, max)
  }

  

  #here we run a quadratic model (x2 and x terms), then return the coefficients
  #to best describe the shape of the patterns
  lm2_func = function(incoming_data){
    
    tryCatch({
      t_dat = as.data.frame(incoming_data)
      colnames(t_dat) = ('y')
      t_dat$x = seq(1:num_folds)
      t_dat$x2 = t_dat$x**2
      linmod = lm(y~x2 + x, data=t_dat)
      return(linmod$coefficients[2:3])
    }, error = function(err) {return(c(NA, NA))})
  }

  #here we run a linear model (x terms), then return the coefficients
  #to best describe the shape of the patterns
  lm1_func = function(incoming_data){
    tryCatch({
      t_dat = as.data.frame(incoming_data)
      colnames(t_dat) = ('y')
      t_dat$x = seq(1:num_folds)
      linmod = lm(y~x, data=t_dat)
      return(linmod$coefficients[2])
    }, error = function(err) {return(NA)})
  }

  if (num_folds>1){
    lm_coef = apply(as.data.frame(Results_DF[ , c(4:(3+num_folds))]), 1, lm2_func)
    Results_DF$quad_x2_coef = lm_coef[1,]
    Results_DF$quad_x1_coef = lm_coef[2,]
    Results_DF$linear_x1_coef = apply(as.data.frame(Results_DF[ , c(4:(3+num_folds))]), 1, lm1_func)
  }
  
  

  if (graphs==T){
    
    tryCatch({
      for (Item in 1:(length(Items))){
      
        if (temp_author_item_list[, c(criterion)][[match(Items[Item], temp_author_item_list[, c(obs.id)])]] == people_to_unmask){
      
          graph_df = subset(Results_DF, (Results_DF[ , c(obs.id)] == Items[Item]) & (Results_DF[ , c(criterion)] != people_to_unmask))
          graph_df_data = graph_df[, c(4:(4+num_folds-1))]
          graph_df_data = t(graph_df_data)
      
          graph = matplot(graph_df_data, type = "l", xlab='Fold', ylab='Accuracy', main=Items[Item])
          legend("bottomleft", legend=graph_df[, c(criterion)], col=seq_len(ncol(graph_df)), cex=0.8, fill=seq_len(ncol(graph_df)))
        }
      }
    }, error = function(err){cat(paste('\nError plotting graph: ', Items[Item], '\n', sep=''))})
  }



  if (unk.in.results==T){
    return(Results_DF)
  } else {
    return(subset(Results_DF, Results_DF[, c(criterion)] != people_to_unmask))
  }


}







#read in the data
data=read.csv('X:/writing/Behn MEH/LIWC2015 Results (texts (96 files)).csv')



#You can sink the output to a file if you want to debug
#sink('x:/test.txt')



  #--data = dataset that you want to use. make sure that all undesired variables are already removed

  #--criterion = the name of the variable that you're trying to predict (e.g., "Author")
  
  #--people_to_unmask = the marker in the "criterion" column that identifies works of unknown authorship
  
  #--obs.id = the variable that is used to denote which text is which (e.g., book titles, etc.)

  #--num_folds = the number of time each text is tested against each author

  #--features_to_drop = the number of features to remove from each extreme at the end of each fold, typically 2 or 3.
    #you can also use 'auto' if you want to strip away all variables by the end of your folds

  #--seed = the random seed that you want to use. important for reproducibility

  #--n.iter = number of times to do the whole thing, typically 5 or more (more is usually better, but will take longer)

  #--verbose = do you want the full-detailed output? useful only if you want to "sink" the output to a file for later inspection
 
  #--graphs = do you want to generate some graphs for your "unmasked" observations?

  #--unk.in.results = do you want to include known-authorship items being unmasked against your unknown authorship items? typically you will not want this
  
  #--unmask.unk.only = do you want to run the analysis on *only* the items marked as 'people_to_unmask'? this can be *way* faster, but you won't get the full output data

Results = Unmask(data=data, criterion="Author", people_to_unmask = 'Dubia ', obs.id='Title', num_folds = 10, features_to_drop = 2, seed=666, n.iter=50, verbose=F, graphs=T, unk.in.results=F, unmask.unk.only=F)
#sink()

write.csv(Results, 'X:/writing/Behn MEH/2016-08-09 - LIWC Unmask Results - 50 iterations.csv', row.names=F, na='')




