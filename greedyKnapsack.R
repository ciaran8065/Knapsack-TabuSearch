tabuSearch3<-function(size = 10, iters = 100, objFunc = NULL, config = NULL,
                      neigh = size, listSize = 9, nRestarts = 10, repeatAll = 1,
                      verbose = FALSE,weights=NULL,values=NULL,limit=NULL)
{
  if (size < 2) { #1 variable isn't enough, need something to compare to.
    stop("error: config too short!")
  }
  if (iters < 2) { #more iterations more likely to find good answer.
    stop("error: not enough iterations!")
  }
  if (listSize >= size) { #You can't have a tabuList with size > the number of possible moves.
    stop("error: listSize too big!")
  }
  if (neigh > size) { #can't check more neighbours than there are
    stop("error: too many neighbours!")
  }
  if (is.null(objFunc)) { #The algorithm requires a user defined objective function
    stop("A evaluation function must be provided. See the objFunc parameter.")
  }
  if(is.null(weights)){
    stop("A vector of weights must be provided")
  }
  if(is.null(values)){
    stop("A vector of values must be provided")
  }
  if(is.null(limit)){
    stop("A maximum weight must be provided")
  }
  if (is.null(config)) { #if no config create a default.
    greedy<-function(w,v){
      vals<-v
      ws<-w
      rs<-vals/ws
      ind<-sort.int(rs,decreasing=T,index.return=T)$ix
      reOrderedWeights<-(ws[ind])
      curW<-0
      Wm<-400
      conf<-numeric(8)
      track<-1
      while(track<=8){
        if(curW+reOrderedWeights[track]>Wm){
          break
        }else{
          curW<-curW+reOrderedWeights[track]
          conf[ind[track]]<-1
          track<-track+1
        }
      }
      return(conf)
    }
    config<-greedy(weights,values)
  }else if (size != length(config)) {
    stop("Length of the starting configuration != size")
  }
  if (repeatAll < 1) {
    stop("error: repeatAll must be > 0")
  }
  iter <- 1 #not to be confused with "iters", iter is used to track the next row/position in the matrices/vectors at which we insert data.
  configKeep <- matrix(, repeatAll * iters * (nRestarts + 3), size) #An empty matrix with rows specified and "size" columns to hold the configuration of binary string at each stage
  eUtilityKeep <- matrix(, repeatAll * iters * (nRestarts + 3),3) #vector of undetermined type, to hold the values of objective function at each itteration
  for (j in 1:repeatAll) { #ALGORITHM START POINT
    if (j > 1) { #if it's not the first iteration through the algorithm we use a random configuration as the user can only specify *initial* config
      config <- matrix(0, 1, size) #a single row matrix of 0's of length "size"
      config[sample(1:size, sample(1:size, 1))] <- 1 #get a random number (from 1:size) of samples from the vector 1:size. Using these as indices set the digits of the config matrix to 1 at these indices.
    }
    
    tabuList <- matrix(0, 1, size) #creating the tabuList as a matrix of 0's with 1 row and "size" columns, it denotes which moves are tabu
    listOrder <- matrix(0, 1, listSize) #This single row matrix holds the order in which certain moves got added to the tabuList, it is used to implement the FIFO method.
    eUtility <- objFunc(config,weights,values,limit) #objective function evaluated at initial condition defined by "config"
    aspiration <- eUtility #initial aspiration value, if the objective function evaluated at a certain configuration is tabu but > this we can cancel it's tabu status
    
    
    preliminarySearch <- function() {
      configKeep[iter, ] <- config #put the current config into the matrix of used configurations in the current row denoted by "iter"
      eUtilityKeep[iter,] <- eUtility #put the value of the current config into the first position of the vector
      
      iter <- iter + 1 #increment the iteration, if it's the first access of this function it is set to 2 because we do 1 outside where we initially get eUtility etc.
      for (i in 2:iters) { #see comment above to explanation of 2:iters
        neighboursEUtility <- matrix(0, size, 3)#1 row matrix to hold the value of each of the neighbours' configurations evaluated at the objective function
        #print(neighboursEUtility)
        #print(config)
        configTemp <- t(matrix(config, size, neigh)) #a matrix where each row is set to config, transpose used as matricies by default fill by column
        #print(configTemp)
        randomNeighbours <- sample(size, neigh) #if neigh<size it chooses the random neighbours as required but if neigh=size it takes all of them
        #print(randomNeighbours)
        diag(configTemp[, randomNeighbours]) <- abs(diag(configTemp[, randomNeighbours]) - 1) #the diagonal of the matrix is inverted so that each row is different by 1 bit, this gives all the possible 1 bit-flip neighbour configurations
        #print(configTemp)
        neighboursEUtility[randomNeighbours,] <- t(apply(configTemp, 1, objFunc,weights=weights,values=values,limit=limit))#"neighboursEUtility" will hold the value of the objective function calculated at each configuration
        #print(neighboursEUtility)
        maxNontaboo <- max(neighboursEUtility[tabuList == 0,1]) #finding the maximum non-tabu value
        #print(maxNontaboo)
        maxTaboo <- max(neighboursEUtility[tabuList == 1,1], 0) #finding the maximum tabu value
        #print(maxTaboo)
        #print(aspiration)
        #here we get that maxTaboo will be 0 which works for the normal ones because 0 would be ignored but here when 0 is good it causes problems.
        
        #determining where to move to next
        #ERROR IN THIS SECTION
        
        move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration[1],  #condition: if the move is tabu AND it's higher than the asp value
                       
                       ifelse(length(which(neighboursEUtility[,1] == maxTaboo)) == 1, #if true: if the maxTabu from neighbour values is unique
                              which(neighboursEUtility[,1] == maxTaboo), #take the value of the position of that neighbour (i.e the next place to move to)
                              sample(which(neighboursEUtility[,1] == maxTaboo), 1)), #if it is not unique take a random one (i.e ties can be broken arbitrarily)
                       #so here is the issue ^^^^^^^^^^^^ we get into here because of the issue with 0 mentioned earlier then we have a sample without specifying what to sample from hence the sampling error.
                       ifelse(length(which(neighboursEUtility[,1] == maxNontaboo & tabuList == 0)) == 1, #if false: (either maxTabu is < maxNonTabu OR maxTabu < aspiration) AND it's unique
                              which(neighboursEUtility[,1] == maxNontaboo & tabuList == 0), #take the value of the position of the maxNonTabu
                              sample(which(neighboursEUtility[,1] == maxNontaboo & tabuList == 0), 1)) #if it's not unique pick a random one
                       
        ) #NEXT MOVE HAS BEEN CHOSEN
        #ERROR IN THIS SECTION
        #print(move)
        if (eUtility[1] >= neighboursEUtility[move,1]) {#if the initial(or previous) evaluated configuration is >= the value of the neighbour config we move to it but DO NOT update aspiration value
          #print("first if")
          tabuList[move] <- 1 #set that move to tabu
          if (sum(tabuList) > listSize) { #in this the default is tabuList in which 10 different moves can be tabu but at max there can only be 9 at one time so that we can always move
            #print("inner first if")
            tabuList[listOrder[1]] <- 0 #set the first move which was made tabu to nonTabu (uses first in first out system)
            listOrder[1:listSize] <- c(listOrder[2:listSize], 0) #discards the first element and moves all the other elements 1 to the left appending a 0 to the end. (like a left shift operation)
          }
          listOrder[min(which(listOrder == 0))] <- move #add the next move found earlier to the first available position in "listOrder"
        }
        else if (neighboursEUtility[move,1] > aspiration[1]){ #if the new move results in a configuration which when evaluated at the objective function is > aspiration value, set the new aspiration value to this value
          #print("else if")
          aspiration <- neighboursEUtility[move,]
        }
        eUtility <- neighboursEUtility[move,] #set the max found value to the value found at this neighbour
        
        config[move] <- abs(config[move] - 1) #update the configuration by switching the bit at the position corresponding to the move
        configKeep[iter, ] <- config #put the configuration used into "configKeep" which is storing the configurations
        eUtilityKeep[iter,] <- eUtility #put the value corresponding to this configuration into "eUtilityKeep"
        iter <- iter + 1 #increment iter
      }
      result = list(aspiration = aspiration, configKeep = configKeep, eUtilityKeep = eUtilityKeep, iter = iter) #putting together the result of the preliminary search
      return(result)
    }
    
    
    if (verbose) #if the user wants to see the steps as the algorithm executes
      cat("Preliminary search stage...\n")
    result <- preliminarySearch() #store result of prelim search, which is a list containing values visible in line 89
    aspiration <- result$aspiration[1] #set the new aspiration value to the value obtained from the preliminary search
    configKeep <- result$configKeep #set the matrix of configurations to the configs gotten from the preliminary search function
    eUtilityKeep <- result$eUtilityKeep #set the vector containing the values of the objective function evaluated at each of the configurations to the values obtained from the preliminary search function
    iter <- result$iter #extract from the result list the number of iterations that occured
    temp_asp <- 0 #temporary value for aspiration
    restarts <- 0 #tracks the number of restarts of the intensification stage
    #print(temp_asp < aspiration[1] & restarts < nRestarts)
    while (temp_asp < aspiration[1] && restarts < nRestarts) { #check if an equal or better solution has been found (when temp_asp >=aspiration) or we run out of restarts
      if (verbose)
        cat("Intensification stage...\n")
      eUtility <- max(eUtilityKeep[,1],na.rm=TRUE) #intensification stage begins with the best solution found thus far
      #print(eUtility)
      temp_asp <- aspiration #assigning the value of the current aspiration to temp_asp
      #print(temp_asp)
      #print(configKeep)
      #print(max(eUtilityKeep[,1]))
      #print(eUtilityKeep[,1] == max(eUtilityKeep[,1]))
      #print(is.na(eUtilityKeep))
      na_rm<-as.matrix(eUtilityKeep[!is.na(eUtilityKeep)],iters,3)
      #print(na_rm)
      config <- configKeep[max(which(na_rm[,1] == max(na_rm[,1]))), ] #takes the configuration from "configKeep" where the position is the position of the most recent, best configuration.
      #print(config)
      result <- preliminarySearch() #re-run the preliminary search
      aspiration <- result$aspiration #obtain results
      configKeep <- result$configKeep
      eUtilityKeep <- result$eUtilityKeep
      iter <- result$iter
      restarts <- restarts + 1 #increment the number of restarts
    }
    if (verbose)
      cat("Diversification stage...\n")
    config <- matrix(0, 1, size) #reset the config matrix
    config[sample(1:size, sample(1:size, 1))] <- 1 #take a random starting config
    eUtility <- objFunc(config,weights,values,limit) #evaluate the objective function at this starting config
    #Now finding the most frequent moves and sets them to tabu
    frequent <- apply(configKeep, 2, function(x) sum(diff(x) != 0)) #For each column in configKeep, "sum(diff(x) !=0)" will count how many times this digit of the config changes.
    tabuList <- as.numeric(rank(frequent, ties.method = "random") > (size - listSize)) #the values in frequent are ranked in order of magnitude, these are then compared to "(size-listSize)" by default this is 1. If frequent is > this 1, that move is set to tabu.
    listOrder <- sample(which(rank(frequent, ties.method = "random") > (size - listSize)), listSize) #A random sample of size "listSize" is taken to be the order of the tabu list, we need an order to use the First in First out method.
    result <- preliminarySearch() #re-run preliminary search in the hopes to find a better solution in a reasonably unexplored area of the solution space
    iter <- result$iter #obtain results
    configKeep <- result$configKeep
    eUtilityKeep <- result$eUtilityKeep
    
  }#end of algorithm, proceed past this point when j>repeatAll
  
  #putting together the end results
  endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ], eUtilityKeep = eUtilityKeep[1:(iter - 1),], iters = iters, neigh = neigh, listSize = listSize, repeatAll = repeatAll,weights=weights,values=values,limit=limit)
  class(endResult) = "tabu"
  return(endResult)
  
}

vals<-c(50,40,30,60,100,150,120,70)
ws<-c(50,30,60,50,50,50,50,200)

eval<-function(conf,weights,values,limit){
  Wm<-limit
  vals<-values
  ws<-weights
  cVal<-sum(conf*vals) #config value
  cW<-sum(conf*ws) #config weight
  lim<-Wm+mean(ws)
  #lim<-Wm+lightest-1 try the -1 so it may handle cases better where everything is the same weigh
  res<-numeric(3)
  res[1]<-ifelse(cW>lim,0,cVal) #if the weight is over the limit return 0, else return the value
  res[2]<-cW
  res[3]<-ifelse(cW>Wm && cW<lim,1,0) #if the weight is greater than Wm AND less than the limit return 1 else return 0
  return(res)
}

summ<-function (res)
{
  v<-res$eUtilityKeep[,3]==0
  exc<-which(res$eUtilityKeep[,3]!=0) #where it is not 0, we can exclude these so that the which.max correctly matches up
  fixed<-res$configKeep[-exc,]
  finalConfig<-fixed[which.max((res$eUtilityKeep[res$eUtilityKeep[,3]==0,])[,1]),]
  ans<-eval(finalConfig,res$weights,res$values,res$limit)
  values[i]<-ans[1]
  cat(" Value found: ",ans[1],"\n","Weight used: ",ans[2],"\n")
  print(rbind(finalConfig))
}

res<-tabuSearch3(size=8,iters=50,objFunc=eval,listSize=4,nRestarts=10,weights=ws,values=vals,limit=400)
summ(res)

