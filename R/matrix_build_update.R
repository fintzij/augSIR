
# Matrix Construction and Update Functions (build_countmat, build_irm, irm_decomp, update_irm, update_eigen, build_tpm, obs_tpm, tpm_seq, emit_seq, updateX, updateW) --------------------------------------

# build_countmat builds the count matrix containing the the numbers of susceptibles and infecteds at event times
build_countmat <- function(X, popsize) {
    initinfec <- sum(X[X[,1]==0,3])
    Xcount <- rbind(c(0, initinfec), X[X[,1]!=0,c(1,3)])
    Xcount[,2] <- cumsum(Xcount[,2])
    
    numsusc <- popsize - cumsum(c(initinfec, X[X[,1]!=0, 3]>0))
    
    Xcount <- cbind(Xcount, numsusc)
    colnames(Xcount)[2] <- "numsick"
    
    return(Xcount)
    
}

# build_irm constructs an array of instantaneous rate matrices, with the bottom right entry denoting the number of infecteds. The top left (n-1)x(n-1) matrix is the irm.
build_irm <- function(Xcount, b, m, a=0, popsize, pop){
    
    if(pop == FALSE){ # construct irm array for an individual
        
        numinf <- 0:(max(Xcount[,2]) + 1)
        irm <- array(0, dim = c(4,4,length(numinf)))
        irm[1,2,] <- b * numinf + a; irm[1,1,] <- -irm[1,2,]
        irm[2,3,] <- m; irm[2,2,] <- -m
        irm[4,4,] <- numinf
        
    } else if(pop == TRUE){ # construct irm array for the population
        
        numinf <- Xcount[,2]
        numsusc <-  Xcount[,3]
        indend <- dim(Xcount)[1]
        
        irm <- array(0, dim = c(3,3, length(numinf) - 1))
        irm[1,2,] <- (b*numinf[1:(indend-1)] + a)*numsusc[1:(indend-1)]; irm[1,1,] <- -irm[1,2,]
        irm[2,3,] <- m*numinf[1:(indend-1)]; irm[2,2,] <- -irm[2,3,]
        
    }
    
    return(irm)
}

# irm_decomp takes an array of path irms and returns an array of eigenvalues, eigenvectors, and the inverse of the eigenvectors
# irm_decomp returns an array, each element of which is another array containing a diagonal matrix of eigenvalues, a matrix of eigenvectors, and 
# the inverse of the matrix of eigenvectors

irm_decomp <- function(pathirm.cur){
    
    irm.decomp <- array(0,dim = c(3,3,3,dim(pathirm.cur)[3]))
    to_reorder <- pathirm.cur[1,1,] > pathirm.cur[2,2,]
    
    for(s in 1:dim(pathirm.cur)[3]){
        
        decomp <- eigen(pathirm.cur[1:3, 1:3, s], symmetric=FALSE)
        
        if(to_reorder[s] == FALSE){
            irm.decomp[,,1,s] <- diag(decomp[[1]])
            irm.decomp[,,2,s] <- decomp[[2]]
            irm.decomp[,,3,s] <- solve(irm.decomp[,,2,s])
            
        } else if(to_reorder[s] == TRUE){
            irm.decomp[,,1,s] <- diag(decomp[[1]][c(2,1,3)])
            irm.decomp[,,2,s] <- decomp[[2]][,c(2,1,3)]
            irm.decomp[,,3,s] <- solve(irm.decomp[,,2,s])
            
        }
    }
    return(irm.decomp)
    
}

# update_irm updates the array of individual level irms, adding a new irm to the end of the array
update_irm <- function(irm, new.numinf, b, m, a, popsize){
    newmat <- diag(0,4)
    newmat[1,2] <- (b * new.numinf + a); newmat[1,1] <- -newmat[1,2]
    newmat[2,3] <- m; newmat[2,2] <- -m
    newmat[4,4] <- new.numinf
    
    irm.new <- array(c(irm,newmat), dim = c(4,4,dim(irm)[3]+1))
    
    return(irm.new)
}

# update_eigen updates the array of eigen decompositions associated with the array of individual level irms, adding a new array at the end
update_eigen <- function(patheigen, pathirm){
    ind <- dim(pathirm)[3]
    
    decomp <- eigen(pathirm[c(1:3),c(1:3),ind], symmetric = FALSE)
    
    if(pathirm[1,1,ind] > pathirm[2,2,ind]){
        decomp$values <- decomp$values[c(2,1,3)]
        decomp$vectors <- decomp$vectors[,c(2,1,3)]
    }
    
    eigen.new <- array(c(patheigen, 
                         array(c(diag(decomp$values), decomp$vectors, solve(decomp$vectors)), dim = c(3,3,3))), 
                       dim = c(3,3,3,dim(patheigen)[4] + 1))   
    
    
    return(eigen.new)    
}


# buildtpm constructs a transition probability matrix using the matrix exponential
buildtpm <- function(values, vectors, inv.vecs, t0, t1){
    vectors%*%diag(exp(diag(values)*(t1-t0)))%*%inv.vecs
}


# tpm_seq constructs an array of transition probability matrices given a sequence of observation times, a matrix of counts in the form of Xcount, and an array of eigen decompositions for irms
# tpm_seq returns a list of length 2. The first element is an array of tpms for the sequence of observation times. The second is a list, each element of which is a stacked array of tpms, where 
# the subsequence of tpms for the subsequence of event times is located in the [,,,1] element, and the subsequence of tpm products is located in the [,,,2] element. 
tpm_seq <- function(Xcount, obstimes, irm.eig){
    tpms <- array(0, dim = c(3,3,length(obstimes)-1)); tpm.seqs <- list()
    
    eventtimes <- Xcount[,1]; numinf <- Xcount[,2]
    
    for(s in 1:(length(obstimes)-1)){
        
        t0 <- obstimes[s]; t1 <- obstimes[s+1]
        timeseq <- c(t0, eventtimes[(eventtimes > t0) & (eventtimes <t1)], t1)
        
        indstart <- sum(eventtimes <= t0); indend <- indstart + length(timeseq)-2 
        inds <- numinf[indstart:indend] + 1 # indices for indexing into array of eigen decompositions
        subseq.length <- length(inds) # number of matrices in tpm subsequence (= length(timeseq) - 1)
        last.ind <- inds[subseq.length] # eigen array index for last tpm in subsequence
        
        tpm.subseq <- array(0, dim = c(3, 3, subseq.length, 2)) # array for subsequence, containing 3x3 submatrices in two sequences, the first seq is tpms, the second is products
        
        tpm.subseq[,, subseq.length, 1] <- tpm.subseq[,, subseq.length, 2] <- buildtpm(values = irm.eig[,,1,last.ind],
                                                                                       vectors = irm.eig[,,2,last.ind],
                                                                                       inv.vecs = irm.eig[,,3,last.ind],
                                                                                       t0 = timeseq[subseq.length],
                                                                                       t1 = t1)
        if(subseq.length > 1){
            for(t in (subseq.length - 1):1){
                ind <- inds[t]
                tpm <- buildtpm(values = irm.eig[,,1,ind],
                                vectors = irm.eig[,,2,ind],
                                inv.vecs = irm.eig[,,3,ind],
                                t0 = timeseq[t],
                                t1 = timeseq[t+1])
                
                tpm.subseq[,, t, 1] <- tpm
                tpm.subseq[,, t, 2] <- tpm %*% tpm.subseq[,, t+1, 2]
            }
        }
        
        tpms[,,s] <- tpm.subseq[,, 1, 2]
        tpm.seqs[[s]] <- tpm.subseq
        
    }
    
    return(list(tpms, tpm.seqs))
    
}



# emit_seq constructs an array of emission matrices given an observation matrix W, the last column of which is the number of augmented trajectories and the binomial sampling probability
emit_seq <- function(W, p){
    emits <- matrix(0, nrow = 3, ncol = dim(W)[1])
    
    for(s in 1:dim(emits)[2]){
        emits[,s] <- dbinom(x = W[s, 2], size = W[s, 3] + c(0,1,0), prob = p)
        
    }
    
    return(emits)
}

# updateW updates an observation matrix. The function takes as arguments an observation matrix with counts, Xcount, and a matrix, W, 
# with columns containing observation times, number of infecteds observed, and number of infecteds based on augmented trajectories. 
updateW <- function(W, Xcount, path){
    if(missing(path)){
        for(k in 1:dim(W)[1]){
            
            W[k,3] <- Xcount[sum(Xcount[,1]<=W[k,1]),2]
            
        }
        
    } else if(missing(Xcount) & !all(path == 0)){ # Important note: when the path is provided, we should be updating Wother
        W[,3] <- W[,3] + ((W[,1] >= path[1]) & (W[,1] < path[2]))
        
    }
    
    
    return(W)
}

# updateX updates the matix of observations X whose columns are event times, subject id, and event codes. Inputs are the original X matrix,
# a vector of times, and the subject id. The function returns an update X matrix, ordered by event time. 
updateX <- function(X, path, j){
    if(all(path == 0) | all(path == Inf)){
        X[X[,2]==j,c(1,3)] <- 0
        
    }  else if(path[1] != Inf & path[2] ==Inf){
        X[X[,2]==j,][1,1] <- path[1]; X[X[,2]==j,][1,3] <- 1
        X[X[,2]==j,][2,1] <- 0 ; X[X[,2]==j,][2,3] <- 0
        
    } else if(all(path!=Inf)){
        X[X[,2]==j,][1,1] <- path[1]; X[X[,2]==j,][1,3] <- 1
        X[X[,2]==j,][2,1] <- path[2] ; X[X[,2]==j,][2,3] <- -1
    }
    
    X <- X[order(X[,1]),]
    
    return(X)
}

########## Resume here. The mistake is likely in the get_Xcount_other and update_Xcount functions
# update_Xcount updates the Xcount.other matrix with a new sample path
update_Xcount <- function(Xcount.other, path){
    
    if(all(path == 0)){
        Xcount <- Xcount.other
        
        Xcount[,3] <- Xcount[,3] + 1 # subject was always susceptible, so we add him to susceptible count
        
    } else if(path[1] == 0 & path[2] != 0){# subject is initially infected 
        
        if(path[2] != Inf){ # and a recovery is observed
            
            ind <- sum(Xcount.other[,1] <= path[2]) + 1
            Xcount <- insertRow(Xcount.other, c(path[2], Xcount.other[ind - 1, 2:3]), ind)
            
            Xcount[1:(ind-1), 2] <- Xcount[1:(ind -1), 2] + 1 # add subject to count of infecteds for times until his recovery. Subject was never susceptible, so no changes to susceptible count.
            
        } else if(path[2] == Inf){
            
            Xcount <- Xcount.other
            
            Xcount[,2] <- Xcount[,2] + 1 # subject was always infected, so we add him to the count of infecteds for all times
            
        }
        
    } else if(all(path != 0)){ # if subject is not initially infected, and and infection is observed
        
        ind1 <- sum(Xcount.other[,1] <= path[1]) + 1 
        
        Xcount <- insertRow(Xcount.other, c(path[1], Xcount.other[ind1 - 1, 2:3]), ind1)
        
        Xcount[1:(ind1 - 1), 3] <- Xcount[1:(ind1 - 1), 3] + 1 # add subject to the count of susceptibles for the times when he was susceptible
        
        if(path[2] != Inf){ # a recovery is observed
            
            ind2 <- sum(Xcount[,1] <= path[2]) + 1
            Xcount <- insertRow(Xcount, c(path[2], Xcount[ind2 - 1, 2:3]), ind2)
            
            Xcount[ind1:(ind2-1), 2] <- Xcount[ind1:(ind2-1), 2] + 1 # add subject to the count of infecteds for the times he was infected
            
        } else if(path[2] == Inf){ # no recovery is observed
            
            Xcount[ind1:nrow(Xcount), 2] <- Xcount[ind1:nrow(Xcount), 2] + 1 # add subject to the count of infecteds for the times he was infected
        }
        
    }
    
    return(Xcount)
    
}

# get_Wother obtains the observation matrix excluding subject j
get_W_other <- function(W.cur, path){
    Wother <- W.cur
    
    if(!all(path == 0)){
        Wother[,3] <- Wother[,3] - ((Wother[,1] >= path[1]) & (Wother[,1] <= path[2]))
        
    }
    
    return(Wother)
} 

# get_Xcount_other obtains the count matrix excluding subject j
get_Xcount_other <- function(Xcount, path){
    
    if(all(path == 0)){
        Xcount.other <- Xcount
        
        Xcount.other[,3] <- Xcount.other[,3] - 1 # subject was always susceptible, so we remove him from susceptible count
        
    } else if(path[1] == 0 & path[2] != 0){# subject is initially infected and a recovery is observed
        
        Xcount.other <- Xcount[Xcount[,1] != path[2],] # remove the row corresponding to the recovery 
        
        Xcount.other[,2] <- Xcount.other[,2] - (Xcount.other[,1] < path[2]) # remove subject from the infecteds count
        
    } else if(all(path != 0)){ # if subject is not initially infected, and and infection is observed
        
        Xcount.other <- Xcount[!(Xcount[,1] %in% path), ] # remove the rows corresponding to the event times for the current path. note that Inf is never in the event times, so we are not falsely removing an unobserved recovery
        
        Xcount.other[,3] <- Xcount.other[,3] - (Xcount.other[,1] < path[1]) # remove subject j from the susceptibles count for the times when he was susceptible
        Xcount.other[,2] <- Xcount.other[,2] - ((Xcount.other[,1] > path[1]) & (Xcount.other[,1] < path[2])) # remove subject j from the infecteds count when he was infected
    }
    
    return(Xcount.other)    
}
