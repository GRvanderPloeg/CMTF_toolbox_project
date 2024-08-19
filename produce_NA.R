produce_NA = function(data,
                      mechanism = "MCAR", #c("MCAR", "MAR", "MNAR"),
                      perc.missing = 0.5,
                      self.mask = NULL, #c("sym","upper","lower")
                      idx.incomplete = NULL,
                      idx.covariates = NULL,
                      weights.covariates = NULL,
                      by.patterns = FALSE,
                      patterns = NULL,
                      freq.patterns = NULL,
                      weights.patterns = NULL,
                      use.all = FALSE,
                      logit.model = "RIGHT",#c("RIGHT","LEFT","MID","TAIL")
                      seed = NULL)
{

  if (!is.null(seed)){
    set.seed(seed)
  }

  stopifnot((mechanism == "MCAR" & !(is.null(perc.missing))) |
              (mechanism %in% c("MAR", "MNAR")))

  if (is.matrix(data)) {
    data <- data.frame(data)
  }

  if (!is.null(self.mask)){
    self.mask <- tolower(self.mask)
  }

  if (mechanism == "MCAR") {
    return(produce_MCAR(data, perc.missing, idx.incomplete))
  } else {

    # temporary fix to handle factors
    # (transform them to numeric and revert the conversion in the end)
    orig.data <- data
    vars_factor <- colnames(data)[!sapply(data, is.numeric)]
    if (length(vars_factor)==1){
      levels_factor <- list(gdata::mapLevels(x=data[,vars_factor]))
    }
    if (length(vars_factor)>1){
      levels_factor <- sapply(data[,vars_factor], FUN = function(x) gdata::mapLevels(x=x))
    }
    data[,vars_factor] <- sapply(data[,vars_factor], as.integer)
    data <- as.data.frame(data)


    # end temporary fix

    #check if there are non-numeric variables
    # if (any(vapply(data, function(x) !(is.numeric(x) | is.factor(x)), logical(1)))) {
    #  data <- as.data.frame(sapply(data, as.numeric))
    # }

    if (by.patterns==TRUE) {

      if(is.null(patterns)){
        if(sum(is.na(data))==0){
          patterns <- mice::ampute.default.patterns(length(data[1,]))
          weights.patterns <- mice::ampute.default.weights(patterns,mechanism)
          freq.patterns <- mice::ampute.default.freq(patterns)
          #check if there are categorical variables
          if (sum(sapply(data, FUN=is.factor))!= 0){
            for (i in which(sapply(data, FUN=is.factor))){
              half1 <- patterns[,1:i]
              half2 <- patterns[,(i+1):length(patterns[1,])]
              patterns <- cbind(half1, matrix(rep(patterns[,i],times=length(levels(as.factor(data[,i])))-1),ncol=length(levels(as.factor(data[,i])))-1, byrow=FALSE),half2)

              half1 <- weights.patterns[,1:i]
              half2 <- weights.patterns[,(i+1):length(weights.patterns[1,])]
              weights.patterns <- cbind(half1, matrix(rep(weights.patterns[,i],times=length(levels(as.factor(data[,i])))-1),ncol=length(levels(as.factor(data[,i])))-1,byrow=FALSE),half2)
              data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
            }
          }
        }else{

          if (sum(sapply(data, FUN=is.factor))!= 0){
            data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
          }


          #Code taken from miss.compare package, accessed: May 10
          MD_patterns <- mice::md.pattern(data, plot = F)
          data_names <- colnames(data)
          MD_patterns <- MD_patterns[, data_names]
          MD_patterns <- MD_patterns[-c(1, nrow(MD_patterns)), ]
          index <- as.numeric(rownames(MD_patterns)) >= 0.05*sum(apply(is.na(data),1,function(x) sum(x)>=1))
          patterns <- MD_patterns[index, ]

          if (is.null(rownames(patterns))) {
            freq.patterns <- 1
          } else {
            totrows <- as.numeric(rownames(patterns))
            freq.patterns <- totrows/sum(totrows)
          }
          # End of code taken from miss.compare package, accessed: May 10

          if (is.null(dim(patterns))){
            patterns <- matrix(patterns,nrow = 1)
          }
          weights.patterns <- logit.bypatterns(data,patterns,mechanism) #this will deal with categorical variables internally


        }

      }
      else{
        if(is.null(weights.patterns)){
          if(sum(is.na(data))==0){
            weights.patterns  <- mice::ampute.default.weights(patterns,mechanism)

            if (sum(sapply(data, FUN=is.factor))!= 0){

              if(length(patterns[1,])==length(data[1,])){
                for (i in which(sapply(data, FUN=is.factor))){

                  half1 <- patterns[,1:i]
                  half2 <- patterns[,(i+1):length(patterns[1,])]
                  patterns <- cbind(half1, matrix(rep(patterns[,i],times=length(levels(as.factor(data[,i])))-1),ncol =length(levels(as.factor(data[,i])))-1, byrow=FALSE),half2)

                  half1 <- weights.patterns[,1:i]
                  half2 <- weights.patterns[,(i+1):length(weights.patterns[1,])]
                  weights.patterns <- cbind(half1, matrix(rep(weights.patterns[,i],times=length(levels(as.factor(data[,i])))-1),ncol=length(levels(as.factor(data[,i])))-1,byrow=FALSE),half2)
                  data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
                }
              }
              data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
            }

          }else{
            if (sum(sapply(data, FUN=is.factor))!= 0){
              if(length(patterns[1,])==length(data[1,])){
                for (i in which(sapply(data, FUN=is.factor))){

                  half1 <- patterns[,1:i]
                  half2 <- patterns[,(i+1):length(patterns[1,])]
                  patterns <- cbind(half1, matrix(rep(patterns[,i],times=length(levels(as.factor(data[,i])))-1), ncol = length(levels(as.factor(data[,i])))-1, byrow=FALSE),half2)
                }
              }
              data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
            }
            weights.patterns <- logit.bypatterns(data,patterns,mechanism)
          }
        }

        if(is.null(freq.patterns)){
          if(sum(is.na(data))==0){

            freq.patterns <- mice::ampute.default.freq(patterns)

            if (sum(sapply(data, FUN=is.factor))!= 0){
              if(length(patterns[1,])==length(data[1,])){
                for (i in which(sapply(data, FUN=is.factor))){

                  half1 <- patterns[,1:i]
                  half2 <- patterns[,(i+1):length(patterns[1,])]
                  patterns <- cbind(half1, matrix(rep(patterns[,i],times=length(levels(as.factor(data[,i])))-1),ncol = length(levels(as.factor(data[,i])))-1, byrow=FALSE),half2)

                  half1 <- weights.patterns[,1:i]
                  half2 <- weights.patterns[,(i+1):length(weights.patterns[1,])]
                  weights.patterns <- cbind(half1, matrix(rep(weights.patterns[,i],times=length(levels(as.factor(data[,i])))-1),ncol=length(levels(as.factor(data[,i])))-1,byrow=FALSE),half2)
                  data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
                }
              }
              data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
            }
            #perc.missing <- perc.missing / sum(apply(patterns, 1, FUN = function(x) sum(x==0))*freq.patterns)
          }else{

            if (sum(sapply(data, FUN=is.factor))!= 0){
              if(length(patterns[1,])==length(data[1,])){
                for (i in which(sapply(data, FUN=is.factor))){

                  half1 <- patterns[,1:i]
                  half2 <- patterns[,(i+1):length(patterns[1,])]
                  patterns <- cbind(half1, matrix(rep(patterns[,i],times=length(levels(as.factor(data[,i])))-1),ncol = length(levels(as.factor(data[,i])))-1, byrow=FALSE),half2)

                  half1 <- weights.patterns[,1:i]
                  half2 <- weights.patterns[,(i+1):length(weights.patterns[1,])]
                  weights.patterns <- cbind(half1, matrix(rep(weights.patterns[,i],times=length(levels(as.factor(data[,i])))-1),ncol=length(levels(as.factor(data[,i])))-1,byrow=FALSE),half2)
                  data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
                }
              }

              data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))

            }

            MD_patterns <- mice::md.pattern(data, plot = F)
            data_names <- colnames(data)
            MD_patterns <- MD_patterns[, data_names]
            MD_patterns <- MD_patterns[-c(1, nrow(MD_patterns)), ]


            ##copied from stackoverflow (https://stackoverflow.com/questions/41779719/find-rows-of-matrix-which-contain-rows-of-another-matrix)
            ind <- match(MD_patterns, patterns)
            rows <- ind %% nrow(patterns)
            m <- matrix(rows, nrow = nrow(MD_patterns))
            matchRows <- apply(m, 1, duplicated, incomparables = NA)
            patternsrows <- which(colSums(matchRows)==ncol(patterns)-1)

            totrows <- as.numeric(rownames(MD_patterns[patternsrows,]))
            freq.patterns <- totrows/sum(totrows)

            #perc.missing <- perc.missing / sum(apply(patterns, 1, FUN = function(x) sum(x==0))*freq.patterns)
          }
        }

      }

      if (!use.all){
        complete_data <- data[apply(is.na(data),1,function(x) (sum(x)>=1)==0),]

        incomplete_data <- data[apply(is.na(data),1,function(x) (sum(x)>=1)>=1),]
        perc.missing <- (perc.missing - dplyr::if_else(nrow(incomplete_data)==0, 0, mean(is.na(incomplete_data)))*nrow(incomplete_data)/nrow(orig.data)) * nrow(orig.data)/nrow(complete_data)

        idx.patterns.var <- which(apply(patterns, 2, function(x) sum(x==0)>0))
        perc.missing <- perc.missing  / sum(apply(patterns, 1, FUN = function(x) sum(x==0))*freq.patterns) * (length(idx.patterns.var)/ncol(orig.data))

        amputed <- mice::ampute(complete_data, prop = perc.missing,
                                patterns = patterns,freq=freq.patterns,
                                mech = mechanism, weights = weights.patterns,
                                type = logit.model, bycases = FALSE)

        tmp <- amputed$amp
        if (length(vars_factor) > 0){
          for (i in 1:length(vars_factor)){
            data[,vars_factor[[i]]] <- as.factor(data[,vars_factor[[i]]])
            complete_data[,vars_factor[[i]]] <- as.factor(complete_data[,vars_factor[[i]]])
            incomplete_data[,vars_factor[[i]]] <- as.factor(incomplete_data[,vars_factor[[i]]])
            tmp[,vars_factor[[i]]] <- as.factor(tmp[,vars_factor[[i]]])

            gdata::mapLevels(x=data[,vars_factor[[i]]]) <- levels_factor[[i]]
            gdata::mapLevels(x=tmp[,vars_factor[[i]]]) <- levels_factor[[i]]
            gdata::mapLevels(x=complete_data[,vars_factor[[i]]]) <- levels_factor[[i]]
            gdata::mapLevels(x=incomplete_data[,vars_factor[[i]]]) <- levels_factor[[i]]
          }
        }
        idx_newNA <- matrix(rep(FALSE, prod(dim(incomplete_data))), nrow = nrow(incomplete_data), ncol = ncol(incomplete_data))

        idx_newNA <- rbind(idx_newNA, is.na(tmp))
        data.incomp <- rbind(incomplete_data, tmp)

      }
      if (use.all) {
        idx.patterns.var <- which(apply(patterns, 2, function(x) sum(x==0)>0))
        perc.missing <- perc.missing  / sum(apply(patterns, 1, FUN = function(x) sum(x==0))*freq.patterns)

        not.missing <- as.data.frame(imputeMean(data))
        idx_newNA <- matrix(rep(FALSE, prod(dim(orig.data))), nrow = nrow(orig.data), ncol = ncol(orig.data))

        if (sum(sapply(data, FUN=is.factor))!= 0){
          data <- data.frame(mltools::one_hot(data.table::as.data.table(data)))
          not.missing <- data.frame(mltools::one_hot(data.table::as.data.table(not.missing)))
        }

        not.missing <- as.matrix(not.missing)
        amputed <- mice::ampute(not.missing, prop = perc.missing*mean(!is.na(orig.data)),
                                patterns = patterns,freq=freq.patterns,
                                mech = mechanism, weights = weights.patterns,
                                type = logit.model, bycases = FALSE)
        data.incomp <- amputed$amp
        if (length(vars_factor) > 0){
          for (i in 1:length(vars_factor)){
            data[,vars_factor[[i]]] <- as.factor(data[,vars_factor[[i]]])
            data.incomp[,vars_factor[[i]]] <- as.factor(data.incomp[,vars_factor[[i]]])

            gdata::mapLevels(x=data[,vars_factor[[i]]]) <- levels_factor[[i]]
            gdata::mapLevels(x=data.incomp[,vars_factor[[i]]]) <- levels_factor[[i]]
          }
        }

        idx_newNA <- pmax(idx_newNA, as.matrix(is.na(data.incomp)))
        idx_newNA <- apply(idx_newNA, c(1,2), as.logical)
        data.incomp[is.na(data)] <- NA #re-storing original missing data

      }



      return(list("data.init" = data,
                  "data.incomp" = data.incomp,
                  "idx_newNA" = idx_newNA))


    } else{ #by.patterns==FALSE

      #code adapted from mice package, accessed: May 10
      if (!is.null(weights.covariates)) {
        if (is.vector(weights.covariates) && (length(weights.covariates) / ncol(data)) %% 1 == 0) {
          weights.covariates <- matrix(weights.covariates, length(weights.covariates) / ncol(data), byrow = TRUE)
        } else if (is.vector(weights.covariates)) {
          stop("Length of weight vector does not match #variables", call. = FALSE)
        } else if (!is.matrix(weights.covariates) && !is.data.frame(weights.covariates)) {
          stop("Weights matrix should be a matrix", call. = FALSE)
        }
        #end of code adapted from mice package, accessed: May 10


        if(!is.null(idx.covariates) & any(idx.covariates!=(weights.covariates>0))){

          stop("Weights.covariates and idx.covariates must agree", call. = FALSE)

        }

        # # have to make sure the weights.covariates and the covariates selected match each other
        if(is.null(idx.covariates)){
          idx.covariates <- (weights.covariates >= 0)
        }
      }

      #code adapted from mice package, accessed: May 10
      if (!is.null(idx.covariates)) {
        if (is.vector(idx.covariates) && (length(idx.covariates) / ncol(data)) %% 1 == 0) {
          idx.covariates <- matrix(idx.covariates, length(idx.covariates) / ncol(data), byrow = TRUE)
        } else if (is.vector(idx.covariates)) {
          stop("Length of idx.covariates vector does not match #variables", call. = FALSE)
        } else if (!is.matrix(idx.covariates) && !is.data.frame(idx.covariates)) {
          stop("idx.covariates should be a matrix", call. = FALSE)
        }
        #end of code adapted from mice package, accessed: May 10
      }

      if (!(mechanism == "MCAR")){
        if (length(vars_factor) == 0){
          return(produce_MAR_MNAR(data, mechanism, perc.missing, self.mask, idx.incomplete, idx.covariates, weights.covariates, logit.model))
        }
        tmp <- produce_MAR_MNAR(data, mechanism, perc.missing, self.mask, idx.incomplete, idx.covariates, weights.covariates, logit.model)

        if (length(vars_factor) > 0){
          for (i in 1:length(vars_factor)){
            tmp$data.init[,vars_factor[[i]]] <- as.factor(tmp$data.init[,vars_factor[[i]]])
            tmp$data.incomp[,vars_factor[[i]]] <- as.factor(tmp$data.incomp[,vars_factor[[i]]])

            gdata::mapLevels(x=tmp$data.init[,vars_factor[[i]]]) <- levels_factor[[i]]
            gdata::mapLevels(x=tmp$data.incomp[,vars_factor[[i]]]) <- levels_factor[[i]]
          }
        }
        return(list("data.init" = tmp$data.init,
                    "data.incomp" = tmp$data.incomp,
                    "idx_newNA" = tmp$idx_newNA))

      }
    }

  }

}
