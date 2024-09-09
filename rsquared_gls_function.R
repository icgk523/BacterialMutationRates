#' Remove random effects from all.vars
#' 
#' @keywords internal
#' 
all_vars_merMod <- function(formula.) {
  
  if(!any(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)
  
  if(inherits(formula., "formula.cerror"))
    
    gsub(" " , "", unlist(strsplit(formula., "~~"))) else {
      
      n <- rownames(attr(terms(formula.), "factors"))
      
      if(any(grepl("\\|", n)))
        
        all.vars(lme4::nobars(formula.)) else
          
          all.vars(formula.)
      
    }
  
}

#' Get vector of untransformed variables
#' 
#' @keywords internal
#' 
all_vars_notrans <- function(formula.) {
  
  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)
  
  if(inherits(formula., "formula")) {
    
    if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)
    
    formula. <- all_vars_trans(formula.)
    
    if(any(grepl(":", formula.))) {
      
      idx <- which(grepl(":", formula.))
      
      for(i in idx) formula.[i] <- paste(sapply(strsplit(formula.[i], ":"), stripTransformations), collapse = ":")
      
      for(j in (1:length(formula.))[-idx]) formula.[j] <- stripTransformations(formula.[j])
      
    } else {
      
      formula. <- sapply(formula., stripTransformations)
      
    }
    
  } else formula. <- unlist(strsplit(formula., " ~~ "))
  
  ret <- gsub("(,.*)", "", formula.)
  
  return(ret)
  
}

#' Get vector of transformed variables
#' 
#' @keywords internal
#' 
all_vars_trans <- function(formula., smoothed = FALSE) {
  
  if(!all(class(formula.) %in% c("formula", "formula.cerror"))) formula. <- formula(formula.)
  
  if(inherits(formula., "formula")) {
    
    if(formula.[[3]] == 1) ret <- deparse(formula.[[2]]) else {
      
      if(any(grepl("\\|", formula.))) formula. <- lme4::nobars(formula.)
      
      ret <- c(rownames(attr(terms(formula.), "factors"))[1], labels(terms(formula.)))
      
      if(smoothed == FALSE) ret <- gsub("(.*)\\,.*", "\\1", gsub("s\\((.*)\\).*", "\\1", ret)) 
      
      # else {
      
      # ret <- gsub("(s\\(.*),.*", "\\1", ret)
      # 
      # if(any(grepl("s\\(", ret))) ret <- sapply(ret, function(x) 
      #   ifelse(grepl("s\\(", x) & !grepl("\\)", x), paste0(x, ")"), x))
      
      # }
      
      # ret <- gsub("(,.*)", "", ret)
      
    }
    
    return(ret)
    
  } else unlist(strsplit(formula., " ~~ "))
  
}

#' Captures output table
#' 
#' @keywords internal
#' 
captureTable <- function(g, row.names = FALSE) {
  
  g1 <- capture.output(print(g, row.names = row.names))
  
  if(all(g1 == "data frame with 0 columns and 0 rows")) 
    
    g1 <- "No independence claims present. Tests of directed separation not possible."
  
  g1 <- paste0(g1, "\n")
  
  return(g1)
  
}

#' Bind data.frames of differing dimensions
#'
#' From: https://stackoverflow.com/a/31678079
#' 
#' @param ... data.frames to be bound, separated by commas
#' 
#'  @keywords internal
#'   
cbind_fill <- function(...) {
  
  nm <- list(...) 
  
  dfdetect <- grepl("data.frame|matrix", unlist(lapply(nm, function(cl) paste(class(cl), collapse = " ") )))
  
  vec <- data.frame(nm[!dfdetect])
  
  n <- max(sapply(nm[dfdetect], nrow)) 
  
  vec <- data.frame(lapply(vec, function(x) rep(x, n)))
  
  if (nrow(vec) > 0) nm <- c(nm[dfdetect], list(vec))
  
  nm <- lapply(nm, as.data.frame)
  
  do.call(cbind, lapply(nm, function (df1) 
    
    rbind(df1, as.data.frame(matrix(NA, ncol = ncol(df1), nrow = n-nrow(df1), dimnames = list(NULL, names(df1))))) )) 
  
}

#' Transform variables based on model formula and store in new data frame
#' 
#' @keywords internal
#' 
dataTrans <- function(formula., data) {
  
  notrans <- all_vars_merMod(formula.)
  
  if(inherits(formula., "formula.cerror")) notrans <- gsub(".*\\((.*)\\)", "\\1", notrans)
  
  trans <- all_vars_trans(formula.)
  
  trans <- unlist(strsplit(trans, "\\:"))
  
  trans <- trans[!duplicated(trans)]
  
  if(any(grepl("scale\\(.*\\)", trans))) {
    # 
    # trans[which(grepl("scale(.*)", trans))] <- notrans[which(grepl("scale(.*)", trans))]
    # 
    warning("`scale` applied directly to variable. Use argument `standardize = TRUE` instead.", call. = FALSE)
    
  }
  
  if(any(!notrans %in% trans)) {
    
    for(k in 1:length(notrans)) {
      
      if(is.factor(data[, notrans[k]])) next else 
        
        if(grepl("scale(.*)", trans[k])) data[, notrans[k]] <- scale(data[, notrans[k]]) else
          
          data[, notrans[k]] <-
            
            sapply(data[, notrans[k]], function(x) eval(parse(text = gsub(notrans[k], x, trans[k]))))
        
    }
    
  }
  
  colnames(data) <- notrans
  
  return(data)
  
}

#' Get ANOVA results from `merMod`
#'
#' @keywords internal
#'
getAnova <- function(model, test.statistic = "F", test.type = "III") {
  
  if(inherits(model, "glmmTMB")) test.statistic = "Chisq"
  
  krp <- as.data.frame(car::Anova(model, test.statistic = test.statistic, type = test.type))
  
  ct <- summary(model)$coefficients
  
  colnames(ct)[2] <- "Std.Error"
  
  ret <- do.call(rbind, lapply(1:nrow(krp), function(i) {
    
    if(rownames(krp)[i] %in% rownames(ct)) {
      
      cbind.data.frame(
        ct[i, 1:2, drop = FALSE],
        DF = krp[i, 3], 
        Crit.Value = krp[i, 1], 
        P = krp[i, ncol(krp)],
        row.names = NULL
      )
      
    } else {
      
      data.frame(
        Estimate = NA,
        Std.Error = NA,
        DF = krp[i, 3], 
        Crit.Value = krp[i, 1], 
        P = krp[i, ncol(krp)]
      )
      
    }
    
  } ) )
  
  # ret <- cbind.data.frame(
  #   ct[, 1:2], 
  #   DF = krp[, 3], 
  #   Crit.Value = krp[, 1], 
  #   P = krp[, ncol(krp)]
  # )
  
  names(ret)[ncol(ret)] <- "Pr(>|t|)"
  
  rownames(ret) <- rownames(krp)
  
  return(ret)
  
}

#' Get random effects from lme
#' 
#' @keywords internal
#' 
findbars.lme <- function(model) {
  
  rand <- model$call$random
  
  sapply(rand, function(i) {
    
    i = gsub(".*\\|(.*)", "\\1", as.character(i)[2])
    
    strsplit(gsub(" " , "", i), "\\/")[[1]]
    
  } )
  
}

#' Get data from model list
#' 
#' @keywords internal
#' 
GetData <- function(modelList) {
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)
  
  modelList <- removeData(modelList, formulas = 1)
  
  data.list <- lapply(modelList, GetSingleData)
  
  data.list <- data.list[!sapply(data.list, is.null)]
  
  data.list <- unname(data.list)
  
  if(all(sapply(data.list, class) == "comparative.data"))
    
    data <- data.list[[1]] else 
      
      data <- do.call(cbind_fill, data.list)
  
  data <- data[, !duplicated(colnames(data), fromLast = TRUE)]
  
  # colnames(data) <- gsub(".*\\((.*)\\).*", "\\1", colnames(data))
  
  data <- as.data.frame(data)
  
  rownames(data) <- 1:nrow(data)
  
  return(data)
  
}

#' Get data from one model
#' 
#' @keywords internal
#' 
GetSingleData <- function(model) {
  
  dat <- data.frame()
  
  switch(class(model)[1],
         "phylolm" = {
           stop("Please provide `data =` argument to `psem`.", call. = FALSE)
         },
         
         "phyloglm" = {
           stop("Please provide `data =` argument to `psem`.", call. = FALSE)
         },
         
         "lm" ={
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },
         
         "negbin" = {
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },
         
         "Sarlm" = {
           dat <- eval(getCall(model)$data, environment(formula(model)))
         },
         
         "glm" = {
           dat <- model$data
         },
         
         "glmmPQL" = {
           dat <- model$data
         },
         
         "pgls" = {
           dat <- model$data
         },
         
         "lmerMod" = {
           dat <- lme4::getData(model) #model@frame
         },
         
         "glmerMod" = {
           dat <- lme4::getData(model) #model@frame
         },
         
         "lmerModLmerTest" = {
           dat <- lme4::getData(model) #model@frame
         },
         
         "glmmTMB" = {
           dat <- model$frame 
         },
         
         "gls" = {
           dat <- nlme::getData(model)
         },
         
         "lme" = {
           dat <- nlme::getData(model)
         },
         "gam" = {
           dat <- eval(getCall(model)$data, environment(formula(model)))
         }
         
  )
  
  return(dat)
  
}

#' Obtain (observation-level) random effects from a generalized linear mixed model
#' 
#' RE = "all" all random effects are reported
#' RE = "RE" just group effects are reported
#' RE = "OLRE" just observation-level effects are reported
#' 
#' @keywords internal
#' 
GetOLRE <- function(sigma, model, X, data, RE = c("all", "RE", "OLRE")) {
  
  if(class(model) %in% c("lmerMod", "glmerMod")) {
    
    if(is.null(X)) X <- model.matrix(model)
    
    rand <- sapply(lme4::findbars(formula(model)), function(x) as.character(x)[3])
    
    rand <- rand[!duplicated(rand)] 
    
  } 
  
  # else if(class(model) %in% c("lme", "glmmPQL")) { }
  
  idx <- sapply(sapply(strsplit(rand, "\\:"), function(x) gsub("\\(|\\)", "", x)), function(x) {
    
    length(unique(data[, x])) == nrow(data)
    
  } )
  
  sigma.names <- unlist(names(sigma)) # unlist(strsplit(names(sigma), "\\."))
  
  idx. <- sapply(sigma.names, function(x) !any(x %in% rand[idx]))
  
  if(RE == "RE") 
    
    out <- sapply(sigma[idx.], function(i) {
      
      if(all(rownames(i) %in% colnames(X))) X. <- X else
        
        X. <- do.call(cbind, model.matrix(model, type = "randomListRaw")) 
      
      Z <- as.matrix(X.[, rownames(i), drop = FALSE])
      
      sum(rowSums(Z %*% i) * Z) / nrow(X.)
      
    } ) else if(RE == "OLRE") {
      
      if(all(idx == FALSE)) out <- 0 else {
        
        out <- sapply(sigma[idx], function(i) {
          
          Z <- as.matrix(X[, rownames(i), drop = FALSE])
          
          sum(rowSums(Z %*% i) * Z) / nrow(X)
          
        } ) } } else if(RE == "all")
          
          out <- sapply(sigma, function(i) {
            
            Z <- as.matrix(X[, rownames(i), drop = FALSE])
            
            sum(rowSums(Z %*% i) * Z) / nrow(X)
            
          } )
  
  if(length(out) == 0) out <- 0
  
  return(out)
  
}

#' Get random effects variance-covariance from lme
#' 
#' @keywords internal
#' 
GetVarCov <- function(model) {
  
  vc <- try(getVarCov(model), silent = TRUE)
  
  if(any(class(vc) == "try-error")) {
    
    vc <- nlme::VarCorr(model)
    
    v <- suppressWarnings(as.numeric(vc[, 1]))
    
    names(v) <- gsub(" =", "", rownames(vc))
    
    vm <- as.list(na.omit(v[-length(v)]))
    
    vl <- lapply(1:length(vm), function(i) matrix(vm[[i]], dimnames = list(names(vm)[i], names(vm)[i])))
    
    names(vl) <- names(which(is.na(v)))
    
    vl
    
  } else list(vc)
  
}

#' Assess significance
#' 
#' @keywords internal
#' 
isSig <- function(p) {
  
  ifelse(p > 0.01 & p < 0.05, "*",
         ifelse(p > 0.001 & p <= 0.01, "**",
                ifelse(p <= 0.001, "***", "")))
  
}

#' Recompute P-values using Kenward-Rogers approximation
#' 
#' @keywords internal
#' 
# KRp <- function(model, vars, data, intercepts = FALSE) {
# 
#   # if(any(grepl("\\*", all_vars_notrans(formula(model)))) & !all(grepl("\\*", vars))) {
# 
#     f <- all_vars_trans(formula(model))
# 
#     model <- update(model, as.formula(paste(f[1], " ~ ", paste(f[-1], collapse = " + "), " + ", paste(onlyBars(formula(model)), collapse = " + "))))
# 
#   # }
# 
#   out <- data.frame()
#   
#   for(x in vars) { #sapply(vars, function(x) {
# 
#     reduceModel <- update(model, as.formula(paste(". ~ . -", x)))
# 
#     if(nobs(model) != nobs(reduceModel)) stop("Different sample sizes for `KRmodcomp`. Remove all NAs and re-run")
#     
#     kr <- try(pbkrtest::KRmodcomp(model, reduceModel), silent = TRUE)
# 
#     if(class(kr) == "try-error") 
#       
#       stop("Cannot obtain P-values from `lmerMod` using `pbkrtest::KRmodcopm`. Consider fitting using `nlme::lme`") else {
#         
#         d <- round(kr$stats$ddf, 2)
#   
#         p <- kr$stats$p.valueU
#   
#         out <- rbind(out, data.frame(d, p))
#         
#       }
# 
#   } # )
# 
#   if(intercepts == TRUE) {
# 
#     reduceModelI <- update(model, as.formula(paste("~ . - 1")), data = data)
# 
#     krI <- try(pbkrtest::KRmodcomp(model, reduceModelI), silent = TRUE)
#     
#     if(class(krI) == "try-error") 
#       
#       stop("Cannot obtain P-values from `lmerMod` using `pbkrtest::KRmodcomp`. Consider re-fitting using `nlme::lme`")else {
#         
#         dI <- krI$stats$ddf
#     
#         pI <- krI$stats$p.valueU
#     
#         out <- rbind(data.frame(d = dI, p = pI), out)
#         
#       }
#       
#   }
#   
#   return(out)
# 
# }

#' Get list of formula from a `psem` object
#' 
#' @keywords internal
#' 
listFormula <- function(modelList, formulas = 0) {
  
  modelList <- removeData(modelList, formulas)
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList <- list(modelList)
  
  fList <- lapply(modelList, function(i) if(any(class(i) %in% c("formula.cerror"))) i else formula(i) )
  
  fList <- lapply(fList, lme4::nobars)
  
  return(fList)
  
}

#' Get number of observations from a model
#' 
#' @keywords internal
#' 
nObs <- function(object, ...) if(any(class(object) %in% c("phylolm", "phyloglm", "Sarlm"))) length(fitted(object)) else nobs(object, ...)

#' Get random effects from merMod
#' 
#' @keywords internal
#' 
onlyBars <- function(formula., slopes = TRUE) {
  
  f <- lme4::findbars(formula.)
  
  if(slopes == TRUE) paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ") else {
    
    # paste(sapply(f, function(x) paste0("(1 ", gsub(".*(\\|.*)", "\\1", f), ")")), collapse = "+")
    
    f <- f[sapply(f, function(x) grepl("1\\||1 \\|", deparse(x)))]
    
    paste(sapply(f, function(x) paste0("(", deparse(x), ")")), collapse = " + ")
    
  }
  
}

#' Do not print attributes with custom functions
#' 
#' @keywords internal
#' 
#' @method print attr
#' 
print.attr <- function(x, ...) {
  
  attributes(x) <- NULL
  
  noquote(x)
  
}

#' Remove data from the model list
#' 
#' formulas = 0, keep everything
#' formulas = 1, remove all formulas including correlated errors
#' formulas = 2, remove only formula but keep correlated errors
#' formulas = 3, remove correlated errors but keep formula
#' 
#' @keywords internal
#' 
removeData <- function(modelList, formulas = 0) {
  
  remove <- c("character", "matrix", "data.frame", "SpatialPointsDataFrame", "comparative.data")
  
  if(formulas == 1) remove <- c(remove, "formula", "formula.cerror")
  
  if(formulas == 2) remove <- c(remove, "formula")
  
  if(formulas == 3) remove <- c(remove, "formula.cerror")
  
  modelList[!sapply(modelList, function(x) any(class(x) %in% remove))]
  
}

#' Strip transformations
#' 
#' @keywords internal
#' 
stripTransformations <- function(x) {
  
  x <- gsub(".*\\((.*)\\).*", "\\1", x)
  
  gsub(" ", "", gsub("(.*)\\+.*", "\\1", x))
  
}

#' Get Response Name as a Character
#' 
#' @keywords internal
#' 
get_response <- function(mod) {
  
  mod <- removeData(mod)
  
  f <- lapply(mod, formula)
  
  r <- lapply(f, function(x) x[[2]])
  
  return(as.character(r))
  
}

#' Get Left-hand side of formulae
#' 
#' @keywords internal
#' 
getLHS <- function(formulaList){
  sapply(formulaList, function(x) as.character(x[[2]]))
}

#' Get Right-hand side of formulae
#' 
#' @keywords internal
#' 
getRHS <- function(formulaList){
  
  rhs <- sapply(formulaList, function(x) all.vars(x)[-1])
  
  unique(do.call(c, rhs))
  
}

#' Operator for non-overlap in sets
#' 
#' @keywords internal
#' 
"%not_in%" <- function(x, y) x[!x %in% y]


#' Get a sorted psem object in DAG order
#' 
#' @description Takes a [psem] object, pulls out the
#' DAG, and then sorts the psem object into the order
#' of the DAG (from exogenous to terminal endogenous
#' variable) for use by other functions. Note: removes
#' correlated errors.
#'
#' @param object A fit [psem] object
#' @param keepdata Defaults to TRUE. Should the
#' data with the psem be included in the returned
#' object?
#'
#' @return A new [psem] object, without the data.
#' @export
getSortedPsem <- function(object, keepdata = TRUE){
  #first, remove data
  dat <- object$data
  object <- removeData(object, formulas = 1)
  
  #Now, get formulae
  formulaList <- listFormula(object)
  lhs <- getLHS(formulaList)
  
  names(object)<- lhs
  
  #sort some dags so we do things in the right order  
  object_dag <- getDAG(formulaList)
  sorted_dag <- sortDag(object_dag, formulaList)
  lhs_sorted <- colnames(sorted_dag)
  lhs_sorted <- lhs_sorted[which(lhs_sorted %in% lhs)]
  
  #Sort the object
  object <- object[lhs_sorted]
  
  #should we include the data?
  if(keepdata) object$data <- dat
  
  #return
  return(object)
}






#' R-squared for linear regression
#'
#' Returns (pseudo)-R^2 values for all linear, generalized linear, and
#' generalized linear mixed effects models.
#'
#' For mixed models, marginal R2 considers only the variance by the fixed
#' effects, and the conditional R2 by both the fixed and random effects.
#' 
#' For generalized additive models fit to gaussian distribution, the function
#' returns ths adjusted-R2. For all other distributions, it returns the proportion
#' of deviance explained.
#'
#' For GLMs (\code{glm}), supported methods include: \itemize{
#' \item\code{mcfadden} 1 - ratio of likelihoods of full vs. null models
#'
#' \item\code{coxsnell} McFadden's R2 but raised to 2/N. Upper limit is < 1
#'
#' \item\code{nagelkerke} Adjusts Cox-Snell R2 so that upper limit = 1. The
#' DEFAULT method
#'
#' } For GLMERs fit to Poisson, Gamma, and negative binomial distributions
#' (\code{glmer}, \code{glmmPQL}, \code{glmer.nb}), supported methods include
#' \itemize{ \item\code{delta} Approximates the observation variance based on
#' second-order Taylor series expansion. Can be used with many families and
#' link functions
#'
#' \item\code{lognormal} Observation variance is the variance of the log-normal
#' distribution
#'
#' \item\code{trigamma} Provides most accurate estimate of the observation
#' variance but is limited to only the log link. The DEFAULT method
#'
#' } For GLMERs fit to the binomial distribution (\code{glmer},
#' \code{glmmPQL}), supported methods include: \itemize{
#' \item\code{theoretical} Assumes observation variance is pi^2/3
#'
#' \item\code{delta} Approximates the observation variance as above. The
#' DEFAULT method
#'
#' }
#'
#' @param modelList a regression, or a list of structural equations.
#' @param method The method used to compute the R2 value (See Details)
#' 
#' @return Returns a \code{data.frame} with the response, its family and link,
#' the method used to estimate R2, and the R2 value itself. Mixed models also
#' return marginal and conditional R2 values.
#' @author Jon Lefcheck <lefcheckj@@si.edu>
#' @references Nakagawa, Shinichi, Paul CD Johnson, and Holger Schielzeth. 
#' "The coefficient of determination R 2 and intra-class correlation coefficient 
#' from generalized linear mixed-effects models revisited and expanded." 
#' Journal of the Royal Society Interface 14.134 (2017): 20170213.
#' @examples
#'
#'   \dontrun{
#'     # Create data
#'     dat <- data.frame(
#'       ynorm = rnorm(100),
#'       ypois = rpois(100, 100),
#'       x1 = rnorm(100),
#'       random = letters[1:5]
#'     )
#'
#'     # Get R2 for linear model
#'     rsquared(lm(ynorm ~ x1, dat))
#'
#'     # Get R2 for generalized linear model
#'     rsquared(glm(ypois ~ x1, "poisson", dat))
#'
#'     rsquared(glm(ypois ~ x1, "poisson", dat), method = "mcfadden") # McFadden R2
#'
#'     # Get R2 for generalized least-squares model
#'     rsquared(gls(ynorm ~ x1, dat))
#'
#'     # Get R2 for linear mixed effects model (nlme)
#'     rsquared(nlme::lme(ynorm ~ x1, random = ~ 1 | random, dat))
#'
#'     # Get R2 for linear mixed effects model (lme4)
#'     rsquared(lme4::lmer(ynorm ~ x1 + (1 | random), dat))
#'
#'     # Get R2 for generalized linear mixed effects model (lme4)
#'     rsquared(lme4::glmer(ypois ~ x1 + (1 | random), family = poisson, dat))
#'
#'     rsquared(lme4::glmer(ypois ~ x1 + (1 | random), family = poisson, dat), method = "delta")
#'
#'     # Get R2 for generalized linear mixed effects model (glmmPQL)
#'     rsquared(MASS::glmmPQL(ypois ~ x1, random = ~ 1 | random, family = poisson, dat))
#'     
#'     # Get R2 for generalized additive models (gam)
#'     rsquared(mgcv::gam(ynorm ~ x1, dat))
#'   }
#'
#' @export
#' 
rsquared <- function(modelList, method = NULL) {
  
  if(!all(class(modelList) %in% c("psem", "list"))) modelList = list(modelList)
  
  modelList <- removeData(modelList, formulas = 1)
  
  evaluateClasses(modelList)
  
  ret <- lapply(modelList, function(i) {
    
    if(all(class(i) %in% c("lm", "pgls"))) r <- rsquared.lm(i) else
      
      if(all(class(i) %in% c("gls"))) r <- rsquared.gls(i) else
        
        if(all(class(i) %in% c("glm", "lm", "negbin"))) r <- rsquared.glm(i, method) else
          
          # if(any(class(i) %in% c("phylolm", "phyloglm"))) r <- rsquared.phylolm(i)
          
          if(all(class(i) %in% c("lme"))) r <- rsquared.lme(i) else
            
            if(all(class(i) %in% c("lmerMod", "merModLmerTest", "lmerModLmerTest"))) r <- rsquared.merMod(i) else
              
              if(any(class(i) %in% c("glmerMod"))) r <- rsquared.glmerMod(i, method) else
                
                if(any(class(i) %in% c("glmmTMB"))) r <- rsquared.glmmTMB(i) else
                  
                  if(any(class(i) %in% c("glmmPQL"))) r <- rsquared.glmmPQL(i, method) else
                    
                    if(any(class(i) %in% c("Sarlm"))) r <- rsquared.Sarlm(i) else
                      
                      if(any(class(i) %in% c("gam"))) r <- rsquared.gam(i) else
                        
                        r <- list(family = "gaussian", link = "identity", method = "none", R.squared = NA)
                      
                      ret <- do.call(data.frame, r)
                      
                      ret <- data.frame(Response = all_vars_merMod(formula(i))[1], ret)
                      
                      # if(ncol(ret) != 5) ret[, ncol(ret) + 1] <- NA
                      
                      return(ret)
                      
  } )
  
  if(length(unique(sapply(ret, ncol))) != 1) {
    
    nc <- max(sapply(ret, ncol))
    
    ret <- lapply(ret, function(i)
      
      if(ncol(i) == nc) i else {
        
        data.frame(
          Response = i$Response,
          family = i$family,
          link = i$link,
          method = i$method,
          Marginal = i$R.squared,
          Conditional = NA
        )
        
      }
      
    )
    
  }
  
  ret <- do.call(rbind, ret)
  
  # ret[, which(sapply(ret, is.numeric))] <- round(ret[, which(sapply(ret, is.numeric))], 2)
  
  rownames(ret) <- NULL
  
  return(ret)
  
}

#' R^2 for lm objects
#' 
#' @keywords internal
#' 
rsquared.lm <- function(model)
  
  list(family = "gaussian", link = "identity", method = "none", R.squared = summary(model)$r.squared)

#' R^2 for gls objects
#' 
#' @keywords internal
#' 
rsquared.gls <- function(model) {
  
  X <- model.matrix(eval(model$call$model), GetData(model))
  
  sigmaF <- var(as.vector(model$coefficients %*% t(X)))
  
  sigmaE <- var(resid(model))
  
  list(family = "gaussian", link = "identity", method = "none", R.squared = sigmaF / (sigmaF + sigmaE))
  
}

#' R^2 for glm objects
#' 
#' @keywords internal
#' 
rsquared.glm <- function(model, method = "nagelkerke") {
  
  if(is.null(method)) method <- "nagelkerke"
  
  link <- family(model)$link
  
  family. <- family(model)$family
  
  family. <- gsub("(.*)\\(.*\\)", "\\1", family.)
  
  if(family. == "gaussian") {
    
    method <- "none"
    
    r <- summary(lm(model$formula, model$data))$r.squared
    
  }
  
  if(method == "mcfadden") {
    
    r <- 1 - (model$deviance / model$null.deviance)
    
  }
  
  if(method %in% c("coxsnell", "nagelkerke")) {
    
    n <- nobs(model)
    
    r <- 1 - exp((model$deviance - model$null.deviance) / n)
    
    if(method == "nagelkerke")
      
      r <- r / (1 - exp(-1 * model$null.deviance / n))
    
  }
  
  
  list(family = family., link = link, method = method, R.squared = r)
  
}

#' R^2 for phylolm objects
#' 
#' @keywords internal
#' 
# rsquared.phylolm <- function(model) {
# 
#   family. <- ifelse(class(model) == "phylolm", "gaussian", model$method)
# 
#   link <- ifelse(class(model) == "phylolm", "identity", "?")
# 
#   list(family = family., link = NA, method = "none", R.squared = NA)
# 
# }

#' R^2 for merMod objects
#' 
#' @keywords internal
#' 
rsquared.merMod <- function(model) {
  
  X <- model.matrix(model)
  
  sigmaF <- var(as.vector(lme4::fixef(model) %*% t(X)))
  
  sigma <- unclass(lme4::VarCorr(model))
  
  sigmaL <- sum(sapply(1:length(sigma), function(i) {
    
    sigma. <- sigma[[i]]
    
    if(all(rownames(sigma.) %in% colnames(X))) X. <- X else
      
      X. <- do.call(cbind, model.matrix(model, type = "randomListRaw")) 
    
    Z <- as.matrix(X.[, rownames(sigma.), drop = FALSE])
    
    sum(rowSums((Z %*% sigma.) * Z))/nrow(X.)
    
  } ) )
  
  sigmaE <- attr(lme4::VarCorr(model), "sc")^2
  
  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)
  
  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)
  
  list(family = "gaussian", link = "identity", method = "none", Marginal = mar, Conditional = con)
  
}

#' R^2 for lme objects
#' 
#' @keywords internal
#' 
rsquared.lme <- function(model) {
  
  X <- model.matrix(eval(model$call$fixed)[-2], droplevels(model$data))
  
  sigmaF <- var(as.vector(nlme::fixef(model) %*% t(X)))
  
  sigma <- GetVarCov(model)
  
  sigmaL <- sum(sapply(sigma, function(i) {
    
    if(all(rownames(i) %in% colnames(X))) X. <- X else
      
      X. <- model.matrix(model$modelStruct$reStruct,
                         data = model$data[rownames(model$fitted), , drop = FALSE]) 
    
    Z <- as.matrix(X.[, rownames(i), drop = FALSE])
    
    sum(rowSums((Z %*% i) * Z))/nrow(X.)
    
  } ) )
  
  sigmaE <- summary(model)$sigma^2
  
  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)
  
  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)
  
  list(family = "gaussian", link = "identity", method = "none", Marginal = mar, Conditional = con)
  
}

#' R^2 for glmer objects
#' 
#' @keywords internal
#' 
rsquared.glmerMod <- function(model, method = "trigamma") {
  
  if(is.null(method)) method <- "trigamma"
  
  link <- family(model)$link
  
  family. <- family(model)$family
  
  family. <- gsub("(.*)\\(.*\\)", "\\1", family.)
  
  if(family. == "Negative Binomial") rsquared.negbin(model, method) else {
    
    X <- model.matrix(model)
    
    sigmaF <- var(as.vector(lme4::fixef(model) %*% t(X)))
    
    sigma <- unclass(lme4::VarCorr(model))
    
    data <- GetData(model)
    
    if(family. == "poisson") {
      
      if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")
      
      sigmaL <- sum(GetOLRE(sigma, model, X, data, RE = "RE"))
      
      sigmaE <- sum(GetOLRE(sigma, model, X, data, RE = "OLRE"))
      
      rand <- onlyBars(formula(model))
      
      f <- paste(all_vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = TRUE))
      
      nullmodel <- suppressWarnings(lme4::glmer(formula(f), family = poisson(link = link), data = data))
      
      # lambda <- attr(VarCorr(model), "sc")^2
      
      lambda <- exp(fixef(nullmodel)[1] + (sigmaL + sigmaE)/2)
      
      omega <- 1
      
      if(link == "mu^0.5") sigmaD <- 0.25 * omega else {
        
        if(link == "log") {
          
          nu <- omega / lambda
          
          if(method == "delta") sigmaD <- nu
          
          if(method == "lognormal") sigmaD <- log(1 + nu)
          
          if(method == "trigamma") sigmaD <- trigamma(1/nu)
          
        } else if(link == "sqrt") {
          
          method <- "delta"
          
          if(method == "delta") sigmaD <- 0.25*omega else stop("Unsupported method!")
          
        }
        
        else stop("Unsupported link function!")
        
      }
      
    } else
      
      if(family. == "binomial") {
        
        if(method == "trigamma") method <- "theoretical"
        
        if(!method %in% c("theoretical", "delta")) stop("Unsupported method!")
        
        sigmaL <- sum(GetOLRE(sigma, model, X, data, RE = "all"))
        
        sigmaE <- 0
        
        if(method == "theoretical") sigmaD <- pi^2/3
        
        if(method == "delta") {
          
          rand <- onlyBars(formula(model))
          
          f <- paste(all_vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))
          
          nullmodel <- suppressWarnings(lme4::glmer(formula(f), family = binomial(link = link), data = data))
          
          vt <- sum(unlist(lme4::VarCorr(nullmodel)))
          
          pmean <- plogis(as.numeric(fixef(nullmodel)) - 0.5 * vt *
                            tanh(as.numeric(fixef(nullmodel)) * (1 + 2 * exp(-0.5 * vt))/6))
          
          sigmaD <- 1/(pmean * (1 - pmean))
          
        }
        
      } else
        
        if(family. == "Gamma") {
          
          if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")
          
          sigmaL <- sum(GetOLRE(sigma, model, X, data, RE = "all"))
          
          sigmaE <- 0
          
          lambda <- attr(lme4::VarCorr(model), "sc")^2
          
          omega <- 1
          
          if(link == "log") {
            
            nu <- omega / lambda
            
            if(method == "delta") sigmaD <- 1/nu
            
            if(method == "lognormal") sigmaD <- log(1 + 1/nu)
            
            if(method == "trigamma") sigmaD <- trigamma(nu)
            
          } else
            
            if(link == "inverse") {
              
              if(method != "delta") method <- "delta"
              
              nu <- omega / lambda
              
              sigmaD <- 1/(nu * lambda^2)
              
            } else stop("Unsupported link function!")
          
        } else stop("Unsupported family!")
    
    mar <- (sigmaF) / (sigmaF + sigmaL + sigmaD + sigmaE)
    
    con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaD + sigmaE)
    
    list(family = family., link = link, method = method, Marginal = mar, Conditional = con)
    
  }
  
}

#' R^2 for negbin objects
#' 
#' @keywords internal
#' 
rsquared.negbin <- function(model, method = "trigamma") {
  
  if(is.null(method)) method <- "trigamma"
  
  link <- family(model)$link
  
  family. <- family(model)$family
  
  family. <- gsub("(.*)\\(.*\\)", "\\1", family.)
  
  X <- model.matrix(model)
  
  sigmaF <- var(as.vector(fixef(model) %*% t(X)))
  
  sigma <- unclass(lme4::VarCorr(model))
  
  rand <- sapply(lme4::findbars(formula(model)), function(x) as.character(x)[3])
  
  rand <- rand[!duplicated(rand)]
  
  data <- GetData(model)
  
  idx <- sapply(strsplit(rand, "\\:"), function(x) {
    
    length(unique(data[, x])) == nrow(data)
    
  } )
  
  sigma.names <- strsplit(names(sigma), "\\.")
  
  idx. <- sapply(sigma.names, function(x) !any(x %in% rand[idx]))
  
  sigmaL <- sum(sapply(sigma[idx.], function(i) {
    
    Z <- as.matrix(X[, rownames(i), drop = FALSE])
    
    sum(rowSums((Z %*% i) * Z))/nrow(X)
    
  } ) )
  
  if(family. == "Negative Binomial") {
    
    rand <- onlyBars(formula(model))
    
    f <- paste(all_vars_trans(formula(model))[1], " ~ 1 + ", onlyBars(formula(model), slopes = FALSE))
    
    nullmodel <- suppressWarnings(lme4::glmer.nb(formula(f), family = negative.binomial(link = link), data))
    
    # nullmodel <- update(model, formula(paste(". ~ 1 +", onlyBars(formula(model)))))
    
    lambda <- as.numeric(exp(fixef(nullmodel) + 0.5 * sum(unlist(lme4::VarCorr(nullmodel)))))
    
    theta <- lme4::getME(model, "glmer.nb.theta")
    
    if(link == "log") {
      
      if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")
      
      nu <- (1/lambda) + (1/theta)
      
      if(method == "delta") sigmaE <- nu
      
      if(method == "lognormal") sigmaE <- log(1 + nu)
      
      if(method == "trigamma") sigmaE <- trigamma(nu^(-1))
      
    } else if(link == "sqrt") {
      
      method <- "delta"
      
      if(method == "delta") sigmaE <- 0.25*(1 + (lambda / theta)) else stop("Unsupported method!")
      
    } else stop("Unsupported link function!")
    
  }
  
  mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)
  
  con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)
  
  list(family = family., link = link, method = method, Marginal = mar, Conditional = con)
  
}

# R^2 for glmmTMB
#
# @keywords internal
#
rsquared.glmmTMB <- function(model) {
  
  list(family = family(model)$family, 
       link = family(model)$link, 
       method = "none", 
       Marginal = performance::r2(model)$R2_marginal, 
       Conditional = performance::r2(model)$R2_conditional)
  
}

#' R^2 for glmmPQL objects
#' 
#' @keywords internal
#' 
rsquared.glmmPQL <- function(model, method = "trigamma") {
  
  if(is.null(method)) method <- "trigamma"
  
  link <- model$family$link
  
  family. <- model$family$family
  
  X <- model.matrix(eval(model$call$fixed)[-2], droplevels(model$data))
  
  sigmaF <- var(as.vector(nlme::fixef(model) %*% t(X)))
  
  sigma <- nlme::VarCorr(model)
  
  sigmaL <- sum(as.numeric(sigma[!grepl("=|Residual", rownames(sigma)), 1]))
  
  data <- GetData(model)
  
  if(family. %in% c("gaussian")) l <- rsquared.lme(model) else {
    
    if(family. %in% c("poisson", "quasipoisson")) {
      
      f <- paste(all_vars_trans(formula(model))[1], " ~ 1")
      
      if(family. == "poisson")
        
        nullmodel <- suppressWarnings(MASS::glmmPQL(formula(f), random = model$call$random, family = poisson(link = link), data = data, verbose = FALSE)) else
          
          nullmodel <- suppressWarnings(MASS::glmmPQL(formula(f), random = model$call$random, family = quasipoisson(link = link), data = data, verbose = FALSE))
        
        lambda <- as.numeric(exp(fixef(nullmodel) + 0.5 * sum(unlist(GetVarCov(nullmodel)))))
        
        if(family. == "poisson") omega <- 1 else omega <- as.numeric(sigma[nrow(sigma), 1])
        
        if(link == "mu^0.5") sigmaE <- 0.25 * omega else {
          
          if(link == "log") {
            
            if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")
            
            nu <- omega / lambda
            
            if(method == "delta") sigmaE <- nu
            
            if(method == "lognormal") sigmaE <- log(1 + (nu))
            
            if(method == "trigamma") sigmaE <- trigamma(1/nu)
            
          } else stop("Unsupported link function!")
          
        }
        
    } else if(family. %in% c("binomial", "quasibinomial")) {
      
      if(method == "trigamma") method <- "delta"
      
      if(!method %in% c("theoretical", "delta")) stop("Unsupported method!")
      
      if(method == "theoretical") sigmaE <- sigmaD <- pi^2/3
      
      if(method == "delta") {
        
        f <- paste(all_vars_trans(formula(model))[1], " ~ 1")
        
        if(family. == "binomial")
          
          nullmodel <- suppressWarnings(MASS::glmmPQL(formula(f), random = model$call$random, family = binomial(link = link), data = data, verbose = FALSE)) else
            
            nullmodel <- suppressWarnings(MASS::glmmPQL(formula(f), random = model$call$random, family = quasibinomial(link = link), data = data, verbose = FALSE))
          
          vt <- sum(unlist(GetVarCov(nullmodel)))
          
          pmean <- plogis(as.numeric(fixef(nullmodel)) - 0.5 * vt *
                            tanh(as.numeric(fixef(nullmodel)) * (1 + 2 * exp(-0.5 * vt))/6))
          
          sigmaE <- 1/(pmean * (1 - pmean))
          
      }
      
    } else if(family. %in% c("Gamma")) {
      
      if(link == "log") {
        
        if(!method %in% c("delta", "lognormal", "trigamma")) stop("Unsupported method!")
        
        nu <- 1 / as.numeric(nlme::VarCorr(model)[nrow(nlme::VarCorr(model)), 1])
        
        if(method == "delta") sigmaE <- 1 / nu
        
        if(method == "lognormal") sigmaE <- log(1 + 1/nu)
        
        if(method == "trigamma") sigmaE <- trigamma(nu)
        
      } else stop("Unsupported link function!")
      
    } else stop("Unsupported family!")
    
    mar <- (sigmaF) / (sigmaF + sigmaL + sigmaE)
    
    con <- (sigmaF + sigmaL) / (sigmaF + sigmaL + sigmaE)
    
    l <- list(family = family., link = link, method = method, Marginal = mar, Conditional = con)
    
  }
  
  return(l)
  
}

#' R^2 for Sarlm objects
#' 
#' @keywords internal
#'
rsquared.Sarlm <- function(model) {
  
  list(family = "gaussian", link = "identity", method = "none", 
       R.squared = summary(model, Nagelkerke = TRUE)$NK)
  
}

#' R^2 for gam objects
#' 
#' @keywords internal
#'
rsquared.gam <- function(model) {
  
  link <- family(model)$link
  
  family. <- family(model)$family
  
  if(family. == "gaussian") r <- summary(model)$r.sq else r <- summary(model)$dev.expl
  
  list(family = family., link = link, method = "none", R.squared = r)
  
}
