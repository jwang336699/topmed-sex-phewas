#' Performs a generic PheWAS (logistic regression coefficient test, Bonferroni correction) on 
#' data extracted from PIC-SURE.
#'
#' @param dt A \code{data.frame} produced by running a PIC-SURE query
#' @param vdict An \code{data.frame} as produced by the "get_multiIndex_variablesDict"
#' function found in utils.R
#' @param varset A logical vector containing an entry for each row of vdict such that 
#' an entry is true iff it is a desired variable for analysis
#' @param dependent A character singleton containing the (simplified) name of the dependent variable
#' @param dep_cat A boolean singleton which is true iff dependent is categorical
#' @param dep_keep A rule taking values of dependent to a boolean; used to remove unwanted entries
#' @param dep_classes A rule taking values of dependent to a boolean; used for dichotomizing data
#'
#' @return An enriched \code{data.frame} which contains p-values, functions of p-values, and other
#' useful information about each variable appended to vdict
#' @import questionr R_lib/utils.R
#' @importFrom questionr odds.ratio
#' @importFrom R_lib/utils.R parsing_varNames
#' @export
phewas <- function(dt, vdict, varset, dependent, dep_cat, dep_keep, dep_classes) {
  # reshape variable names to eliminate spaces
  vdict[["df_name"]] <- parsing_varNames(vdict[["name"]])
  
  # separate variables into categorical and continuous
  cat_vars <- vdict[, "categorical"] == TRUE
  categorical <- vdict[cat_vars & varset, ][["df_name"]]
  continuous <- vdict[!cat_vars & varset, ][["df_name"]]
  
  # eliminate the dependent variable from the lists of independent variables
  dependent <- vdict[vdict[["simplified_name"]] == dependent,][["df_name"]]
  if (dep_cat) {
    categorical <<- categorical[-which(categorical == dependent)]
  } else {
    continuous <<- continuous[-which(continuous == dependent)]
  }
  
  # filter and cast dependent variable according to inputted rules
  dt <- dt[dep_keep(dt[[dependent]]),]
  dt$classified <- ifelse(dep_classes(dt[[dependent]]), 1, 0)
  dependent <- 'classified'
  
  # define statistical test function
  # test <- function(set, ind, dep, isCat) {
  #   ret <- list(NA)
  #   length(ret) <- 3
  #   if (isCat) {
  #     tbl <- table(set[[ind]], set[[dep]])
  #     ret[[1]] <- fisher.test(tbl)[[1]]
  #     ors <- lapply(2:nrow(tbl), function(x) { odds.ratio(rbind(tbl[1,], tbl[x,]))[,1:3] })
  #     ors <- Reduce(function(a,b) { rbind(a,b) }, ors, c())
  #     scaled_ors <- vapply(ors[,1], function(x) { if (x > 1) {x} else {1/x} }, numeric(1))
  #     chosen <- which.max(scaled_ors)
  #     ret[[2]] <- ors[chosen,1]
  #     ret[[3]] <- paste(ors[chosen,2], ors[chosen,3], sep = '-')
  #   } else {
  #     select <- set[[dep]] == '1'
  #     c1 <- set[select,]
  #     c2 <- set[!select,]
  #     results <- wilcox.test(c1[[ind]], c2[[ind]], conf.int = T)
  #     ret[[1]] <- results[[2]]
  #     ret[[2]] <- round(results[[9]], 5)
  #     ret[[3]] <- paste(vapply(results[[8]], function(x) { round(x, 5) }, numeric(1)), collapse = '-')
  #   }
  # }
  
  # define statistical test
  anova_model <- function(data, dependent_var, independent_var) {
    ret <- list(NA)
    length(ret) <- 4
    model <- glm(as.formula(paste(dependent_var, "~ 1 +", independent_var)),
                 data = data,
                 family = binomial(link="logit"))
    # model_reduced <- glm(as.formula(paste(dependent_var, "~ 1")),
    #                      data = data,
    #                      family = binomial(link="logit"))
    # ret[[1]] <- anova(model, model_reduced, test =  "LRT")[2, "Pr(>Chi)"]
    ret[[1]] <- summary(model)$coefficients[2,4]
    ret[[2]] <- exp(coef(model))[2]
    ret[[3]] <- suppressMessages(paste(round(confint(model)[2,1], 5), round(confint(model)[2,2], 5), sep = '-'))
    ret[[4]] <- nrow(data)
    
    return(ret)    
  }
  
  results = list()
  errors =  list()
  warnings = list()
  
  # perform all tests
  for (independent in c(categorical_varnames, continuous_varnames)) {
    int <- na.omit(dt[, c(dependent, independent)])
    int <- int[int[[independent]] != '',]
    tryCatch({
      results[[independent]] <- anova_model(int, dependent, independent)
      #                error_list[[independent_varname]] <- NA
      #               warning_list[[independent_varname]] <- NA
    },
    error = function(e) {
      print(paste("error", independent))
      errors[[independent]] <<- e
    },
    warning = function(w) {
      print(paste("warning", independent))
      warnings[[independent]] <<- w                 
    }
    )
  }
  
  # bind into a dataframe with columns: name, p, odds ratio, odds ratio CI, n
  aggr <- as.data.frame(cbind(names(results), Reduce(function(a,b) { rbind(a,b) }, results, c())))
  colnames(aggr) <- c('df_name', 'pvalues', 'odds_ratio', 'or_CI', 'n')
  for (name in names(aggr)) { aggr[[name]] <- unlist(aggr[[name]]) }
  aggr[['log_pvalues']] <- -log10(aggr[['pvalues']])
  
  # integrate with original variable dictionary, adjust using Bonferroni
  vdict_enhanced <- dplyr::left_join(vdict, aggr, by="df_name")
  vdict_enhanced$adj_pvalues <- p.adjust(vdict_enhanced$pvalues, method="bonferroni")
  vdict_enhanced$log_adj_pvalues <- -log10(vdict_enhanced$adj_pvalues)
  
  return(vdict_enhanced)
}