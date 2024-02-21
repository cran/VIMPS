#' @importFrom stats as.formula predict reorder
#' @importFrom utils read.csv write.csv write.table

#*******************************************************************************
#
# The functions to drive execution
#
#*******************************************************************************

#*******************************************************************************
# Default ML models used in simulation
#*******************************************************************************

# Generates predictions for samples using a full random forest tree. This is
# utilized in prediction thresholds and evaluating performance on the full
# set of variables.
normal_model = function(df_t, df_T, formula, mtry=NULL, min.node.size=NULL) {
  # Train model
  model = ranger::ranger(as.formula(formula), data=df_t, probability=TRUE, mtry=mtry, min.node.size=min.node.size)
  # Find out which column has probability predictions for class 1
  class_1_col = match("1", colnames(model$predictions))
  # Get the class 1 probability predictions for the testing data
  yh = predict(model, data=df_T)$predictions[,class_1_col]
  return(yh)
}

# Generates predictions for samples using just a tree one level deep. This is
# utilized in evaluating performance when removing variables/domains or
# replacing them with knockoffs
one_tree_model = function(df_t, df_T, formula, mtry=NULL, min.node.size=NULL) {
  # Train model
  model = ranger::ranger(as.formula(formula), data=df_t, probability=TRUE, num.trees=1, mtry=mtry, min.node.size=min.node.size)
  # Find out which column has probability predictions for class 1
  class_1_col = match("1", colnames(model$predictions))
  # Get the class 1 probability predictions for the testing data
  yh = predict(model, data=df_T)$predictions[,class_1_col]
  return(yh)
}


#*******************************************************************************
# Use Yoden's J score to find the threshold to use
#*******************************************************************************

# We need to wrap the predictions in a try-catch block so we needed to create
# an individual function for it because of how R syntax works
calc_roc = function(yh, y) {
  tryCatch(
    {
      pred = ROCR::prediction(yh, y)
      return(pred)
    },
    error=function(e) {
      return(NULL)
    }
  )
}

# Use x-fold validation with yoden's j score to determine the threshold for
# classifying a sample as positive or negative
calc_threshold = function(dat, dep_var, num_folds=10, model=normal_model, mtry=NULL, min.node.size=NULL) {
  folds = caret::createFolds(dat[,dep_var], k=num_folds)
  scores = c()

  for(k in 1:num_folds) {
    # Get our train and test splits for a given fold
    df_t = dat[-folds[[k]], ]
    df_T = dat[folds[[k]], ]

    # Train model and get results
    formula = paste(dep_var, " ~ .")
    yh = model(df_t, df_T, formula, mtry, min.node.size)

    # Get our tpr, fpr, and threshold values
    pred = calc_roc(yh, df_T[,dep_var])

    # Catch case where test split has no positive outcomes
    if(is.null(pred)) {
      next
    }

    perf = ROCR::performance(pred, "tpr", "fpr")
    tpr = perf@y.values[[1]]
    fpr = perf@x.values[[1]]
    thresholds = perf@alpha.values[[1]]

    # Calc yoden's J
    j_scores = tpr-fpr
    max_idx = which(j_scores==max(j_scores))
    # Take the min because sometimes we can have infinity as one of the max vals
    yoden_j = min(thresholds[max_idx])
    scores = append(scores, yoden_j)
  }

  message("Thresholds for each k fold: ")
  message(scores)

  # Get the average score and return it
  threshold = list(avg = mean(scores), all = scores)
  message("Average threshold: ")
  message(threshold$avg)

  return (threshold)
}


#*******************************************************************************
# Calc the model metrics using all the variables
#*******************************************************************************

calc_var_metrics = function(dat, dep_var, threshold, iterations=500, model=normal_model, mtry=NULL, min.node.size=NULL) {
  preds = matrix(nrow=0, ncol=nrow(dat))

  for(i in 1:iterations) {
    # Do bootstrap resampling
    dat_bs = dat[sample(nrow(dat), replace=TRUE),]

    # Create stratified train/test splits w/o VOI
    idx_tr = caret::createDataPartition(dat_bs[,dep_var], p=0.8, list=FALSE)
    df_t = dat_bs[idx_tr,]
    df_T= dat_bs[-idx_tr,]

    # Train model and get results
    formula = paste(dep_var, " ~ .")
    yh = model(df_t, df_T, formula, mtry, min.node.size)

    # Create vector to store predictions
    vec = rep(NA, nrow(dat))
    rows = as.integer(rownames(df_T))

    # Collect the results
    for(i in 1:length(rows)) {
      pred = yh[i]
      row_idx = rows[i]
      vec[row_idx] = pred
    }

    # Add the results to the matrix
    preds = rbind(preds, vec)
  }

  # Get the average prediction for each variable
  mus = colMeans(preds, na.rm=TRUE)
  # Label each sample positive/negative on cut off
  yh = mus > threshold

  # Computer our confusion matrix
  mtx = table(dat[, dep_var], yh)
  tn = mtx[1,1]
  fp = mtx[1,2]
  fn = mtx[2,1]
  tp = mtx[2,2]

  # Compute the metrics and save them
  acc = (tp+tn)/length(dat[, dep_var])
  sens = tp/(tp+fn)
  spec = tn/(tn+fp)

  message(paste("Accuracy:", acc))
  message(paste("Sensitivity:", sens))
  message(paste("Specificity:", spec))

  return (list(acc=acc, sens=sens, spec=spec))
}


#*******************************************************************************
# Calc knock-offs
#*******************************************************************************

# Wrapper script for client-side operation
gen_kos = function(dat, dep_var, iterations=100, ko_path=NULL) {
  # Set-up the folder to store the KO values
  if(is.null(ko_path)) {
    ko_path = paste(tempdir(), "kos/", sep="/")
  }
  if(dir.exists(ko_path)) {
    unlink(ko_path, recursive=TRUE)
  }
  dir.create(ko_path)

  message(paste("Creating", iterations, "kos in", ko_path))

  for(i in 1:iterations) {
    gen_ko(dat, dep_var, i, ko_path)
  }
}

# regular script for server-side operation
gen_ko = function(dat, dep_var, iteration, ko_path) {
  # Set-up data in format needed for knock-off package
  X = dat[,!(names(dat) %in% c(dep_var))]
  X = data.matrix(X)

  # Create knock-off variables
  Xk = knockoff::create.second_order(X)

  # Store results to disk
  fout = paste(ko_path, "ko_set_", iteration, ".csv", sep="")
  write.table(Xk, fout, row.names=FALSE, col.names=TRUE, sep=",")
}


#*******************************************************************************
# Calc KO Probs
#*******************************************************************************

# Wrapper function for client-side operation
calc_ko_probs = function(dat, dep_var, doms, ko_path=NULL, results_path=NULL, iterations=500, model=one_tree_model, mtry=NULL, min.node.size=NULL) {
  message(paste("Calculating ", iterations, " rounds of probabilities using knockoffs"))

  # Set default ko_path if none provided
  if(is.null(ko_path)) {
    ko_path = paste(tempdir(), "kos/", sep="/")
  }

  num_kos=length(list.files(ko_path))

  # Set-up the folder to store the results
  if(is.null(results_path)) {
    results_path = paste(tempdir(), "results/", sep="/")
  }
  if(dir.exists(results_path)) {
    unlink(results_path, recursive=TRUE)
  }
  dir.create(results_path)

  num_vois = length(unique(doms$domain))

  for(i in 1:num_vois) {
    calc_ko_prob(dat, doms, ko_path, results_path, i, dep_var, iterations, num_kos, model)
  }
}

calc_ko_prob = function(dat, doms, ko_path, results_path, voi_idx, dep_var, iterations, num_kos, model, mtry=NULL, min.node.size=NULL) {
  # Only keep the variables we have in our dictionary
  vars = c(doms$variable, dep_var)
  dat = dat[,names(dat) %in% vars]

  # Get our domains and VOIs
  dom = unique(doms$domain)[voi_idx]
  vois = doms[doms$domain==dom,]$var

  # Set up our storage for the predicted probabilities
  ko_pred = matrix(nrow=0, ncol=nrow(dat))

  for(i in 1:iterations) {
    # Do bootstrap resampling
    dat_bs = dat[sample(nrow(dat), replace=TRUE),]

    # Create knock-off variables and assign knock-off over top regular VOI
    ko_file = paste(ko_path, "/ko_set_", sample.int(num_kos, 1, replace=TRUE), ".csv", sep="")
    Xk = read.csv(ko_file)
    dat_bs[,vois] = Xk[,vois]

    # Create stratified train/test splits w/ knock-off var
    idx_tr = caret::createDataPartition(dat_bs[,dep_var], p=0.8, list=FALSE)
    df_t = dat_bs[idx_tr,]
    df_T= dat_bs[-idx_tr,]

    # Train model and get results
    formula = paste(dep_var, " ~ .")
    yh = model(df_t, df_T, formula, mtry, min.node.size)

    # Create vector to store predictions
    vec = rep(NA, nrow(dat))
    rows = as.integer(rownames(df_T))

    for(i in 1:length(rows)) {
      pred = yh[i]
      row_idx = rows[i]
      vec[row_idx] = pred
    }

    # Add the results to the matrix
    ko_pred = rbind(ko_pred, vec)
  }

  fout = paste(results_path, dom, "_ko.csv", sep="")
  write.table(ko_pred, fout, row.names=FALSE, col.names=FALSE, sep=",")
}


#*******************************************************************************
# Calc dom Probs
#*******************************************************************************

# Client-side wrapper function
calc_dom_probs = function(dat, dep_var, doms, results_path=NULL, iterations=500, model=one_tree_model, mtry=NULL, min.node.size=NULL) {
  message(paste("Calculating", iterations, "rounds of probabilities by domain removal"))

  # Set-up the folder to store the results
  if(is.null(results_path)) {
    results_path = paste(tempdir(), "results/", sep="/")
  }
  if(dir.exists(results_path)) {
    unlink(results_path, recursive=TRUE)
  }
  dir.create(results_path)

  num_vois = length(unique(doms$domain))

  for(i in 1:num_vois) {
    calc_dom_prob(dat, doms, results_path, i, dep_var, iterations, model)
  }
}

calc_dom_prob = function(dat, doms, results_path, voi_idx, dep_var, iterations, model, mtry=NULL, min.node.size=NULL) {
  # Only keep the variables we have in our dictionary
  vars = c(doms$variable, dep_var)
  dat = dat[,names(dat) %in% vars]

  # Get our domains and VOIs
  dom = unique(doms$domain)[voi_idx]
  vois = doms[doms$domain==dom,]$var

  # Set up our storage for the predicted probabilities
  dom_pred = matrix(nrow=0, ncol=nrow(dat))

  for(i in 1:iterations) {
    # Do bootstrap resampling
    dat_bs = dat[sample(nrow(dat), replace=TRUE),]

    # Remove the variable of interest
    dat_bs = dat_bs[, !(names(dat_bs) %in% vois)]

    # Create stratified train/test splits w/ knock-off var
    idx_tr = caret::createDataPartition(dat_bs[,dep_var], p=0.8, list=FALSE)
    df_t = dat_bs[idx_tr,]
    df_T= dat_bs[-idx_tr,]

    # Train model and get results
    formula = paste(dep_var, " ~ .")
    yh = model(df_t, df_T, formula, mtry, min.node.size)

    # Create vector to store predictions
    vec = rep(NA, nrow(dat))
    rows = as.integer(rownames(df_T))

    for(i in 1:length(rows)) {
      pred = yh[i]
      row_idx = rows[i]
      vec[row_idx] = pred
    }

    # Add the results to the matrix
    dom_pred = rbind(dom_pred, vec)
  }

  fout = paste(results_path, dom, "_dom.csv", sep="")
  write.table(dom_pred, fout, row.names=FALSE, col.names=FALSE, sep=",")
}


#*******************************************************************************
# Calc Metrics
#*******************************************************************************

calc_metrics = function(dat, dep_var, cutoff, test_type, metrics, input_path=NULL, output_file=NULL) {
  message("Calculating variable importance")

  if(is.null(input_path)) {
    input_path = paste(tempdir(), "results/", sep="/")
  }

  if(is.null(output_file)) {
    output_file = paste(tempdir(), "/", test_type, "_output.csv", sep="")
  }

  y = dat[,dep_var]

  # Get our list of result files
  pattern = paste("*_", test_type, ".csv", sep="")
  files = list.files(input_path, pattern=pattern)

  results = matrix(nrow=0, ncol=6)
  var_names = c()

  # Iterate over our result files
  for (i in 1:length(files)) {
    # Open our results
    file_name = paste(input_path, files[i], sep="")
    df = read.csv(file_name, stringsAsFactors=FALSE, header=FALSE)

    # Compute the predictions of the results via the cutoff
    mus = colMeans(df, na.rm=TRUE)
    yh = mus > cutoff

    # Computer our confusion matrix
    mtx = table(y, yh)
    tn = mtx[1,1]
    fp = mtx[1,2]
    fn = mtx[2,1]
    tp = mtx[2,2]

    # Compute the metrics and save them
    acc = (tp+tn)/length(y)
    sens = tp/(tp+fn)
    spec = tn/(tn+fp)

    var_names = append(var_names, gsub(substr(pattern, 2, nchar(pattern)), "", files[i]))
    results = rbind(results, c(acc, sens, spec, metrics$acc-acc, metrics$sens-sens, metrics$spec-spec))
  }

  # Convert results to Dataframe and save to disk
  results = data.frame(results)
  colnames(results) = c("acc", "sens", "spec", "delta_acc", "delta_sens", "delta_spec")
  results$variable = var_names
  write.csv(results, output_file, row.names=FALSE)

  return(results)
}


#*******************************************************************************
# Run entire process
#*******************************************************************************

#' @title calc_vimps
#'
#' @description Calculate the variable importance of the domains for a given
#'    dataset
#'
#' @param dat A dataframe of data
#' @param dep_var The dependent variable in the dat
#' @param doms A dataframe of the variables in dat and the domain they belong to
#' @param calc_ko True/False to calculate the knock_off importance
#' @param calc_dom True/False to calculate the domain importance
#' @param num_folds The number of folds to use while calculating the classification threshold for predictions
#' @param num_kos The number of sets of knock off variables to create
#' @param model_all The model to use in full ensemble mode in calculations
#' @param model_subset The model to use sigularly for building ensembles from
#' @param mtry The mtry value to use in the random forests
#' @param min.node.size The min.node.size value to use in the random forests
#' @param iterations Number of trees to build while calculating variable importance
#' @param ko_path Where to store the knock off variable sets
#' @param results_path Where to store the intermediary results for calculating variable importance
#' @param output_file_ko Where to store the results of the knock off variable importance
#' @param output_file_dom Where to store the results of the domain variable importance
#'
#' @return List with 1) Threshold for binary class labeling 2) Model metrics using all variables 3) Model metrics using knock-off variables 4) Variable importance with knock-offs
#'
#' @examples
#' calc_vimps(
#'   data.frame(
#'     X1=c(2,8,3,9,1,4,3,8,0,9,2,8,3,9,1,4,3,8,0,9),
#'     X2=c(7,2,5,0,9,1,8,8,3,9,7,2,5,0,9,1,8,8,3,9),
#'     Y=c(0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1)),
#'  "Y",
#'  data.frame(domain=c('X1','X2'),
#'  variable=c('X1','X2')),
#'  num_folds=2,
#'  num_kos=1,
#'  iterations=50)
#'
#' @export
calc_vimps = function(dat, dep_var, doms, calc_ko=TRUE, calc_dom=FALSE,
                      num_folds=10, num_kos=100, model_all=normal_model,
                      model_subset=one_tree_model, mtry=NULL, min.node.size=NULL,
                      iterations=500, ko_path=NULL, results_path=NULL,
                      output_file_ko=NULL, output_file_dom=NULL) {
  results = list()

  message("")
  message("#*******************************************************************************")
  message("# Calculating thresholds")
  message("#*******************************************************************************")
  threshold = calc_threshold(dat, dep_var, num_folds, model_all, mtry, min.node.size)
  results[['thresholds']] = threshold

  message("")
  message("#*******************************************************************************")
  message("# Calculating variable metrics")
  message("#*******************************************************************************")
  metrics = calc_var_metrics(dat, dep_var, threshold$avg, iterations, model_all, mtry, min.node.size)
  results[['metrics']] = metrics

  message("")
  message("#*******************************************************************************")
  message("# Generating knockoff variables")
  message("#*******************************************************************************")
  gen_kos(dat, dep_var, num_kos, ko_path)

  if(calc_ko) {
    message("")
    message("#*******************************************************************************")
    message("# Calculating knockoff probabilities")
    message("#*******************************************************************************")
    calc_ko_probs(dat, dep_var, doms, ko_path, results_path, iterations, model_subset, mtry, min.node.size)

    message("")
    message("#*******************************************************************************")
    message("# Calculating knockoff metrics")
    message("#*******************************************************************************")
    results_ko = calc_metrics(dat, dep_var, threshold$avg, "ko", metrics, results_path, output_file_ko)
    results[['results_ko']] = results_ko
    results_ko = results_ko[c("variable", "delta_acc")]
    colnames(results_ko) = c("variable", "importance")
    results[['ko_importance']] = results_ko
  }

  if(calc_dom) {
    message("")
    message("#*******************************************************************************")
    message("# Calculating domain probabilities")
    message("#*******************************************************************************")
    calc_dom_probs(dat, dep_var, doms, results_path, iterations, model_subset, mtry, min.node.size)

    message("")
    message("#*******************************************************************************")
    message("# Calculating domain metrics")
    message("#*******************************************************************************")
    results_dom = calc_metrics(dat, dep_var, threshold$avg, "dom", metrics, results_path, output_file_dom)
    results[['results_dom']] = results_dom
    results_dom = results_dom[c("variable", "delta_acc")]
    colnames(results_dom) = c("variable", "importance")
    results[['dom_importance']] = results_dom
  }

  return(results)
}


#*******************************************************************************
# Visualize Results
#*******************************************************************************

#' @title graph_results
#'
#' @description Graph the variable importance results from calc_vimps
#'
#' @param results The results from calc_vimps
#' @param object Which object from results to use for graphing results
#'
#' @return No return value
#'
#' @importFrom ggplot2 ggplot aes geom_point
#'
#' @export
graph_results = function(results, object) {
  df = results[[object]]

  ggplot2::ggplot(data=df, ggplot2::aes_string(x="importance", y=reorder("variable", "importance"))) +
    ggplot2::geom_point() +
    ggplot2::theme_bw() +
    ggplot2::xlab("KO Domain Variable Importance") +
    ggplot2::ylab("")
}

