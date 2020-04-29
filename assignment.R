library(tidyverse)
library(tuneR)
library(tsfeatures)
library(zoo)
library(glmnet)

###
# functions
###

# identifies intervals in a time series 
# where the standard deviation within that window is higher than a threshold
identify_intervals <- function(Y, windowSize=9000, thresholdSize=0.6) {
  Y_scaled = abs(scale(Y))
  testStat = rollapply(Y_scaled, windowSize, mean)
  h = which(testStat > thresholdSize)
  
  intervals <- data.frame(start=numeric(), end=numeric())
  start = 0
  for (i in 1:length(h)) {
    if (start == 0) {
      start = h[i]
    }
    
    if (i+1 > length(h) | h[i+1] != h[i] + 1) {
      intervals <- add_row(intervals, start=max(start-4000,1), end=min(h[i]+5000,length(Y)))
      start = 0
    }
  }
  
  return(intervals)
}

LRclassify = function(waveseq)
{
  maxPos = which.max(waveseq)
  minPos = which.min(waveseq)
  call = ifelse(maxPos < minPos, "Left", "Right")
  return(call)
}

###
# data loading
###

files <- list.files('input/', recursive=T)
wave_seqs <- list()
for (f in files) {
  wave_seqs <- append(wave_seqs, list(readWave(paste0('input/', f))))
}

# extract labels that follow the underscore in the file name
labels <- unlist(lapply(strsplit(files, '/'), function(x) x[[length(x)]]))
labels <- unlist(lapply(strsplit(labels, '_'), '[[', 1))
labels_ind <- lapply(labels, function(x) strsplit(x, '')[[1]])

###
# interval
###

wave_seq_events = list()
for(i in 1:length(wave_seqs)) {
  print(i)
  wave_file = wave_seqs[[i]]
  Y = wave_file@left
  intervals = identify_intervals(Y)
  intervals = intervals %>% filter(if (end != length(Y) || start != 1) end - start > 9200 else T)
  if (nrow(intervals) != length(labels_ind[[i]])) {
    intervals = intervals %>% filter(end - start > 10500)
  }
  waves = apply(intervals, 1, function(x) Y[x[1]:x[2]])
  wave_seq_events[[i]] = waves
}

###
# feature extraction
###

Y_list <- unlist(wave_seq_events, recursive=F)
resamp_Y_list <- lapply(Y_list, function(x) x[seq(1, length(x), by=20)]) # downsample sequences for performance improvement
y <- unlist(lapply(resamp_Y_list, LRclassify))

saveRDS(cbind(Y_list, y), 'event_seqs.rds')

# extract time series features from wave sequences
X <- tsfeatures(resamp_Y_list) %>% select(-frequency, -nperiods, -seasonal_period)

###
# logit regression
###

# 0/1 encode left or right strings
y_logit = ifelse(y == "Right", 1, 0)
X_logit <- as.matrix(X)

# cross-validated search for optimal lambda
lambdas = NULL
for (i in 1:50) {
  cv.out = cv.glmnet(X_logit, y_logit, family='binomial', type.measure='class')
  errors = data.frame(cv.out$lambda, cv.out$cvm)
  lambdas = rbind(lambdas, errors)
}

# select optimal lambda based on mse
lambdas_agg = aggregate(lambdas[,2], list(lambdas$cv.out.lambda), mean)
i = which(lambdas_agg[2] == min(lambdas_agg[2]))
error_dist = lambdas[which(lambdas[1] == lambdas[i,1]),2]
saveRDS(1-error_dist, 'logreg_accuracy.RDS') 
l = ifelse(length(lambdas_agg[i,1]) > 1, lambdas_agg[i,1][1], lambdas_agg[i,1])

# final model
logit_model = glmnet(X_logit, y_logit, family='binomial', lambda=cv.out$lambda.min)

# repeated 5-fold cross validation for individual left/right error classification
# over the whole data set
n = nrow(X_logit)
auc = c()
for (i in 1:50) {
  cvSets = cvTools::cvFolds(n, 5)  # permute all the data, into 5 folds
  
  for (j in 1:5) {
    test_id = cvSets$subsets[cvSets$which == j]
    X_test = X_logit[test_id, ]
    X_train = X_logit[-test_id, ]
    y_test = y_logit[test_id]
    y_train = y_logit[-test_id]
    
    model = glmnet(X_train, y_train, family='binomial', lambda=cv.out$lambda.min)
    res = predict(model, X_test, type='response')
    auc = c(auc, MLmetrics::AUC(res, y_test))
  }
}

saveRDS(f1, 'logreg_auc.rds')
    
# function to extract time series features from sequence
# and run logistic regression model to classify as left or right
leftright_classify <- function(seq) {
  features= tsfeatures(seq) %>% select(-frequency, -nperiods, -seasonal_period)
  prob = predict(logit_model, newx=as.matrix(features), type='response')
  return(ifelse(prob > 0.5, "Right", "Left"))
}

###
# streaming
###

simulate_streaming <- function(sequence, down_sample=20) {
  prev_buffer = c() # previous second
  sd_buffer = c() # event buffer
  events = list()
  j = 1
  
  # downsample time series to increase performance
  for (i in seq(1, length(sequence), by=down_sample)) {
    # always keep one second previous buffer
    if (length(prev_buffer) < 10000/down_sample) {
      prev_buffer = c(prev_buffer, sequence[i])
    } else {
      prev_buffer = c(prev_buffer[-1], sequence[i])
    }
    
    if (length(prev_buffer) > 1 & sd(prev_buffer) > 220) {
      if (length(sd_buffer) == 0) {
        # add previous second to event
        sd_buffer = c(prev_buffer, sequence[i])
      } else {
        sd_buffer = c(sd_buffer, sequence[i])
      }
    # sd has dropped back down so add event to list if there is one
    } else if (length(sd_buffer) > 0) {
      events[[j]] <- sd_buffer
      j = j + 1
      sd_buffer = c()
    }
  }
  
  # if stream ends while sd is above threshold, event could be missed
  if (length(sd_buffer) > 0) {
    events[[j]] <- sd_buffer
  }
  
  # predict labels for event and return sequence
  result = lapply(events, leftright_classify)
  return(result)
}

# predicted label sequences
sequences = lapply(wave_seqs, function(x) x@left)
results = lapply(sequences, simulate_streaming)
results = lapply(results, unlist)
results = lapply(results, function(x) ifelse(x == "Right", "R", "L"))
results = lapply(results, function(x) paste(x, collapse=''))

# truth label sequences
sequence_truth = lapply(wave_seq_events, function(x) unlist(lapply(x, LRclassify)))
sequence_truth = lapply(sequence_truth, function(x) ifelse(x == "Right", "R", "L"))
sequence_truth = lapply(sequence_truth, function(x) paste(x, collapse=''))

# compare sequence accuracies
comparison = data.frame(truth=unlist(sequence_truth), predicted=unlist(results))
max_lengths = pmax(str_length(comparison$truth), str_length(comparison$predicted))
distances = stringdist::stringdist(comparison$truth, comparison$predicted, method="lv")
comparison$acc = 1 - (distances / max_lengths) # normalise levenshtein distance

saveRDS(comparison, 'sequence_accuracy.rds')
