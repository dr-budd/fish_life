## SET UP ----

## data_set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## load libraries
library(caret) ## used to split data evenly
library(glmnetUtils) ## for cva.glmnet
library(glmnet) ## contains the nested_elastic net regression model
library(cocor) ## compare correlations
library(plotmo) ## for plotres
library(tidyverse) ## helpful for data handling and visualization

## specify function to extract the parameters 
## from the best model (by Daniel Freeman)
get_model_params <- function(alpha_crossval_training) {
  alpha <- alpha_crossval_training$alpha
  lambdaMin <- sapply(alpha_crossval_training$modlist, `[[`, "lambda.min")
  lambdaSE <- sapply(alpha_crossval_training$modlist, `[[`, "lambda.1se")
  error <- sapply(alpha_crossval_training$modlist, function(mod) {min(mod$cvm)})
  best <- which.min(error)
  data.frame(alpha = alpha[best], lambdaMin = lambdaMin[best],
             lambdaSE = lambdaSE[best], eror = error[best])
}

## IMPORT DATA ----

## import summary metadata 
metadata_sum <- read.csv("dataFiles/06.00_metadata_sum.csv") %>%
  ## filter for lifespan species
  filter(!is.na(mean_lifespan))

## import CpG oe data
cpg_oe_data_unmod <- read.csv("dataFiles/06.00_cpg_oe_data.csv") %>%
  column_to_rownames(., var="X") %>%
  ## replace NAs with zeros
  replace(is.na(.), 0) %>%
  ## format as matrix
  as.matrix(.)

## filter for lifespan species
cpg_oe_data <- cpg_oe_data_unmod[metadata_sum$organism_name, ]

## IMPORT FC DATA ----

## import CpG O/E data (fc)
fc_oe <- read.csv("dataFiles/06.00_cpg_oe_data_fc.csv") %>%
  column_to_rownames(., var="X") %>%
  ## replace NAs with zeros
  replace(is.na(.), 0) %>%
  ## format as matrix
  as.matrix(.)

## metadata summary for "five copies" genomes (fc)
meta_sum_fc <- read.csv("dataFiles/06.00_metadata_sum_fc.csv")

## NESTED ELASTIC NET ----

## choose a 'lifespan type'
ls_type <- "mean_lifespan"

## print it
paste0("lifespan type specified: ", ls_type)

## STEP 1: split the data into training (70 %) and testing (30 %) data sets ----

## create blank files to store data
nested_elastic_parameters <- NULL
nested_elastic_correlations <- NULL
nested_elastic_metadata_summary <- NULL
nested_elastic_metadata_summary_fc <- NULL

## create table containing all promoter ids to store coefficient information
nested_promoter_coefficients <- data.frame(promoter_id = colnames(cpg_oe_data))

## set seed for create data partition function
set.seed(1)
  
## use the createDataPartition function to split the data 70/30
## 10 times based on the percentiles of y
sample_index_list <- createDataPartition(y = log(metadata_sum[[ls_type]]), 
                                    times = 10, 
                                    p = 0.7,
                                    list = TRUE) 

## create counter
partition_number <- 0

for (sample_index in sample_index_list) {
  
  ## add counter
  partition_number <- partition_number + 1

  ## use the sample_index to take the selected 70% (training)
  training_oe <- data.matrix(cpg_oe_data[sample_index, ])
  
  ## use the inverse sample_index to select the remaining 30% (testing)
  testing_oe <- data.matrix(cpg_oe_data[-sample_index, ])
  
  ## do the same for the phenotype data 
  training_meta <- metadata_sum[sample_index, ]
  testing_meta <- metadata_sum[-sample_index, ]
  
  ## STEP 2: perform cross validation ----
    
  ## set seed for glmnet cross validation 
  set.seed(1234)
  
  ## perform 10-fold cross validation for glmnet to determine optimal alpha and lambda values
  ## (using the training data)
  alpha_crossval_training <- cva.glmnet(x = training_oe,
                                        y = log(training_meta[[ls_type]]),
                                        nfolds = 10, 
                                        family = "gaussian")
  
  ## get model parameters for the best alpha (the alpha that produces the minimum error)
  alpha_min <- get_model_params(alpha_crossval_training)$alpha
  lambda_min <- get_model_params(alpha_crossval_training)$lambdaMin
  lambda_1se <- get_model_params(alpha_crossval_training)$lambdaSE
  mse_min <- get_model_params(alpha_crossval_training)$eror
  
  ## STEP 3: make the model ----
  
  ## fit a glm with elastic net regularisation to the training data
  glmnet_model <- glmnet(training_oe,
                            log(training_meta[[ls_type]]), 
                            family="gaussian", 
                            nlambda=100, 
                            alpha=alpha_min)
  
  ## extract the coefficients for the specified value of lambda
  promoters <- as.matrix(coef(glmnet_model, s=lambda_1se))
  
  ## select the non-zero coefficients 
  promoter_coefficient <- data.frame(promoters) %>%
    rename(!!paste0("partition_", partition_number) := s1) %>%
    rownames_to_column(., var = "promoter_id")
  
  ## count the number of predictive promoters
  num_pred_proms <- nrow(data.frame(promoters) %>%
                           filter(s1 != 0))
  
  ## STEP 4: use the model to predict the data (to determine its performance) ----
  
  ## TRAINING ----
  
  ## use the model and the minimum lambda value to predict log(lifespan) for the training data
  training_meta$log_predicted_lifespan <- predict(glmnet_model, 
                              training_oe, 
                              type="response", 
                              s=lambda_1se)[,'s1']
  
  ## detransform the predicted value
  training_meta$detrans_predicted_lifespan <- exp(training_meta$log_predicted_lifespan)
  
  ## transform the input value
  training_meta$log_lifespan <- log(training_meta[[ls_type]])
  
  ## calculate residuals
  training_meta$log_residuals <- training_meta$log_lifespan - training_meta$log_predicted_lifespan
  
  ## determine the pearson correlation in the training data set
  pearson_training <- cor.test(x = training_meta$log_lifespan, 
                                   y = training_meta$log_predicted_lifespan, 
                                   method = "pearson")

  
  ## create vector for lower and upper conf. intervals
  conf_cols <- c("conf_0", "conf_1")
  
  ## create temporary table containing pearson cor info (training)
  ## NB: this will produce lots of row.names related warnings...
  temp_training_table <- cbind.data.frame(pearson_training[], row.names = NULL) %>%
    ## add col for conf
    cbind(., conf_cols) %>%
    ## add mean squared error 
    mutate(mse = mean((training_meta$log_lifespan - training_meta$log_predicted_lifespan)^2)) %>%
    ## add mean absolute error
    mutate(mae = mean(abs(training_meta$log_lifespan - training_meta$log_predicted_lifespan))) %>%
    ## make 1 row
    pivot_wider(., names_from = "conf_cols",
      values_from = "conf.int") %>%
    ## add suffixes
    rename_with( ~ paste0(.x, "_train"))
  
  ## TESTING ----
  
  ## use the model and the minimum lambda value to predict log(lifespan) for the testing data
  testing_meta$log_predicted_lifespan <- predict(glmnet_model, 
                             testing_oe, 
                             type="response", 
                             s=lambda_1se)[,'s1']
  
  ## detransform the predicted value
  testing_meta$detrans_predicted_lifespan <- exp(testing_meta$log_predicted_lifespan)

  ## transform the input value
  testing_meta$log_lifespan <- log(testing_meta[[ls_type]])
  
  ## calculate residuals
  testing_meta$log_residuals <- testing_meta$log_lifespan - testing_meta$log_predicted_lifespan
  
  ## determine the pearson correlation in the testing data set
  pearson_testing <- cor.test(x = testing_meta$log_lifespan,
                                   y = testing_meta$log_predicted_lifespan,
                                   method = "pearson")

  ## create vector for lower and upper conf. intervals
  conf_cols <- c("conf_0", "conf_1")
  
  ## create temporary table containing pearson cor info (testing)
  temp_testing_table <- cbind.data.frame(pearson_testing[], row.names = NULL) %>%
    ## add col for conf
    cbind(., conf_cols) %>%
    ## add mean squared error 
    mutate(mse = mean((testing_meta$log_lifespan - testing_meta$log_predicted_lifespan)^2)) %>%
    ## add mean absolute error
    mutate(mae = mean(abs(testing_meta$log_lifespan - testing_meta$log_predicted_lifespan))) %>%
    ## make 1 row
    pivot_wider(., names_from = "conf_cols",
                values_from = "conf.int") %>%
    ## add suffixes
    rename_with( ~ paste0(.x, "_test"))
  
  ## make another df with predicted lifespans only (testing only)
  nested_elastic_lifespan <- testing_meta %>%
    select(organism_name, detrans_predicted_lifespan) %>%
    rename(!!paste0(partition_number) := detrans_predicted_lifespan)
  
  ## FOR GENOMES WITH 5 COPIES (FC) ----
  
  ## use the model and the minimum lambda value to predict log(lifespan) for the training data
  meta_sum_fc$log_predicted_lifespan <- predict(glmnet_model, 
                                                fc_oe, 
                                                type="response", 
                                                s=lambda_1se)[,'s1']
  
  ## detransform the predicted value
  meta_sum_fc$detrans_predicted_lifespan <- exp(meta_sum_fc$log_predicted_lifespan)
  
  ## transform the input value
  meta_sum_fc$log_lifespan <- log(meta_sum_fc[[ls_type]])
  
  ## calculate residuals
  meta_sum_fc$log_residuals <- meta_sum_fc$log_lifespan - meta_sum_fc$log_predicted_lifespan
  
  ## add partition number
  meta_sum_fc$partition_number = partition_number
  
  ## create metadata
  nested_elastic_metadata_summary_fc <- rbind(nested_elastic_metadata_summary_fc,
                                              meta_sum_fc)
  
  ## COMPARE TESTING AND TRAINING CORRELATIONS ----
  
  ## t-test residuals (absolute) and pull significance code
  sig_diff_resid <- as.data.frame(t.test(abs(training_meta$log_residuals), 
                                         abs(testing_meta$log_residuals))[], 
                                      row.names=NULL)$p.value[1]
  
  ## compare correlations 
  
  ## create testing df
  cocor_test <- cbind.data.frame(act_ls = testing_meta$log_lifespan, 
                                     testing_meta$log_predicted_lifespan) %>%
    rename(pred_ls = 1)
  
  ## create training df
  cocor_train <- cbind.data.frame(act_ls = training_meta$log_lifespan, 
                                      training_meta$log_predicted_lifespan) %>%
    rename(pred_ls = 1)
  
  ## create vars as per manual
  n1 <- nrow(training_meta)
  n2 <- nrow(testing_meta)
  
  r1.jk <- pearson_training$estimate[[1]]

  r2.hm <- pearson_testing$estimate[[1]]

  ## compare correlations and pull p-value
  sig_diff_cor <- cocor.indep.groups(r1.jk, r2.hm, n1, n2, 
                                     data.name=c("cocor_train", "cocor_test"),
                                     var.labels=c("act_ls", "pred_ls", 
                                                  "act_ls", "pred_ls"))@fisher1925$p.value
  
  ## create a df with model parameters
  temp_parameters <- cbind(ls_type, 
                           partition_number, 
                           alpha_min, 
                           lambda_min, 
                           lambda_1se, 
                           mse_min, 
                           num_pred_proms)
  
  ## combine your tables and the lifespan and seed info
  temp_table <- cbind(ls_type, 
                      partition_number, 
                      temp_training_table,
                      temp_testing_table, 
                      sig_diff_resid, 
                      sig_diff_cor)
  
  ## append for model
  nested_elastic_parameters <- rbind(temp_parameters, nested_elastic_parameters)
  
  ## append for model coefficients
  nested_promoter_coefficients <- left_join(promoter_coefficient, 
                                     nested_promoter_coefficients, 
                                     by = "promoter_id")
  
  ## append for correlations
  nested_elastic_correlations <- rbind(temp_table, nested_elastic_correlations)
  
  ## and for the meta
  training_meta$data_set <- "training"
  testing_meta$data_set <- "testing"
  sum_meta <- rbind(training_meta, testing_meta) %>%
    mutate(partition_number = partition_number)
  
  nested_elastic_metadata_summary <- rbind(nested_elastic_metadata_summary, 
                                           sum_meta)
  
  # ## MAKE CV PLOTS (not used) ----
  # 
  # plot_name <- paste0("figuresTables/IcvPlots/icv_plots_m", partition_number, ".pdf")
  # pdf(plot_name, width = 11.75, height = 7.25)
  # 
  # ## create a 2 x 2 plotting matrix
  # par(mfrow = c(2, 2),
  #     ## reduce margin size
  #     mar=c(3.5,3.5,3.5,3.5),
  #     ## reduce space between axis title, axis labels and axis ticks
  #     mgp=c(2,0.5,0))
  # 
  # ## print a summary of the cv results (including lambda min and labmda 1se)
  # plot(alpha_crossval_training, ann=FALSE)
  # mtext(text="A", xpd=NA, side=3, adj=0, font=2, line=1) ## A, B, C, D label
  # mtext(text="log(Lambda)", xpd=NA, side=1, line=2, font=1, cex=0.9) ## x axis
  # mtext(text="Mean-squared error", xpd=NA, side=2, line=2, font=1, cex=0.9) ## y axis
  # 
  # ## print a summary of the cv results (including lambda min and labmda 1se)
  # minlossplot(alpha_crossval_training, cv.type = c("1se"),
  #             bg="blue", pch=21,
  #             type = "b", ann=FALSE)
  # mtext(text="B", xpd=NA, side=3, adj=0, font=2, line=1)
  # mtext(text="Alpha", xpd=NA, side=1, line=2, font=1, cex=0.9) ## x axis
  # mtext(text="Mean-squared error", xpd=NA, side=2, line=2, font=1, cex=0.9) ## y axis
  # 
  # ## plot coefficients
  # plotres(glmnet_model, predict.s = lambda_1se, which =1, ann=FALSE)
  # mtext(text="C", xpd=NA, side=3, adj=0, font=2, line=1)
  # mtext(text="log(Lambda)", xpd=NA, side=1, line=2, font=1, cex=0.9) ## x axis
  # mtext(text="Coefficients", xpd=NA, side=2, line=2, font=1, cex=0.9) ## y axis
  # 
  # ## plot residuals
  # plotres(glmnet_model, main="",
  #         predict.s = lambda_1se, which=3, id.n=0, ann=FALSE)
  # mtext(text="D", xpd=NA, side=3, adj=0, font=2, line=1)
  # mtext(text="Fitted", xpd=NA, side=1, line=2, font=1, cex=0.9) ## x axis
  # mtext(text="Residuals", xpd=NA, side=2, line=2, font=1, cex=0.9) ## y axis
  # 
  # ## close the file
  # dev.off()
    
}

## export correlations, coefficients and lifespans
write.csv(nested_elastic_parameters, "dataFiles/08.00_nested_elastic_parameters.csv", row.names = FALSE)
write.csv(nested_elastic_correlations, "dataFiles/08.00_nested_elastic_correlations.csv", row.names = FALSE)
write.csv(nested_promoter_coefficients, "dataFiles/08.00_nested_promoter_coefficients.csv", row.names = FALSE)
write.csv(nested_elastic_metadata_summary, "dataFiles/08.00_nested_elastic_metadata_summary.csv", row.names = FALSE)

## fc
write.csv(meta_sum_fc, "dataFiles/08.00_nested_elastic_metadata_summary_fc.csv", row.names = FALSE)

## print session info
sessionInfo()

## END SCRIPT
##

