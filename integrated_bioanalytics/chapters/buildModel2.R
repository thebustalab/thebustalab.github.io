#### buildModel2

    buildModel2 <- function(
        data,
        model_type = c(
            "linear_regression", "random_forest_regression",
            "random_forest_classification", "logistic_regression",
            "contrastive_learning"
        ),
        input_variables = NULL,
        output_variable = NULL,
        fold_cross_validation = 10,
        optimization_parameters = list(
            n_vars_tried_at_split = seq(10,20,5),
            n_trees = seq(100,200,50),
            min_leaf_size = 1
        )
    ) {

        ## Prepare the data

            ## Correct the formula and the data, start the output
                model_type <- match.arg(model_type)
                if (is.null(input_variables) || length(input_variables) == 0) {
                    stop("Please provide at least one input variable.")
                }
                if (is.null(output_variable) || !nzchar(output_variable)) {
                    stop("Please provide a valid output variable.")
                }
                if (!output_variable %in% colnames(data)) {
                    stop("The specified output variable is not present in the data.")
                }
                if (!all(input_variables %in% colnames(data))) {
                    missing_inputs <- input_variables[!input_variables %in% colnames(data)]
                    stop(
                        paste(
                            "The following input variables are not present in the data:",
                            paste(missing_inputs, collapse = ", ")
                        )
                    )
                }
                if (!is.numeric(fold_cross_validation) || fold_cross_validation < 2) {
                    stop("fold_cross_validation must be a numeric value greater than or equal to 2.")
                }
                if (model_type %in% c("random_forest_regression", "random_forest_classification")) {
                    required_optimization_args <- c("n_vars_tried_at_split", "n_trees", "min_leaf_size")
                    if (!all(required_optimization_args %in% names(optimization_parameters))) {
                        stop("optimization_parameters must include n_vars_tried_at_split, n_trees, and min_leaf_size.")
                    }
                    if (any(lengths(optimization_parameters[required_optimization_args]) == 0)) {
                        stop("optimization_parameters entries must contain at least one value.")
                    }
                }
                data <- as.data.frame(data)
                if( any(is.na(data)) ) { stop("Yo - there are missing data in your data set, please deal with that before proceeding.") }
                formula <- paste0("`", output_variable, "` ~ `", paste(input_variables, collapse = "` + `"), "`")
                output <- list()
                output$model_type <- model_type

            ## Various checkes
                    if (dim(data)[1] < length(input_variables)) {
                        stop("Heyyyyy, there are fewer observations in the data than input variables. You should change the number of input variables or make more observations.")
                    }
                    if (!is.numeric(data[[output_variable]]) && !model_type %in% c("random_forest_classification")) {
                        stop("Uh oh! You cannot predict the value of a categorical variable with a regression model. Choose a classification model instead.")
                    }
                    if (model_type == "random_forest_classification" && !is.factor(data[[output_variable]])) {
                        data[[output_variable]] <- as.factor(data[[output_variable]])
                    }

        ## Modeling

            if (model_type == "linear_regression") {

                ## Run the fit and start the output
                    fit <- lm(
                        data = data,
                        formula = formula,
                        x = TRUE, y = TRUE, model = TRUE, qr = TRUE
                    )
                    output$model <- fit

                    input_y <- fit$y[order(as.numeric(names(fit$y)))]
                    total_sum_squares <- sum((input_y - mean(input_y, na.rm = TRUE))^2, na.rm = TRUE)
                    residual_sum_squares <- sum((summary(fit)$residuals)^2, na.rm = TRUE)

                    metrics <- rbind(
                        data.frame(
                            variable = c("r_squared", "total_sum_squares", "residual_sum_squares"),
                            value = c(summary(fit)$r.squared, total_sum_squares, residual_sum_squares),
                            std_err = NA,
                            type = "statistic",
                            p_value = NA,
                            p_value_adj = NA
                        ),
                        data.frame(
                            variable = names(coefficients(fit)),
                            value = as.numeric(coefficients(fit)),
                            std_err = round(as.numeric(summary(fit)$coefficients[,2]), 4),
                            type = "coefficient",
                            p_value = round(summary(fit)$coefficients[,4], 8),
                            p_value_adj = p.adjust(round(summary(fit)$coefficients[,4], 8))
                        )
                    )
                    metrics$value <- round(as.numeric(metrics$value), 4)
                    rownames(metrics) <- NULL
                    output$metrics <- metrics
            }

            if (model_type == "random_forest_regression") {

                ## Create a workflow that has a recipe and a model
                    workflow() %>%
                        add_recipe(
                            recipe( as.formula(formula), data = data ) # %>%
                            # step_normalize(all_numeric()) # %>% step_impute_knn(all_inputs())
                        ) %>%
                        add_model(
                            rand_forest() %>% # specify that the model is a random forest
                            set_args(mtry = tune(), trees = tune(), min_n = tune()) %>% # specify that the `mtry` and `trees` parameters needs to be tuned
                            set_engine("ranger", importance = "impurity") %>% # select the engine/package that underlies the model
                            set_mode("regression") # choose either the continuous regression or binary classification mode
                        ) -> workflow_for_training

                ## Tune the model on the training set
                    tune_results <- tune_grid(
                        workflow_for_training,
                        resamples = vfold_cv(data, v = fold_cross_validation), #CV object
                        grid = expand.grid(
                            mtry = optimization_parameters$n_vars_tried_at_split,
                            trees = optimization_parameters$n_trees,
                            min_n = optimization_parameters$min_leaf_size
                        ), # grid of values to try
                        metrics = metric_set(rmse) # metrics we care about
                    )

                # Check model parameters if you want
                    output$metrics <- collect_metrics(tune_results)
                    names(output$metrics)[names(output$metrics) == "mtry"] <- "n_vars_tried_at_split"
                    names(output$metrics)[names(output$metrics) == "trees"] <- "n_trees"
                    names(output$metrics)[names(output$metrics) == "min_n"] <- "min_leaf_size"
                    if ("n" %in% names(output$metrics)) {
                        output$metrics$fold_cross_validation <- output$metrics$n
                        output$metrics$n <- NULL
                    }

                ## Apply the best parameters to the workflow to creat the output model
                    output$model <- fit(
                        finalize_workflow(workflow_for_training, select_best(tune_results, metric = "rmse")),
                        data
                    )
            }

            if (model_type == "random_forest_classification") {

                ## Create a workflow that has a recipe and a model
                    workflow() %>%
                        add_recipe(
                            recipe( as.formula(formula), data = data ) # %>%
                            # step_normalize(all_numeric()) # %>% step_impute_knn(all_inputs())
                        ) %>%
                        add_model(
                            rand_forest() %>% # specify that the model is a random forest
                            set_args(mtry = tune(), trees = tune(), min_n = tune()) %>% # specify that the `mtry` and `trees` parameters needs to be tuned
                            set_engine("ranger", importance = "impurity") %>% # select the engine/package that underlies the model
                            set_mode("classification") # choose either the continuous regression or binary classification mode
                        ) -> workflow_for_training

                ## Tune the model on the training set
                    tune_results <- tune_grid(
                        workflow_for_training,
                        resamples = vfold_cv(data, v = fold_cross_validation), #CV object
                        grid = expand.grid(
                            mtry = optimization_parameters$n_vars_tried_at_split,
                            trees = optimization_parameters$n_trees,
                            min_n = optimization_parameters$min_leaf_size
                        ), # grid of values to try
                        metrics = metric_set(accuracy) # metrics we care about
                    )

                # Check model parameters if you want
                    output$metrics <- collect_metrics(tune_results)
                    names(output$metrics)[names(output$metrics) == "mtry"] <- "n_vars_tried_at_split"
                    names(output$metrics)[names(output$metrics) == "trees"] <- "n_trees"
                    names(output$metrics)[names(output$metrics) == "min_n"] <- "min_leaf_size"
                    if ("n" %in% names(output$metrics)) {
                        output$metrics$fold_cross_validation <- output$metrics$n
                        output$metrics$n <- NULL
                    }
                    # select_best(tune_results, metric = "rmse")

                ## Apply the best parameters to the workflow to creat the output model
                    output$model <- fit(
                        finalize_workflow(workflow_for_training, select_best(tune_results, metric = "accuracy")),
                        data
                    )
            }

        ## Return a model trained on all the input data
            
            return(output)
    }
