# Import SKM library
library(SKM)

# Load the dataset
load("data/EYT_1.RData", verbose = TRUE)

# Oder Pheno by Env and Lines
Pheno <- Pheno[order(Pheno$Env, Pheno$Line), ]
geno_order <- sort(rownames(Geno))
# Order Geno by Line
Geno <- Geno[geno_order, geno_order]

# Design matrices
Z_g <- model.matrix(~ 0 + Line, data = Pheno)
# Linear kernel
K_g <- Z_g %*% Geno %*% t(Z_g)

# Env design matrix without the first column
X_E <- model.matrix(~ 0 + Env, data = Pheno)[, -1]
K_e <- tcrossprod(X_E) / ncol(X_E)

X_g <- SKM::cholesky(K_g)
K_ge <- K_e * K_g
X_ge <- SKM::cholesky(K_ge)

# Bind all matrices in a single one
X <- cbind(X_E, X_g, X_ge)
# Retrieve the two mixed response variables
y <- Pheno[, c("DTHD", "GY")]

# Folds generation for stratified random cross validation
# it keeps the distribution of environment at each fold
folds <- SKM::cv_random_strata(
  data = Pheno$Env,
  folds_number = 5,
  testing_proportion = 0.25
)

# List for storing the individual predictions at each fold
Predictions <- list(DTHD = data.frame(), GY = data.frame())

for (fold_number in seq_along(folds)) {
  fold <- folds[[fold_number]]

  # Split training and testing data
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training, ]
  y_testing <- y[fold$testing, ]

  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::deep_learning(
    X_training,
    y_training,

    # Tunable hyperparameters
    learning_rate = list(min = 0.001, max = 0.1),
    epochs_number = 50,
    batch_size = 32,
    layers = list(
      list(
        neurons_number = list(min = 20, max = 50),
        activation = "sigmoid",
        dropout = list(min = 0.2, max = 0.4)
      ),
      list(
        neurons_number = list(min = 20, max = 50),
        activation = "tanh",
        dropout = list(min = 0.2, max = 0.4)
      )
    ),

    # Tune configuration parameters
    tune_type = "Bayesian_optimization",
    tune_cv_type = "Random",
    tune_testing_proportion = 0.2,
    tune_folds_number = 2,
    tune_bayes_samples_number = 5,
    tune_bayes_iterations_number = 5,

    # Other algorithm's parameters
    optimizer = "adam",
    shuffle = TRUE,
    early_stop = TRUE,
    early_stop_patience = 10
  )

  # Predict over testing set
  predictions <- predict(model, X_testing)

  for (trait in names(Predictions)) {
    # Save the predictions along with environment and line information
    Predictions[[trait]] <- rbind(
      Predictions[[trait]],
      data.frame(
        Fold = fold_number,
        Line = Pheno$Line[fold$testing],
        Env = Pheno$Env[fold$testing],
        Observed = y_testing[, trait],
        Predicted = predictions[[trait]]$predicted
      )
    )
  }
}

# Errors metrics for prediction performance
DTHDSummaries <- SKM::gs_summaries(Predictions$DTHD)

# Print summaries by line, environment and fold
DTHDSummaries$line
DTHDSummaries$env
DTHDSummaries$fold

# Errors metrics for prediction performance
GYSummaries <- SKM::gs_summaries(Predictions$GY)

# Print summaries by line, environment and fold
GYSummaries$line
GYSummaries$env
GYSummaries$fold
