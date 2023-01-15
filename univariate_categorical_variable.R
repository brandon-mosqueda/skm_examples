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
# Retrieve the categorical response variable
y <- Pheno$GY

# Folds generation for leave one environmet out cross validation
folds <- SKM::cv_leve_one_group_out(as.character(Pheno$Env))

# Data frame for storing the individual predictions at each fold
Predictions <- data.frame()

for (fold_number in seq_along(folds)) {
  fold <- folds[[fold_number]]

  # Split training and testing data
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training]
  y_testing <- y[fold$testing]

  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::generalized_linear_model(
    X_training,
    y_training,

    alpha = c(0, 0.5, 1),

    tune_type = "Grid_search"
  )

  # Predict over testing set
  predictions <- predict(model, X_testing)

  # Save the predictions along with environment and line information
  Predictions <- rbind(
    Predictions,
    data.frame(
      Fold = fold_number,
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = y_testing,
      Predicted = predictions$predicted,
      predictions$probabilities
    )
  )
}

# Errors metrics for prediction performance
Summaries <- SKM::gs_summaries(Predictions)

# Print summaries by line, environment and fold
Summaries$line
Summaries$env
Summaries$fold
