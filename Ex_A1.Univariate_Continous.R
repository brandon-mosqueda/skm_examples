# Import SKM library
library(SKM)

# Load the dataset
load("EYT_1.RData", verbose = TRUE)

# Oder Pheno by Env and Lines
Pheno <- Pheno[order(Pheno$Env, Pheno$Line), ]
geno_order <- sort(rownames(Geno))
# Order Geno by Line
Geno <- Geno[geno_order, geno_order]

# Design matrices
Z_g<- model.matrix(~ 0 + Line, data = Pheno)
K_g <- Z_g %*% Geno %*% t(Z_g) # Linear kernel

# Env design matrix without the first column
X_E<- model.matrix(~ 0 + Env, data = Pheno)[, -1]
K_e <- X_E %*% t(X_E) / ncol(X_E)

X_g <- SKM::cholesky(K_g)
K_ge <- K_e * K_g
# Interaction matrix
X_ge<- SKM::cholesky(K_ge)

# Bind all matrices in a single one
X <- cbind(X_E, X_g, X_ge)

# Retrieve the continuos response variable
y <- as.numeric(Pheno$GY)

# Folds generation for k-fold cross validation
folds <- SKM::cv_kfold(nrow(Pheno), k = 5)

# Data frame for storing the individual predictions at each fold
Predictions <- data.frame()

for (i in seq_along(folds)) {
  fold <- folds[[i]]

  # Split training and testing data
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training]
  y_testing <- y[fold$testing]

  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::random_forest(
    X_training,
    y_training,
    trees_number = c(100, 200),
    node_size = c(5, 2),

    tune_type = "Grid_search",
    tune_cv_type = "k_fold",
    tune_folds_number = 5,
    tune_loss_function = "mse"
  )

  # Predict over testing set
  predictions <- predict(model, X_testing)

  # Save the predictions along with environment and line information
  Predictions <- rbind(
    Predictions,
    data.frame(
      Fold = i,
      Line = Pheno$Line[fold$testing],
      Env = Pheno$Env[fold$testing],
      Observed = y_testing,
      Predicted = predictions$predicted
    )
  )
}

# Errors metrics for prediction performance
Summaries <- SKM::gs_summaries(Predictions)

# Print summaries by line, environment and fold
Summaries$line
Summaries$env
Summaries$fold
