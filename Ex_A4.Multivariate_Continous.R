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

# Put the matrices in the expected format with the model
# to use with each one
X <- list(
  list(x =X_E, model = "FIXED"),
  list(x =K_g, model = "BRR"),
  list(x =K_ge, model = "BRR")
)
# Retrieve the two continuos response variables
y <- SKM::to_matrix(Pheno[, c("DTHD", "DTMT")])

# Folds generation for random cross validation
folds <- SKM::cv_random(
  records_number = nrow(Pheno),
  folds_number = 5,
  testing_proportion = 0.25
)

# List for storing the individual predictions at each fold
Predictions <- list(DTHD = data.frame(), DTMT = data.frame())

for (i in seq_along(folds)) {
  fold <- folds[[i]]

  # Set testing indices to NA in the response variables
  y_na <- y
  y_na[fold$testing, ] <- NA

  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::bayesian_model(
    X,
    y_na,

    iterations_number = 10000,
    burn_in = 5000
  )

  # Predict over testing set
  predictions <- predict(model)

  for (trait in names(Predictions)) {
    # Save the predictions along with environment and line information
    Predictions[[trait]] <- rbind(
      Predictions[[trait]],
      data.frame(
        Fold = i,
        Line = Pheno$Line[fold$testing],
        Env = Pheno$Env[fold$testing],
        Observed = y[fold$testing, trait],
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
DTMTSummaries <- SKM::gs_summaries(Predictions$DTMT)

# Print summaries by line, environment and fold
DTMTSummaries$line
DTMTSummaries$env
DTMTSummaries$fold
