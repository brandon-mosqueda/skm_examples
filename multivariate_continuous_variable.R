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
Lines <- model.matrix(~ 0 + Line, data = Pheno)
# Env design matrix without the first column
Envs <- model.matrix(~ 0 + Env, data = Pheno)[, -1]

Geno <- SKM::cholesky(Geno)
Geno <- Lines %*% Geno
# Interaction matrix
GenoxEnv <- model.matrix(~ 0 + Geno:Envs)

# Put the matrices in the expected format with the model
# to use with each one
X <- list(
  list(x = Envs, model = "FIXED"),
  list(x = Geno, model = "BRR"),
  list(x = GenoxEnv, model = "BRR")
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

for (fold_number in seq_along(folds)) {
  fold <- folds[[fold_number]]

  # Set testing indices to NA in the response variables
  y_na <- y
  y_na[fold$testing, ] <- NA

  # *********** MODEL EVALUATION AND HYPERPARAMETERS TUNING WITH SKM ***********
  model <- SKM::bayesian_model(
    X,
    y_na,

    iterations_number = 100,
    burn_in = 50
  )

  # Predict over testing set
  predictions <- predict(model)

  for (trait in names(Predictions)) {
    # Save the predictions along with environment and line information
    Predictions[[trait]] <- rbind(
      Predictions[[trait]],
      data.frame(
        Fold = fold_number,
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
