# Install devtools if not installed
if (!require ("devtools")) {
  install.packages ("devtools")
}

# Install the most updated version of these dependencies
devtools::install_github("cran/randomForestSRC")
devtools::install_github("gdlc/BGLR-R")
devtools::install_github("rstudio/tensorﬂow")

# Install SKM
devtools::install_github("brandon-mosqueda/SKM")
