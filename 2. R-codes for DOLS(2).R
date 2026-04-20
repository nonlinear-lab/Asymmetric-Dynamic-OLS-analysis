#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 1: set a dorking directory and library
#++++++++++++++++++++++++++++++++++++++++++++++++
setwd("C:/R/fmols")

# load the package
# install.packages("cointReg")
library(cointReg)
mout<- matrix(data=NA,nrow=4,ncol=1)
#++++++++++++++++++++++++++++++++++++++++++++++++
# Step 2: load data
#++++++++++++++++++++++++++++++++++++++++++++++++
my_data <- read.csv("data(4).csv")

for (j in 1:4) {
  if (j == 1) {
    y <- as.matrix(my_data[2:T, 3])
    x <- as.matrix(my_data[2:T, 2])
  } else if (j == 2) {
    y <- as.matrix(my_data[2:T, 4])
    x <- as.matrix(my_data[2:T, 2])
  } else if (j == 3) {
    set.seed(42)
    y <- as.matrix(my_data[2:T, 5])
    x <- as.matrix(my_data[2:T, 2])
  } else {
    set.seed(42)
    y <- as.matrix(my_data[2:T, 6])
    x <- as.matrix(my_data[2:T, 2])
  }
  #++++++++++++++++++++++++++++++++++++++++++++++++
  # Step 3: DOLS analysis
  #++++++++++++++++++++++++++++++++++++++++++++++++
  # 2. DOLS Estimation
  # Create a constant intercept
  intercept <- rep(1, length(y))

  # We use the cointRegD function. You need to specify the number of leads and lags.
  # Here we use 2 leads and 2 lags as an example.
  dols_model <- cointRegD(y = y, x = x, deter = intercept, n.lead = 2, n.lag = 2)
  dols   <- dols_model$theta[2]
  mout[j,1]   <- dols
  print(dols_model)
}
print(mout)
