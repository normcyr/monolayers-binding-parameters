# R script to calculate the binding parameters (and corresponding errors)
# MIP and synergy of proteins upon binding to lipid monolayers
#
# References for the calculations:
# - Calvez et al. (2009) Biochimie 91:718.
# - Calvez et al. (2011) Langmuir 27:1373.
# - Boisselier et al. (2017) Advances in Colloid and Interface Science 243:60.
#
# Web server using these calculations:
# http://www.crchudequebec.ulaval.ca/BindingParametersCalculator/


# Import data from CSV file
df <- read.table("example-data.csv", header = TRUE, sep = ",")

# Create the standard error function
se <- function(x) sqrt(var(x)/length(x))

# Perform the linear regression 
linear_regression <- lm(df$Delta_Pi ~ df$Pi_i)

# Generate a plot for visual representation
plot(df$Pi_i, df$Delta_Pi,
     main = "ΔΠ vs Πi plot", type = "p", pch = 19, cex = 1.5,
     xlab = "Πi (mN/m)", ylab = "ΔΠ (mN/m)",
     xlim = c(0,40), ylim = c(0,25))
abline(linear_regression)

# Get the MIP and synergy
a <- coef(linear_regression)[["df$Pi_i"]]
b <- coef(linear_regression)[["(Intercept)"]]

MIP <- -b / a
synergy <- a + 1

# Basic stats on the data
n <- length(df$Pi_i)
r2 <- (cor(df$Delta_Pi, df$Pi_i))^2
sd_x <- sd(df$Pi_i)
sd_y <- sd(df$Delta_Pi)
ave_x <- mean(df$Pi_i)

sum_sq_x <- sum(df$Pi_i^2)
sum_x_ave_x_sq <- sum((df$Pi_i - ave_x)^2)

# Stats on the linear regression
se2 <- ((n - 1) * sd_y^2 * (1 - r2)) / (n - 2)

var_a <- vcov(linear_regression)[2,2]
var_b <- vcov(linear_regression)[1,1]
cov_ab <- (-ave_x * se2) / sum_x_ave_x_sq

# Determination of the confidence interval on the x-intercept (MIP)
var_I <- ( ((b / a)^2) * ((var_a / a^2) + (var_b / b^2) - ((2 * cov_ab) / (a * b))))
CI_MIP <- 1.96 * sqrt(var_I)

# Determination of the confidence interval on the synergy
CI_synergy = (sd_y * sqrt(1 - r2)) / (sd_x * sqrt(n - 2))
