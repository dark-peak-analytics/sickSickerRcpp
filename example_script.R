# This R script is a replication of 'vignette/run_microsim.rmd' for those
# who would rather run the analysis in an R script.

# It demonstrates the C++ functions defined in 'src/define_microsim.cpp'.

# The model used as an example is the sick-sicker model published elsewhere;
# Krijkamp et al. 2018. Microsimulation modeling for health decision sciences using R: a tutorial
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6349385/.

# We convert the model to C++ and uses the Rcpp package to run as R functions.
# This runs MUCH faster because:
# 1. all individuals run through the model simultaneously (vectorised)
# 2. C++ is faster than R for these operations.

# Of these two things (1) is probably more important!

# Please ensure that you open this project by launching the sickSickerRcpp.Rproj file.
#
# If this is not possible then:
#   
# 1. set the working directory to the project folder (sickSickerRcpp/)
# e.g. setwd("C:/Users/r_a_s/Documents/Projects/DPA_courses/sickSickerRcpp")
# 2. ensure the file 'src/define_microsim.cpp' is in the working directory
# "define_microsim.cpp" %in% list.files("src/") 
# 3. install the "Rcpp" and "RcppArmadillo" packages
#install.packages("Rcpp")
#install.packages("RcppArmadillo")

rm(list = ls())

# Use the `sourceCpp` function from the `Rcpp` package to source the functions 
# defined in the`.cpp` file. and make avalable as R functions:
# MicroSimV_Cpp SampleV_Cpp CostsV_Cpp EffsV_Cpp ProbsV_Cpp SickSickerMicroSim_Cpp
Rcpp::sourceCpp(
  file = paste0(
    here::here(),
    "/src",
    "/define_microsim.cpp"
  )
)

# define individual parameters
n.i   <- 100000                # number of simulated individuals
n.t   <- 30                    # time horizon, 30 cycles
v.n   <- c("H","S1","S2","D")  # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s   <- length(v.n)           # the number of health states
v.M_1 <- rep("H", n.i)         # everyone begins in the healthy state
d.c   <- d.e <- 0.03           # equal discounting of costs and QALYs by 3%
v.Trt <- c("No Treatment", "Treatment") # store the strategy names

# Transition probabilities (per cycle)
p.HD    <- 0.005               # probability to die when healthy
p.HS1   <- 0.15          	     # probability to become sick when healthy
p.S1H   <- 0.5           	     # probability to become healthy when sick
p.S1S2  <- 0.105         	     # probability to become sicker when sick
rr.S1   <- 3             	     # rate ratio of death in sick vs healthy
rr.S2   <- 10            	     # rate ratio of death in sicker vs healthy
r.HD    <- -log(1 - p.HD) 	   # rate of death in healthy
r.S1D   <- rr.S1 * r.HD  	     # rate of death in sick
r.S2D   <- rr.S2 * r.HD  	     # rate of death in sicker
p.S1D   <- 1 - exp(- r.S1D)    # probability to die in sick
p.S2D   <- 1 - exp(- r.S2D)    # probability to die in sicker

# Cost and utility inputs
c.H     <- 2000                # cost of remaining one cycle healthy
c.S1    <- 4000                # cost of remaining one cycle sick
c.S2    <- 15000               # cost of remaining one cycle sicker
c.Trt   <- 12000               # cost of treatment (per cycle)

u.H     <- 1                   # utility when healthy
u.S1    <- 0.75                # utility when sick
u.S2    <- 0.5                 # utility when sicker
u.Trt   <- 0.95                # utility when being treated

# Define starting health state, using numbers instead of characters to identify the health states:
# 1 = H, 2 = S1, 3 = S2, 4 = D
v_M_1 <- rep(1, n.i)

# Create a vector of transition probabilities:
t_p = c(p.HD, p.HS1, p.S1H, p.S1S2, p.S1D, p.S2D)
names(t_p) <- c("p.HD", "p.HS1", "p.S1H", "p.S1S2", "p.S1D", "p.S2D")

# Create a vector containing costs parameters:
c_vec = c(c.H, c.S1, c.S2, c.Trt)
names(c_vec) <- c("c.H", "c.S1", "c.S2", "c.Trt")

# Create a vector containing utilities parameters:
u_vec = c(u.H, u.S1, u.S2, u.Trt)
names(u_vec) <- c("u.H", "u.S1", "u.S2", "u.Trt")

# Run the sick sicker model for each policy:
## No treatment:
# store into an R list object
ls_no_trt_res <- MicroSimV_Cpp(
  v_S_t = v_M_1,
  t_P = t_p,
  v_C = c_vec,
  v_U = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cl = 1,
  d_dC = d.c,
  d_dE = d.e,
  b_Trt = FALSE,
  n_Seed = 1
)

## Treatment:
# store into an R list object
ls_trt_res <- MicroSimV_Cpp(
  v_S_t = v_M_1,
  t_P = t_p,
  v_C = c_vec,
  v_U = u_vec,
  n_I = n.i,
  n_S = n.s,
  n_T = n.t,
  n_Cl = 1,
  d_dC = d.c,
  d_dE = d.e,
  b_Trt = TRUE,
  n_Seed = 1
)

# Run incremental analysis
names(ls_no_trt_res)

## Change in mean costs:
delta_costs <- ls_trt_res$tc_hat - ls_no_trt_res$tc_hat
delta_costs

## Change in mean QALYs:
delta_effects <- ls_trt_res$te_hat - ls_no_trt_res$te_hat
delta_effects

## ICER:
delta_costs / delta_effects

## Net benefits @ £30,000:
lambda <- 30000

### Treatment NB:
ls_trt_res$te_hat * lambda - ls_trt_res$tc_hat

### No treatment NB:
ls_no_trt_res$te_hat * lambda - ls_no_trt_res$tc_hat

### Incremental NB:
delta_effects * lambda - delta_costs