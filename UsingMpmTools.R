install.packages("devtools") # install a package that enables you to install mpmtools from the github website
library(devtools) # load that package

devtools::install_github("BruceKendall/mpmtools") # install the mpmtools package
library(mpmtools) # load the mpmtools package

########

Caswell_Ex_2.1 # an example matrix from Caswell's 2001 book

lambda1(Caswell_Ex_2.1) # dominant eigenvalue of that example matrix

stable_stage(Caswell_Ex_2.1) # eigenvector associated with the dominant eigenvector

########

# Create a demography schedule, with juvenile and senescent age classes:
demog_sched <- data.frame(x = 0:7,
                          sx = c(0.05, 0.2, 0.35, 0.8, 0.9, 0.9, 0.75, 0.4),
                          mx = c(0, 0, 0, 0.5, 1, 3, 3, 1.5))

# Construct a Leslie matrix from this demography schedule:
A1 <- make_Leslie_matrix(demog_sched)

lambda1(A1) # dominant eigenvalue

stable_stage(A1) # eigenvector associated with the dominant eigenvector
