# First we define the initial age distribution n_vec, 
# which can contain any non-negative numbers.
# I'm calling it n_vec instead of n just to remind us that it's a vector.
# Note that in R, we use the letter c to create a vector.
n_vec <- c(1, 1, 1, 1, 1)

# Here is a bar plot of the initial age distribution:
barplot(t(n_vec))

# Next we define the population projection matrix A, 
# where the first row contains fecundity values (any non-negative numbers)
# and the diagonal below the main diagonal contains survival probabilities
# (which must be between 0 and 1). All other entries are zero.
A <- matrix(c(0, 0, 1, 2, 0.1, 
              0.8, 0, 0, 0, 0, 
              0, 0.7, 0, 0, 0, 
              0, 0, 0.6, 0, 0, 
              0, 0, 0, 0.5, 0), nrow = 5, byrow = TRUE)

# To find the age distribution at the next time step, 
# we multiply the previous age distribution vector (n_vec) 
# by the population projection matrix (A). 
# In R, matrix multiplication is written as %*%.
n_vec <- A %*% n_vec

# Here is a bar plot of new initial age distribution:
barplot(t(n_vec))

# If we repeat the multiplication 100 times then the 
# age distribution converges to the stable age distribution.
for(i in 1:100) n_vec <- A %*% n_vec

# Here is a bar plot of new initial age distribution:
barplot(t(n_vec))

# The other way to get the stable age distribution is by calculating
# the right eigenvector of A, corresponding to the dominant eigenvalue:
w <- Re(eigen(A)$vectors[, 1])
# We can normalize w to get the proportions in each age class
# (remember that any scalar multiple of an eigenvector is also an eigenvector).
w_rescaled <- w / sum(w)

# Here is a bar plot of stable age distribution:
barplot(w_rescaled)

# The reproductive value is the left eigenvector of A, 
# corresponding to the dominant eigenvalue:
v <- Re(eigen(t(A))$vectors[, 1])
# We need to rescale v so it can be interpreted as a reproductive value 
# (remember that any scalar multiple of an eigenvector is also an eigenvector).
# The rescaling factor is 1 / (t(v) * w), where t means transpose.
v_rescaled <- v / c(t(v) %*% w)

# Here is a bar plot of the reproductive value as a function of age:
barplot(v_rescaled)

# The population increases because the dominant eigenvalue 
# of matrix A is greater than 1:
Re(eigen(A)$values[1])

# Now try changing the entries of A to see the effects on 
# w, v and the population growth rate.





