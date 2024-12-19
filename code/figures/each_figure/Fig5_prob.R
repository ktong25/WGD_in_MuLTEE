##### Calculate the probability of finding k changes of a specific chromosome in n total changes

# Binomial probability
# Equivalent to the binomial probability of randomly selecting from 16 objects n times and getting k times of object A (each time p = 1/16)
n <- 53
p <- 1/16
k <- 1:11
probabilities <- dbinom(k, n, p)
names(probabilities) <- k
probabilities
# 1            2            3            4            5            6            7            8 
# 0.1155213511 0.2002370086 0.2269352764 0.1891127304 0.1235536505 0.0658952803 0.0294959826 0.0113067933 
# 9           10           11 
# 0.0037689311 0.0011055531 0.0002881138 
p <- 1/32  # if consider gain/loss of a chromosome as separate events
probabilities <- dbinom(k, n, p)
names(probabilities) <- k
probabilities
# 1            2            3            4            5            6            7            8 
# 3.177862e-01 2.665304e-01 1.461618e-01 5.893622e-02 1.863145e-02 4.808116e-03 1.041389e-03 1.931609e-04 
# 9           10           11 
# 3.115499e-05 4.421998e-06 5.576127e-07 
prob_11_or_more <- sum(dbinom(11:53, n, p))  # dbinom(11, n, p) # sum(dbinom(0:53, n, p))
# 6.276178e-07
prob_10_or_more <- sum(dbinom(10:53, n, p))  # dbinom(10, n, p)
# 5.049616e-06

##### Probability of observing six or more pairs with identical karyotype changes, each with changes in no more than two chromosomes, assuming equal probability for gains and losses of 16 chromosomes
# We have 16 different lights, each light may change to red or white (thus 2 colors) with equal probability. 
# We have 23 streets each with 16 different lights. Now we let each street change its light colors and each street can change no more than 2 lights. 
# What is the probability of observing 6 pairs of streets, where each pair changes in the same way?

# Calculate the number of ways each strain can change
ways_to_change_strain <- 1 + (choose(16, 1) * 2) + (choose(16, 2) * 2^2)  # change no more than two chromosomes
prob_same_change <- 1 / ways_to_change_strain
prob_not_same_change <- (ways_to_change_strain - 1) / ways_to_change_strain

# Function to calculate probability for exactly k pairs
calculate_probability_for_k_pairs <- function(k) {
  strains_involved <- 2 * k
  remaining_strains <- 23 - strains_involved
  combinations_for_k_pairs <- 1
  
  for (i in seq(k, 1, by=-1)) {
    combinations_for_k_pairs <- combinations_for_k_pairs * (choose(remaining_strains + 2 * i, 2) / 2)
  }
  
  if (remaining_strains > 1) {
    other_pairs_k <- choose(remaining_strains, 2)
  } else {
    other_pairs_k <- 0
  }
  
  probability_k_pairs <- combinations_for_k_pairs * (prob_same_change^k) * (prob_not_same_change^other_pairs_k)
  return(probability_k_pairs)
}

# Calculate probability for 6 or more pairs
prob_6_or_more <- 0

# We start from 6 pairs up to 11 pairs
for (k in 6:11) {
  prob_6_or_more <- prob_6_or_more + calculate_probability_for_k_pairs(k)
}

# Print the total probability
print(prob_6_or_more)
# 8.241817e-06
