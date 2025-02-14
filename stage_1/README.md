# 1. Function for translating DNA to Protein
translate_dna <- function(dna_sequence) {
  codon_table <- list(
    "ATA" = "I", "ATC" = "I", "ATT" = "I", "ATG" = "M",
    "ACA" = "T", "ACC" = "T", "ACG" = "T", "ACT" = "T",
    "AAC" = "N", "AAT" = "N", "AAA" = "K", "AAG" = "K",
    "AGC" = "S", "AGT" = "S", "AGA" = "R", "AGG" = "R",
    "CTA" = "L", "CTC" = "L", "CTG" = "L", "CTT" = "L",
    "CCA" = "P", "CCC" = "P", "CCG" = "P", "CCT" = "P",
    "CAC" = "H", "CAT" = "H", "CAA" = "Q", "CAG" = "Q",
    "CGA" = "R", "CGC" = "R", "CGG" = "R", "CGT" = "R",
    "GTA" = "V", "GTC" = "V", "GTG" = "V", "GTT" = "V",
    "GCA" = "A", "GCC" = "A", "GCG" = "A", "GCT" = "A",
    "GAC" = "D", "GAT" = "D", "GAA" = "E", "GAG" = "E",
    "GGA" = "G", "GGC" = "G", "GGG" = "G", "GGT" = "G",
    "TCA" = "S", "TCC" = "S", "TCG" = "S", "TCT" = "S",
    "TTC" = "F", "TTT" = "F", "TTA" = "L", "TTG" = "L",
    "TAC" = "Y", "TAT" = "Y", "TAA" = "*", "TAG" = "*",
    "TGC" = "C", "TGT" = "C", "TGA" = "*", "TGG" = "W"
  )
  
  protein <- paste(sapply(seq(1, nchar(dna_sequence) - 2, by=3),
                          function(i) codon_table[[substr(dna_sequence, i, i+2)]],
                          USE.NAMES = FALSE), collapse = "")
  return(protein)
}

# 2. Logistic Growth Function
logistic_growth <- function(time, K, r, lag, exp_duration) {
  if (time < lag) {
    return(0.1)  
  } else if (time < lag + exp_duration) {
    return(0.1 * exp(r * (time - lag)))  
  } else {
    return(K / (1 + ((K - 0.1) / 0.1) * exp(-r * (time - (lag + exp_duration)))))
  }
}

simulate_growth_curve <- function(K=100, r=0.5, t_max=100) {
  lag <- sample(2:10, 1)
  exp_duration <- sample(5:15, 1)
  time <- seq(0, t_max, length.out=100)
  population <- sapply(time, function(t) logistic_growth(t, K, r, lag, exp_duration))
  return(data.frame(Time=time, Population=population))
}

# Generate 100 Growth Curves
generate_growth_curves <- function(n=100) {
  library(dplyr)
  growth_data <- bind_rows(lapply(1:n, function(i) simulate_growth_curve()), .id = "Curve_ID")
  return(growth_data)
}

# 3. Determine Time to Reach 80% of Carrying Capacity
time_to_reach_80_percent <- function(df, K=100) {
  threshold <- 0.8 * K
  reached_time <- df[df$Population >= threshold, ][1, "Time"]
  return(reached_time)
}

# 4. Hamming Distance Calculation
hamming_distance <- function(str1, str2) {
  max_len <- max(nchar(str1), nchar(str2))
  str1 <- sprintf("%-*s", max_len, str1)  
  str2 <- sprintf("%-*s", max_len, str2)
  sum(strsplit(str1, "")[[1]] != strsplit(str2, "")[[1]])
}

# Example Usage:
dna_sequence <- "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
print(translate_dna(dna_sequence))

growth_df <- simulate_growth_curve()
print(head(growth_df))

 all_curves_df <- generate_growth_curves()
print(head(all_curves_df))

print(time_to_reach_80_percent(growth_df))


print(hamming_distance("abirbenamor", "abirbnam"))

