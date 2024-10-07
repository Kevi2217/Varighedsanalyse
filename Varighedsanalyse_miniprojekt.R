library(ggplot2)
library(dplyr)
library(tidyr)

# Load necessary libraries
library(ggplot2)

# Loop through each column in the dataset
for (i in names(melanoma30)) {
  # Check if the column is numeric before plotting
  if (is.numeric(melanoma30[[i]])) {
    
    # Determine the bin width dynamically based on the range of the variable
    bin_width <- (max(melanoma30[[i]], na.rm = TRUE) - min(melanoma30[[i]], na.rm = TRUE)) / 30
    
    # Histogram
    p_hist <- ggplot(melanoma30, aes_string(x = i)) +
      geom_histogram(binwidth = bin_width) +
      ggtitle(paste("Histogram of", i)) +
      xlab(i) +
      ylab("Frequency") +
      theme_minimal()
    print(p_hist)
    
    # Boxplot
    p_box <- ggplot(melanoma30, aes_string(y = i)) +
      geom_boxplot() +
      ggtitle(paste("Boxplot of", i)) +
      ylab(i) +
      theme_minimal()
    print(p_box)
    
  } else if (is.factor(melanoma30[[i]]) || is.character(melanoma30[[i]])) {
    # For categorical variables, we create bar plots
    p_bar <- ggplot(melanoma30, aes_string(x = i)) +
      geom_bar() +
      ggtitle(paste("Bar Plot of", i)) +
      xlab(i) +
      ylab("Count") +
      theme_minimal()
    print(p_bar)
    
  } else {
    message(paste("Skipping", i, ": Not a numeric or categorical variable."))
  }
}




  