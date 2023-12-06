# Read the input file
setwd("/Users/gabri/Library/CloudStorage/Box-Box/CPAE_manuscript/new_analyses/project_annotate_functions")
input_file <- "GBM.inputfile.txt"
text_lines = readLines(input_file)

lines_starting_with_MF <- character(0)
lines_starting_with_K <- character(0)
lines_starting_with_slash <- character(0)

# Loop through each line and categorize them
for (line in text_lines) {
  if (startsWith(line, "MG")) {
    lines_starting_with_MF <- c(lines_starting_with_MF, line)
  } else if (startsWith(line, "K")) {
    lines_starting_with_K <- c(lines_starting_with_K, line)
  } else if (startsWith(line, "/")) {
    lines_starting_with_slash <- c(lines_starting_with_slash, line)
  }
}

# Convert the categorized lines to data frames
df_MF <- data.frame(Line_MF = lines_starting_with_MF)

lines_between_MF_and_slash <- character(0)

# Initialize a flag to track whether we are between "MF" and "/"
between_MF_and_slash <- FALSE

# Assuming you have already read the text file into the 'text_lines' vector

# Initialize an empty list to store the lines between "MF" and "/"
lines_between_MF_and_slash_list <- list()

# Initialize a variable to store the current "MF" line
current_MF_line <- NULL

# Loop through each line and extract lines between "MF" and "/"
for (line in text_lines) {
  if (startsWith(line, "MG")) {
    current_MF_line <- line
    lines_between_MF_and_slash_list[[line]] <- character(0)
  } else if (startsWith(line, "/")) {
    current_MF_line <- NULL
  } else if (!is.null(current_MF_line)) {
    lines_between_MF_and_slash_list[[current_MF_line]] <- 
      c(lines_between_MF_and_slash_list[[current_MF_line]], line)
  }
}

# Access the list elements by name (MF lines)
names(lines_between_MF_and_slash_list)
# Function to check if there are no instances of "K" followed by five numbers in a vector
check_K_followed_by_5_numbers <- function(vector) {
  no_K_followed_by_5_numbers <- sapply(vector, function(str) !grepl("\\bK[0-9]{5}\\b", str))
  return(all(no_K_followed_by_5_numbers))
}

# Apply the function to each vector in the list
results <- sapply(lines_between_MF_and_slash_list, check_K_followed_by_5_numbers)
lines_between_MF_and_slash_list = lines_between_MF_and_slash_list[-which(results == TRUE)]
lines_between_MF_and_slash_list

# Initialize an empty list to store the data frames
data_frames_list <- list()

# Loop through the list and create data frames
for (mf_line_name in names(lines_between_MF_and_slash_list)) {
  # Extract lines between "MF" and "/"
  lines_between_MF_and_slash <- lines_between_MF_and_slash_list[[mf_line_name]]
  
  # Initialize an empty vector to store the extracted strings
  extracted_strings <- character(0)
  
  # Define a regular expression pattern to match strings starting with "K" followed by five numbers
  pattern <- "\\bK[0-9]{5}\\b"
  
  # Loop through lines and extract matching strings
  for (line in lines_between_MF_and_slash) {
    matches <- regmatches(line, gregexpr(pattern, line))
    if (length(matches[[1]]) > 0) {
      extracted_strings <- c(extracted_strings, matches[[1]])
    }
  }
  
  # Create a data frame
  df <- data.frame(MF_Line_Name = mf_line_name, Extracted_Strings = extracted_strings)
  
  # Append the data frame to the list
  data_frames_list[[mf_line_name]] <- df
}

# Access the data frames in the list by name
data_frames_list[["MF_line_name"]]



# Combine the individual data frames into a single data frame
combined_df <- do.call(rbind, data_frames_list)

# Rename the columns if needed
colnames(combined_df) <- c("MF_Line_Name", "Extracted_Strings")
write.csv(combined_df,"gbm_tab.txt")
