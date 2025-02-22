# Darko Gavrilovic 
# BIOL 448 RP R Script
# Fibropapillomatosis in Green Sea Turtles (Chelonia Mydas)
# Initial Data Set Import/Modification


################################## Set Working Directory #############################

setwd("C:/users/gavri/OneDrive/Desktop/BIOL 448 RP/FP R Scripts")


############################### Working Packages/Libraries ###########################

library(tidyverse)
library(ggplot2)
library(dplyr)
library(vegan)
library(viridis)
library(data.table)
library(readr)

############################### Import Data Set #####################################

turtle_data <- read_csv("all-activity-data.csv")


######################### Filter for Green Turtles Only ############################

cm_turtles <- turtle_data %>% filter(species == "CM")


#################### Create a Column for Julian Time in Data Set ####################


# Convert Date to "Date" format
cm_turtles$date_yyyy.mm.dd <- as.Date(cm_turtles$date_yyyy.mm.dd, format="%Y-%m-%d")

str(cm_turtles$date_yyyy.mm.dd)
# Date[1:32982], format: "2006-03-02" "2006-03-05" "2006-03-07" "2006-03-07" "2006-03-12" "2006-03-15" ...
# date data is already in "Date" form, no modification required for Julian function

# Get julian dates relative to Jan 1st, 1970...
cm_turtles$julian_time <- as.numeric(cm_turtles$date_yyyy.mm.dd)


# Here is the conversion for normal date....
# as.Date(xxxxx, origin = "1970-01-01")

cm_turtles <- cm_turtles %>%
  mutate(julian_time = julian(date_yyyy.mm.dd))
# created "julian_time" column


######################### Modify FP presence column to 1-0 Binary ###################

cm_turtles <- cm_turtles %>%
  mutate(fibropapillomas_binary = ifelse(is.na(fibropapillomas), 0, as.integer(fibropapillomas)))


# This creates an additional column "fibropapillomas_binary"
# where fp present = 1, absent = 0 and NA values become 0. 


######################## See which years contain FP observations ####################

years_with_FP <- cm_turtles %>%
  mutate(year = format(as.Date(date_yyyy.mm.dd), "%Y")) %>%
  filter(fibropapillomas_binary == 1) %>%
  group_by(year) %>%
  summarise(count = n()) %>%
  arrange(year)

years_with_FP
# # A tibble: 8 × 2
# year  count
# <chr> <int>
# 1 2010      3
# 2 2013      1
# 3 2014      1
# 4 2018      1
# 5 2019      3
# 6 2021     49
# 7 2022     71
# 8 2023     77

# We see that FP is first recorded in 2010, although Sarah is confident in data only between years 2013-2023


# This creates a function to sum "TRUE" Values in Fibropapillomas column

sum_true_fp <- function(cm_turtles) {
  # Ensure the 'fibropapillomas' column is logical or factor
  if (!is.logical(cm_turtles$fibropapillomas) && !is.factor(cm_turtles$fibropapillomas)) {
    stop("The 'fibropapillomas' column should be logical or factor.")
  }
  
  # Sum TRUE values in the 'fibropapillomas' column, excluding NAs
  true_count <- sum(cm_turtles$fibropapillomas == TRUE, na.rm = TRUE)
  
  return(true_count)
}


# Call the function....

true_fp_count <- sum_true_fp(cm_turtles)

true_fp_count
# [1] 206
# Returns a value of 206 turtles that have FP


# unique turtle ID to see how many turtles total 


################## Filter the Data Set to include 2013-2023 entries ##################

filtered_cm_turtles <- cm_turtles %>%
  filter(year >= 2010 & year <= 2023)  # Keep only rows within the year range 2013-2023


################## Filtered Data for Turtles ONLY WITH FP ############################

fp_only_data <- filtered_cm_turtles %>%
  filter(fibropapillomas_binary == 1)

# This also corresponds with out FP counts from earlier (not counting 2010, now there are 204 turtles recorded with FP)

# Double check filtering using only fibropapillomas column TRUE/FALSE
fp_only_data2 <- filtered_cm_turtles %>%
  filter(fibropapillomas == TRUE)
# Also shows 204 turtles.... seems to correspond.


####################### Filter for first FP Observations ONLY #######################


### ONLY TURTLES WITH FP, BUT REMOVE DUPLICATE OBSERVATIONS ###

# Filter for all turtles observed with FP over the years, plus, remove any duplicate observation
# If dont remove dublicates, pseudoreplication...

fp_only_ID <- filtered_cm_turtles %>%
  filter(fibropapillomas_binary == 1) %>%  # Step 1: Keep only turtles with fibropapillomas
  distinct(turtle_id, .keep_all = TRUE)   # Step 2: Keep only the first observation per Tag


# ORIGINALLY: I was going to remove all returning turtles, with or without FP, to 1) ensure 
# each turtle is only counted once (for stats), avoid overrepresentation of returning individuals, 
# BUT, I would lose longitudinal data which would show disease progression and which turtles develop FP over years. 

# INSTEAD: I Will only remove turtles with FP that return, to avoid pseudoreplication in the target demographic, 
# although pseudoreplication may occur in the disease-free turtles, this will allow me to track disease aqusition over time AND, 
# identify more risk factors + may be more ecologically relevent (emergence vs.occurence)


############################ Remove All Returning Turtles by ID ########################


# Make a Data set where you filter out any turtles without FP

# fp_only_ID     ## This Data Set Keeps only first observation of turtles with FP, filters out returning turtles that are observed again with FP

fp_only_first_obs <- filtered_cm_turtles %>%
  filter(fibropapillomas_binary == 1) %>%  # Keep only turtles with FP
  group_by(turtle_id) %>%  # Group by turtle ID
  arrange(year, .by_group = TRUE) %>%  # Arrange by year within each turtle ID group
  slice(1) %>%  # Keep only the first occurrence per turtle
  ungroup()  # Remove grouping to restore a normal dataframe

# Ensure NA turtle IDs with FP are kept
fp_only_first_obs <- filtered_cm_turtles %>%
  filter(fibropapillomas_binary == 1 & is.na(turtle_id)) %>%
  bind_rows(fp_only_first_obs)


# Make a data set that only contains turtles without FP

filtered_no_fp <- filtered_cm_turtles %>%
  filter(fibropapillomas_binary == 0)  # Keep only turtles without FP

# Now filter to remove any returning turtles that have been identified (those with turtle IDs)
# Leave all NA turtles (ie. new turtles)

filtered_no_fp_unique <- filtered_no_fp %>%
  distinct() %>%  # Remove exact duplicates
  group_by(turtle_id) %>%  # Group by turtle_id
  arrange(year, .by_group = TRUE) %>%  # Arrange by year within each turtle_id
  slice(1) %>%  # Keep only the first observation per turtle
  ungroup()  # Ungroup for a normal dataframe

# Ensure NA values are retained
filtered_no_fp_unique <- filtered_no_fp %>%
  filter(is.na(turtle_id)) %>%
  bind_rows(filtered_no_fp_unique)

# View the final dataset
filtered_no_fp_unique


# Now combine the two data sets above, to get a cohesive data set without returning turtles

final_combined_data <- bind_rows(filtered_no_fp_unique, fp_only_first_obs)




########################### FP Occurrence Per Year + Plot #############################


# Summarize the data to count fibropapillomas occurrences per year

fp_year_summary <- final_combined_data %>%
  filter(fibropapillomas_binary == 1) %>%  # Filter only rows where fibropapillomas_binary = 1
  group_by(year) %>%                       # Group data by year
  summarise(fp_count = n())                # Count occurrences of fibropapillomas

# Same as summary created earlier that shows how many turtles per year have FP, BUT allows for easier plotting

# PLOT #

ggplot(fp_year_summary, aes(x = year, y = fp_count)) +
  geom_line(color = "blue", size = 1) +       # Line plot
  geom_point(color = "red", size = 2) +      # Points on the line
  # geom_smooth(method="glm") +
  labs(
    title = "Occurrences of Fibropapillomas Over Time",
    x = "Year",
    y = "Number of Turtles with Fibropapillomas"
  ) +
  theme_minimal()


######################### FP Occurrence Within Years (Season) #####################

# Create a separate data set for each year

# 2019 Data 
data_2019 <- final_combined_data %>%
  filter(year == 2019)  

# 2021 Data 
data_2021 <- final_combined_data %>%
  filter(year == 2021)  

# 2022 Data 
data_2022 <- final_combined_data %>%
  filter(year == 2022)  

# 2023 Data 
data_2023 <- final_combined_data %>%
  filter(year == 2023) 


# 2019 Data Does not have enough entries...
# 2020 Has no FP Occurences...


############# Creating Plots Based on Occurence in Binned Julian Values (30 days) ############


###### Start with 2021 Data ########


# Group by julian bin, ie. all occurences in a week binned together for a smoother plot
binned_2021_julian <- data_2021 %>%
  mutate(julian_bin = cut(julian_time, breaks = seq(18694, 18991, by = 30))) %>%
  group_by(julian_bin) %>%
  summarise(fp_count = sum(fibropapillomas_binary == 1, na.rm = TRUE), .groups = "drop")

binned_2021_julian_plot <- ggplot(binned_2021_julian, aes(x = julian_bin, y = fp_count)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Fibropapillomas Occurrences Throughout 2021",
    x = "Julian Time (Months)",
    y = "Number of Turtles with Fibropapillomas"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

binned_2021_julian_plot

### Now by numerical month for visualization

library(lubridate)

binned_2021 <- data_2021 %>%
  mutate(month = month(date_yyyy.mm.dd)) %>%  # Extract numeric month (1-12)
  group_by(month) %>%
  summarise(fp_count = sum(fibropapillomas_binary == 1, na.rm = TRUE), .groups = "drop")

binned_2021_plot <- ggplot(binned_2021, aes(x = month, y = fp_count)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Fibropapillomas Occurrences Throughout 2021",
    x = "Month",  # Change label to "Month"
    y = "Number of Turtles with Fibropapillomas"
  ) +
  scale_x_continuous(breaks = 1:12) +  # Display months 1 to 12 on the x-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

binned_2021_plot


############## How about 2022 ###############3

binned_2022_julian <- data_2022 %>%
  mutate(julian_bin = cut(julian_time, breaks = seq(19049, 19336, by = 30))) %>%
  group_by(julian_bin) %>%
  summarise(fp_count = sum(fibropapillomas_binary == 1, na.rm = TRUE), .groups = "drop")


binned_2022_julian_plot <- ggplot(binned_2022_julian, aes(x = julian_bin, y = fp_count)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Fibropapillomas Occurrences Throughout 2022",
    x = "Julian Time (Months)",
    y = "Number of Turtles with Fibropapillomas"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

binned_2022_plot


### Numerical Visualization ####

binned_2022 <- data_2022 %>%
  mutate(month = month(date_yyyy.mm.dd)) %>%  # Extract numeric month (1-12)
  group_by(month) %>%
  summarise(fp_count = sum(fibropapillomas_binary == 1, na.rm = TRUE), .groups = "drop")

binned_2022_plot <- ggplot(binned_2022, aes(x = month, y = fp_count)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Fibropapillomas Occurrences Throughout 2022",
    x = "Month",  # Change label to "Month"
    y = "Number of Turtles with Fibropapillomas"
  ) +
  scale_x_continuous(breaks = 1:12) +  # Display months 1 to 12 on the x-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

binned_2022_plot



######## How about 2023 #########

binned_2023_julian <- data_2023 %>%
  mutate(julian_bin = cut(julian_time, breaks = seq(19421, 19710, by = 30))) %>%
  group_by(julian_bin) %>%
  summarise(fp_count = sum(fibropapillomas_binary == 1, na.rm = TRUE), .groups = "drop")

binned_2023_julian_plot <- ggplot(binned_2023_julian, aes(x = julian_bin, y = fp_count)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Fibropapillomas Occurrences Throughout 2023",
    x = "Julian Time (Months)",
    y = "Number of Turtles with Fibropapillomas"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

binned_2023_julian_plot


### Numerical Representation ###


binned_2023 <- data_2023 %>%
  mutate(month = month(date_yyyy.mm.dd)) %>%  # Extract numeric month (1-12)
  group_by(month) %>%
  summarise(fp_count = sum(fibropapillomas_binary == 1, na.rm = TRUE), .groups = "drop")

binned_2023_plot <- ggplot(binned_2023, aes(x = month, y = fp_count)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red", size = 2) +
  labs(
    title = "Fibropapillomas Occurrences Throughout 2023",
    x = "Month",  # Change label to "Month"
    y = "Number of Turtles with Fibropapillomas"
  ) +
  scale_x_continuous(breaks = 1:12) +  # Display months 1 to 12 on the x-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

binned_2023_plot



########################### Occurence With Julian Date (continuous) ##############################


##### This Shows the Trend alot better...######

######## OLD PLOT #######

# zero_binary_data_2021 <- data_2021 %>%
#   mutate(fibropapillomas_binary = ifelse(is.na(fibropapillomas_binary), 0, fibropapillomas_binary))


# unbinned_2021_plot <- ggplot(zero_binary_data_2021, aes(x = julian_time, y = fibropapillomas_binary)) +
#   geom_point(color = "red", size = 2) +  # Red points for data points
#   geom_smooth(method = "loess", color = "blue", size = 1) +  # Trendline (LOESS smooth)
#   labs(
#    title = "Fibropapillomas Occurrences Throughout 2021",
#    x = "Julian Time (Days)",
#    y = "Fibropapillomas Occurrence"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
# unbinned_2021_plot

###############




## 2021 ##


unbinned_2021_plot <- ggplot(data_2021, aes(x = julian_time, y = fibropapillomas_binary)) +
  geom_point(color = "red", size = 2) +  # Red points for data points
  geom_smooth(method = "gam", color = "blue", size = 1) +  # Trendline (LOESS smooth)
  labs(
    title = "Fibropapillomas Occurrences Throughout 2021",
    x = "Julian Time (Days)",
    y = "Fibropapillomas Occurrence"
  ) +
  scale_x_continuous(
    breaks = seq(min(data_2021$julian_time), max(data_2021$julian_time), by = 30),  # Adjust step size
    labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%b")  # Label with month abbreviation
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for better readability

unbinned_2021_plot




## 2022 ##

unbinned_2022_plot <- ggplot(data_2022, aes(x = julian_time, y = fibropapillomas_binary)) +
  geom_point(color = "red", size = 2) +  # Red points for data points
  geom_smooth(method = "gam", color = "blue", size = 1) +  # Trendline (LOESS smooth)
  labs(
    title = "Fibropapillomas Occurrences Throughout 2022",
    x = "Julian Time (Days)",
    y = "Fibropapillomas Occurrence"
  ) +
  scale_x_continuous(
    breaks = seq(min(data_2022$julian_time), max(data_2022$julian_time), by = 30),  # Adjust step size
    labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%b")  # Label with month abbreviation
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for better readability

unbinned_2022_plot




## 2023 ###

unbinned_2023_plot <- ggplot(data_2023, aes(x = julian_time, y = fibropapillomas_binary)) +
  geom_point(color = "red", size = 2) +  # Red points for data points
  geom_smooth(method = "gam", color = "blue", size = 1) +  # Trendline (LOESS smooth)
  labs(
    title = "Fibropapillomas Occurrences Throughout 2023",
    x = "Julian Time (Days)",
    y = "Fibropapillomas Occurrence"
  ) +
  scale_x_continuous(
    breaks = seq(min(data_2023$julian_time), max(data_2023$julian_time), by = 30),  # Adjust step size
    labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%b")  # Label with month abbreviation
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for better readability

unbinned_2023_plot




#################### Fibropapillomatosis and Turtle Size ##########################



##### How is CCL and CCW Changing Over Time #####

# Optional: Calculating average CCL per year
# avg_ccl_per_year <- turtle_data %>%
#   group_by(year) %>%
#   summarise(avg_ccl = mean(ccl, na.rm = TRUE))


# Create plot to see if temporal changes in carapace length
ggplot(final_combined_data, aes(x = as.integer(year), y = ccl_min_avg_cm)) +
  geom_point(color = "red", size = 2) +  # Red points for data points
  geom_line(color = "blue", size = 1) +  # Blue line showing temporal change
  geom_smooth(method = "loess", color = "black", size = 1, linetype = "dashed") +  # LOESS smoothing to show trend
  labs(
    title = "Temporal Change in Curved Carapace Length (CCL) Over the Years",
    x = "Year",
    y = "Average Curved Carapace Length (CCL)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for better readability

# No significant differences in CCL over years, dashed line is best fit. 

ggplot(final_combined_data, aes(x = as.integer(year), y = ccw_max_avg_cm)) +
  geom_point(color = "red", size = 2) +  # Red points for data points
  geom_line(color = "blue", size = 1) +  # Blue line showing temporal change
  geom_smooth(method = "loess", color = "black", size = 1, linetype = "dashed") +  # LOESS smoothing to show trend
  labs(
    title = "Temporal Change in Curved Carapace Width (CCW) Over the Years",
    x = "Year",
    y = "Average Curved Carapace Length (CCW)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate labels for better readability

# Same with carapace width, although, with greater variability





######### Perform Logistic Regressions/Visualizations for CCL with FP Occurence ###########

#### Remove Outliers ####

Q1ccl <- quantile(final_combined_data$ccl_min_avg_cm, 0.25, na.rm = TRUE)
Q3ccl <- quantile(final_combined_data$ccl_min_avg_cm, 0.75, na.rm = TRUE)
IQRccl <- Q3ccl - Q1ccl

# Define outliers as values outside 1.5*IQR above Q3 or below Q1
outlier_lowerccl <- Q1ccl - 1.5 * IQRccl
outlier_upperccl <- Q3ccl + 1.5 * IQRccl

# Filter out rows with extreme outliers in CCL
cleaned_dataccl <- final_combined_data %>%
  filter(!is.na(ccl_min_avg_cm)) %>%
  filter(ccl_min_avg_cm >= outlier_lowerccl & ccl_min_avg_cm <= outlier_upperccl)

# Logistic regression model to examine size (carapace length) and FP occurrence
ccl_vs_FP <- glm(fibropapillomas_binary ~ ccl_min_avg_cm, data = cleaned_dataccl, family = "binomial")

summary(ccl_vs_FP)

# p-value = 0.674


## Scatterplot for FP vs CCL ##

ggplot(cleaned_dataccl, aes(x = ccl_min_avg_cm, y = fibropapillomas_binary)) +
  geom_jitter(width = 0.2, height = 0.1, alpha = 0.5, color = "red") +  # Adding jitter for visibility
  labs(
    title = "Turtle Size vs. Fibropapillomatosis Occurrence",
    x = "Curved Carapace Length (cm)",
    y = "Fibropapillomas Occurrence (1 = Yes, 0 = No)"
  ) +
  theme_minimal()


## Box Plot for CCL ###

# Create a boxplot comparing CCL by FP status (1 = FP, 0 = no FP)
ggplot(cleaned_dataccl, aes(x = factor(fibropapillomas_binary), y = ccl_min_avg_cm)) +
  geom_boxplot(aes(fill = factor(fibropapillomas_binary)), color = "black", 
               outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Boxplot of Curved Carapace Length (CCL) by FP Occurrence",
    x = "Fibropapillomatosis (FP) Status",
    y = "Curved Carapace Length (CCL) in cm",
    fill = "FP Status"
  ) +
  scale_x_discrete(labels = c("No FP", "Has FP")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 14),
        axis.title = element_text(size = 13)) 




## Remove Outliers for CCW ##

Q1ccw <- quantile(final_combined_data$ccw_max_avg_cm, 0.25, na.rm = TRUE)
Q3ccw <- quantile(final_combined_data$ccw_max_avg_cm, 0.75, na.rm = TRUE)
IQRccw <- Q3ccw - Q1ccw

# Define outliers as values outside 1.5*IQR above Q3 or below Q1
outlier_lowerccw <- Q1ccw - 1.5 * IQRccw
outlier_upperccw <- Q3ccw + 1.5 * IQRccw

# Filter out rows with extreme outliers in CCL
cleaned_dataccw <- final_combined_data %>%
  filter(!is.na(ccw_max_avg_cm)) %>%
  filter(ccw_max_avg_cm >= outlier_lowerccw & ccw_max_avg_cm <= outlier_upperccw)


### CCW vs FP Linear Regression ###

ccw_vs_FP <- glm(fibropapillomas_binary ~ ccw_max_avg_cm + julian_time, data = cleaned_dataccw, family = "binomial")

summary(ccw_vs_FP)

# p-value = 0.792, carapace width does not significantly affect the likelihood of FP occurrence in this model.
# julian time (the time of observation) is highly statistically significant in predicting FP occurrence.


# Scatterplot for CCW vs FP

ggplot(cleaned_dataccw, aes(x = ccw_max_avg_cm, y = fibropapillomas_binary)) +
  geom_jitter(width = 0.2, height = 0.1, alpha = 0.5, color = "red") +  # Adding jitter for visibility
  labs(
    title = "Turtle Size vs. Fibropapillomatosis Occurrence",
    x = "Curved Carapace Width (cm)",
    y = "Fibropapillomas Occurrence (1 = Yes, 0 = No)"
  ) +
  theme_minimal()


# Box Plot for CCW

# Create a boxplot comparing CCW by FP status (1 = FP, 0 = no FP)
ggplot(cleaned_dataccw, aes(x = factor(fibropapillomas_binary), y = ccw_max_avg_cm)) +
  geom_boxplot(aes(fill = factor(fibropapillomas_binary)), color = "black", 
               outlier.colour = "red", outlier.shape = 16, outlier.size = 2) +
  labs(
    title = "Boxplot of Curved Carapace Width (CCW) by FP Occurrence",
    x = "Fibropapillomatosis (FP) Status",
    y = "Curved Carapace Width (CCW) in cm",
    fill = "FP Status"
  ) +
  scale_x_discrete(labels = c("No FP", "Has FP")) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        title = element_text(size = 14),
        axis.title = element_text(size = 13)) 



## Correlation tests for CCL/CCW and FP Occurrence ##

# Since FP presence is binary (1/0), I performed point-biserial correlation. 
# Like Pearson correlation but between a continuous and a binary variable

## For CCL/FP ##

cor.test(final_combined_data$ccl_min_avg_cm, final_combined_data$fibropapillomas_binary)

# No Correlation...


## For CCW/FP ##

cor.test(final_combined_data$ccw_max_avg_cm, final_combined_data$fibropapillomas_binary)

# No Correlation


##### Chi-squared test of independence ####


# Create small, medium, large categories based on ccl in data set

final_combined_data <- final_combined_data %>%
  mutate(
    ccl_category = case_when(
      ccl_min_avg_cm < 100 ~ "Small", 
      ccl_min_avg_cm >= 100 & ccl_min_avg_cm < 110 ~ "Medium", 
      ccl_min_avg_cm >= 110 ~ "Large",
      TRUE ~ NA_character_))

contingency_table <- table(final_combined_data$ccl_category, final_combined_data$fibropapillomas_binary)

chi_square_result <- chisq.test(contingency_table)
chi_square_result

# Differences not significant (p-value = 0.5991)

# testing changes in R, push to Git repo

