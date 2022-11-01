# Corinne Jackson
# Last Edited: 2022-10-28
# BINF*6210 Assignment 2

# Can we use ND1 sequence data on the order Diptera to sort flies into their respective genera?
# https://en.wikipedia.org/wiki/Fly

# PART 1: Installing Packages ----

# First, we need to install and open the libraries required for this project

install.packages("rentrez")
library(rentrez)
install.packages("seqinr")
library(seqinr)
install.packages("BiocManager")
BiocManager::install("Biostrings")
library(Biostrings)
install.packages("tidyverse")
library(tidyverse)
install.packages("caTools")
library(caTools)
install.packages("caret")
library(caret)
install.packages("randomForest")
library(randomForest)
install.packages("RColorBrewer")
library(RColorBrewer)

# PART 2: Fetching Records ----

# ND1 is located in the mitochondrial genome from base pair 3,307 to 4,262 (length = 955bp)
# I will filter for sequences that are at least half the length of the ND1 sequence, and allow for sequences that may
# be slightly longer than the gene sequence
diptera_ND1 <- entrez_search(db = "nuccore", term = "Diptera[ORGN] AND ND1[GENE] AND 500:1200[SLEN]", retmax = 200)
maxHits <- diptera_ND1$count
# Set our retmax argument to maxHits.
diptera_ND1 <- entrez_search(db = "nuccore", term = "Diptera[ORGN] AND ND1[GENE] AND 500:1200[SLEN]", retmax = maxHits, use_history=TRUE)

# We need to get the summaries for the 1041 hits in batches because entrez_summary will throw a 414 error otherwise
chunk_size <- 200 # We will do chunks of 200 entries at a time
ids_chunked <- split(diptera_ND1$ids, ceiling(seq_along(diptera_ND1$ids)/chunk_size)) # This splits the entries into chunks
diptera_ids <- vector(mode="character", 0) # Initialize the vector that will contain the IDs after the batching is done

# Loop through each chunk and grab the ID of each entry in each chunk. We store the ID of each entry then append
# it to the vector variable we previously initialized
i= 1
while(i <= length(ids_chunked)){
  diptera_ND1_temp <- entrez_summary(db = "nuccore", id = ids_chunked[[i]])
  ids_temp <- extract_from_esummary(diptera_ND1_temp, "uid")
  diptera_ids <- c(diptera_ids, ids_temp)
  i = i + 1
}

# Now that we have extracted all the IDs, we will have to get the sequence data associated with those in chunks
ids_chunked <- split(diptera_ids, ceiling(seq_along(diptera_ids)/chunk_size)) # This splits the entries into chunks of 200
diptera_seq <- vector(mode="character", 0) # Initialize the vector that will contain the sequences after the batching is done

# Loop through each chunk of IDs and grab the corresponding sequence. We store the sequence of each ID then append
# it to the vector variable we previously initialized
i= 1
while(i <= length(ids_chunked)){
  temp_seq <- entrez_fetch(db = "nuccore", id = ids_chunked[[i]], rettype = "fasta")
  diptera_seq <- c(diptera_seq, temp_seq)
  i = i + 1
}

# Save the file of sequences
write(diptera_seq, "diptera_seq.fasta", sep = "\n")

# PART 3: Quality Checks and Filtering ----

# Convert the sequences to a string set then make a dataframe
stringSet <- readDNAStringSet("diptera_seq.fasta")
df_dipteraND1 <- data.frame(ND1_Title = names(stringSet), ND1_Sequence = paste(stringSet))

# Here we will count the length of each sequence in the dataframe by looping through them and adding them
# to a numeric vector
seq_length <- vector(mode="numeric", 0)

i = 1
while(i <= length(df_dipteraND1$ND1_Sequence)){
  seq_length <- c(seq_length, nchar(df_dipteraND1$ND1_Sequence[i]))
  i = i + 1
}

# Let's see the distribution of sequence lengths; we should not have any below 500 or any above 1200nt
hist(seq_length,
        main = "Distribution of Sequence Lengths",
        breaks = 10,
        xlab = "Sequence Length (nt)",
        ylab = "Number of Sequences",
        ylim = c(0,600),
        col = brewer.pal(n = 5, name = "YlGnBu"),
        las=2)

# Based on the distribution, we can see that the majority of the sequences used in this analysis are
# incomplete cds's and are actually around 50% of the total length of ND1

# Next we will work on extracting the genus under Diptera that each sequence actually belongs to
# Make a new column to keep the sequence title but grab the genus name
df_dipteraND1$Genus_Name <- word(df_dipteraND1$ND1_Title, 2L)
# Rearrange the columns.
df_dipteraND1 <- df_dipteraND1[, c("ND1_Title", "Genus_Name", "ND1_Sequence")]

# Now we can take a look at how many occurrences of each genus there are
genus_counts <- df_dipteraND1 %>% 
  group_by(Genus_Name) %>%
  summarise(count=n())

view(genus_counts)

# From here we can see that the majority of the genera are only represented a few times, which is not enough for us 
# to include in our classification, so we need to remove these records from our dataframe.

# We will keep 4 genera that occur the most: Anastrepha, Bactrocera, Drosophila, and Melaloncha
genus_filtered_df <-df_dipteraND1[df_dipteraND1$Genus_Name == "Anastrepha" | df_dipteraND1$Genus_Name == "Bactrocera" | df_dipteraND1$Genus_Name == "Drosophila" | df_dipteraND1$Genus_Name == "Melaloncha", ]

# Notably, there severe class imbalance between the genera which we will need to address before we train the model. 
# This will come later after all the sequences have been filtered and features have been added.

# Now we will remove sequence gaps and Ns on the ends and within the sequences.
# Filtering out sequences with a high degree (>5%) of internal N content
genus_filtered_df2 <- genus_filtered_df %>%
  mutate(ND1_Sequence2 = str_remove(ND1_Sequence, "^[-N]+")) %>%
  mutate(ND1_Sequence2 = str_remove(ND1_Sequence2, "[-N]+$")) %>%
  mutate(ND1_Sequence2= str_remove_all(ND1_Sequence2, "-+")) %>%
  filter(str_count(ND1_Sequence2, "N") <= (0.05 * str_count(ND1_Sequence)))

# Because the majority of our sequences are approximately half the length of the complete ND1 gene, we will filter
# our data so that we only keep sequences of similar length:

# We are assigning a variable name to hold the information about the first quartile sequence length and third
#quartile length so that we can exclude the data that falls above and below these points, respectively.
q1 <- quantile(nchar(genus_filtered_df2$ND1_Sequence2), probs = 0.25, na.rm = TRUE)
q1

q3 <- quantile(nchar(genus_filtered_df2$ND1_Sequence2), probs = 0.75, na.rm = TRUE)
q3

# Now we can use our calculated q1 and q3 lengths. We are keeping sequences with length greater than or equal 
# to the first quartile length value and also sequences with length equal to or less than the third quartile value to
# constrain length variability. This is important for model training as outliers can negatively affect model
# performance and increase error variance.
genus_filtered_df3 <- genus_filtered_df2  %>%
  filter((str_count(ND1_Sequence2) >= q1 & str_count(ND1_Sequence2) <= q3))

# Next we will convert our genus-filtered data to a data frame with the filtered sequence data
dfSeq <- as.data.frame(genus_filtered_df3[,c("Genus_Name", "ND1_Sequence2")])

# Let's clear up some space now that we've done all that processing...
rm(seq_length, genus_counts2, genus_filtered_df3, diptera_ND1, diptera_seq, diptera_ids, maxHits, temp_seq, chunk_size, ids_chunked, diptera_seq, stringSet, df_dipteraND1, genus_counts, genus_filtered_df, genus_filtered_df2, diptera_ND1_temp, i, ids_temp, q1, q3)

# PART 4: Adding Features ----

# Convert the nucleotides to a DNAStringSet (class) so that we can use functions from the Biostrings package
dfSeq$ND1_Sequence2 <- DNAStringSet(dfSeq$ND1_Sequence2)

# Calculating the nucleotide frequencies and appending onto the dataframe
dfSeq <- cbind(dfSeq, as.data.frame(letterFrequency(dfSeq$ND1_Sequence2, letters = c("A", "C","G", "T"))))

# I am not going to calculate the proportion of each nucleotide as a feature because this information
# would become redundant as I've already included the counts of each nucleotide. Redundant features
# can slow down the time it takes to train the model and ultimately affect model performance.

# Adding dinucleotide frequency (k-mers of length 2) to account for variability in sequence lengths.
dfSeq <- cbind(dfSeq, as.data.frame(dinucleotideFrequency(dfSeq$ND1_Sequence2, as.prob = TRUE)))

# Adding trinucleotide frequency (k-mers of length 3)
dfSeq <- cbind(dfSeq, as.data.frame(trinucleotideFrequency(dfSeq$ND1_Sequence2, as.prob = TRUE)))
df <- dfSeq[-c(2)]
df <- data.frame(df)

# PART 5: Building the Model ----

# To begin, I will train a model using all of the features, then I will remove features to see
# if it can maintain performance.

# Because the up_sample function requires the classes to be factorized, I'll factorize them right
# off the bat so I don't have to do it in every subsequent use of the function
df$Genus_Name <- as.factor(df$Genus_Name)

# Set seed so that this analysis can be reproduced
set.seed(6210)

# I will partition the dataframe into a training set and a test set. I will use an 80:20 train:test ratio...
split <- sample.split(df, SplitRatio = 0.80)
train <- subset(df, split == "TRUE")
test <- subset(df, split == "FALSE")

# As previously mentioned, we have severe class imbalance present in our data set. Class imbalance will
# affect the performance of our model, so we need to correct it before training by oversampling.
# I am choosing oversampling as opposed to undersampling so there is no data loss.
# I will only perform oversampling of the training data to avoid overfitting of the model.
up_train <- upSample(x = train[-1], y = train$Genus_Name)

# Fitting Random Forest to the train dataset
classifier_RF <- randomForest(x = up_train[-85],
                             y = up_train$Class,
                             ntree = 50)
classifier_RF

# Let's test our model on the testing dataset and see how well it makes its predictions on unseen data
y_pred <- predict(classifier_RF, newdata = test[-1])

# We can look at the confusion matrix to see how well it performed
conf_matrix <- confusionMatrix(data = y_pred, reference = test[, 1])
conf_matrix 

# Perfect performance

# PART 6: Building the Second Model ----

# Since we still have perfect performance with unseen data, lets see if it performs as well if we only train with 
# the nucleotide counts and 2-mer features

# Removing the 3-mers from the training and test sets
test2 <- test[,1:21]
up_train2 <- up_train[,1:20]
up_train2$Class <- up_train$Class

# We will build another random forest model using this updated data and the corresponding response data
classifier_RF2 <- randomForest(x = up_train2[-21],
                             y = up_train2$Class,
                             ntree = 50)

# Let's evaluate our model using the testing dataset and see how well it makes its predictions on unseen data
y_pred2 <- predict(classifier_RF2, newdata = test2[-1])

# We can look at the confusion matrix to see how well it performed
conf_matrix2 <- confusionMatrix(data = y_pred2, reference = test2[, 1])
conf_matrix2

# Perfect performance again

# PART 7: Building the Third Model ----

# We still have perfect performance with unseen data, so lets see if it performs as well if we only provide it with 
# the nucleotide counts only

# Removing the 2-mers from the training and test sets
test3 <- test2[,1:5]
up_train3 <- up_train[,1:4]
up_train3$Class <- up_train$Class

# We will build another random forest model using this updated data and the corresponding response data
classifier_RF3 <- randomForest(x = up_train3[-5],
                               y = up_train3$Class,
                               ntree = 50)

# Let's test our model on the testing dataset and see how well it makes its predictions on unseen data
y_pred3 <- predict(classifier_RF3, newdata = test3[-1])

# We can look at the confusion matrix to see how well it performed
conf_matrix3 <- confusionMatrix(data = y_pred3, reference = test3[, 1])
conf_matrix3

# Again we have perfect performance, indicating that only nucleotide counts are required for a random forest
# model to accurately predict the genus that a sequence belongs to, at least between the four genera included
# in this analysis.

# Given that only the nucleotide counts are required, let's see which nucleotide(s) are most important for
# the model's performance:

plot(sort(classifier_RF3$importance),
     main = "Importance of Nucleotide Counts for Genus Discrimination",
     xlab="Nucleotide", 
     ylab="Mean Decrease in Gini Coefficient",
     pch = 16,
     col = brewer.pal(n = 5, name = "Dark2"),
     cex = 2,
     xaxt= "n")
axis(1, at = c(1.0, 2.0, 3.0, 4.0), labels = c("T", "C", "G", "A"), tck = 0)

# A and G are much more important than T and C; more research into any biological significance of these nucleotides
# in the ND1 sequence would be required to determine if there is any reason for this.
