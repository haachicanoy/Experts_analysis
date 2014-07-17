options(warn=-1)

## Text mining example

# ========================================================================= #

policy.HTML.page <- readLines("http://policy.unt.edu/policy/3-5")
length(policy.HTML.page)

policy.HTML.page[186:202]

id.1 <- 3+which(policy.HTML.page=="                    TOTAL UNIVERSITY        </div>")
id.2 <- id.1+5
text.data <- policy.HTML.page[id.1:id.2]
rm(policy.HTML.page, id.1, id.2)

text.data

td.1 <- gsub(pattern="<p>", replacement="", x=text.data,
             ignore.case=T, perl=F, fixed=F, useBytes=F)

td.2 <- gsub(pattern="</p>", replacement="", x=td.1,
             ignore.case=T, perl=F, fixed=F, useBytes=F)

text.d <- td.2; rm(text.data, td.1, td.2)

library(tm)
txt <- VectorSource(text.d); rm(text.d)
text.corpus <- Corpus(txt); rm(txt)
inspect(text.corpus)

text.corpus <- tm_map(text.corpus, tolower); inspect(text.corpus)
text.corpus <- tm_map(text.corpus, removePunctuation); inspect(text.corpus)
text.corpus <- tm_map(text.corpus, removeNumbers); inspect(text.corpus)
text.corpus <- tm_map(text.corpus, removeWords, stopwords("english")); inspect(text.corpus)

library(SnowballC)
text.corpus <- tm_map(text.corpus, stemDocument)
detach("package:SnowballC")
inspect(text.corpus)

text.corpus <- tm_map(text.corpus, stripWhitespace); inspect(text.corpus)
text.corpus <- tm_map(text.corpus, PlainTextDocument)

tdm <- TermDocumentMatrix(text.corpus)
inspect(tdm[1:20,])

#find which words were used most
findFreqTerms(x=tdm, lowfreq=8, highfreq=Inf)

#find words which associate together
findAssocs(x=tdm, term="computing", corlimit=0.6)

#only the common terms (the smaller the percentage, the fewer terms will be retained)
tdm.common.60 <- removeSparseTerms(x=tdm, sparse=0.60)
tdm.common.20 <- removeSparseTerms(x=tdm, sparse=0.20)
tdm
tdm.common.60
tdm.common.20

inspect(tdm.common.20)

library(wordcloud)
m <- as.matrix(tdm)
# calculate the frequency of words
v <- sort(rowSums(m), decreasing=TRUE)
myNames <- names(v)
d <- data.frame(word=myNames, freq=v)
wordcloud(d$word, d$freq, min.freq=3)

(freq.terms <- findFreqTerms(tdm, lowfreq=8))

library(graph)
library(Rgraphviz)
plot(tdm, term=freq.terms, corThreshold=0.1, weighting=T)

# Data mining with our data

wk_dir = "F:/CIAT/CWR_project/Experts_evaluation"
path   = paste(wk_dir,"/Downloaded_data_(2014-04-28)/clean",sep="")

maps_data   = read.csv(paste(path,"/maps_corrected.csv",sep=""),header=T)
scores_data = read.csv(paste(path,"/scores_corrected.csv",sep=""),header=T)
scores_text = read.csv(paste(path,"/sc_text_corrected.csv",sep=""),header=T)
fps_data    = read.csv(paste(wk_dir,"/allpriorities_2014-04-21.csv",sep=""),header=T)
fps_modf    = read.csv(paste(wk_dir,"/modpriorities_2014-05-15.csv",sep=""),header=T)


names(scores_text)

crops <- unique(as.character(scores_text$Crop_code))
as.character(scores_text$Comments_gap_analysis)

crop_select <- function(x, ...)
{
  require(tm)
  require(wordcloud)
  information <- scores_text[scores_text$Crop_code==paste(x),]
  gap.comments <- as.character(information$Comments_gap_analysis)
  tax.addition <- as.character(information$Additional_taxa)
  gap.notes    <- as.character(information$Notes)
  
  ### Comments about Gap Analysis
  # Remove line breaks
  gap.comm <- gsub(pattern="\n", replacement=" ", x=gap.comments,
                   ignore.case=T, perl=F, fixed=F, useBytes=F)
  comm.txt <- VectorSource(gap.comm)
  text.corpus <- Corpus(comm.txt)
  
  text.corpus <- tm_map(text.corpus, tolower)
  text.corpus <- tm_map(text.corpus, removePunctuation)
  text.corpus <- tm_map(text.corpus, removeNumbers)
  text.corpus <- tm_map(text.corpus, removeWords, stopwords("english"))
  
  library(SnowballC)
  text.corpus <- tm_map(text.corpus, stemDocument)
  detach("package:SnowballC")
  
  text.corpus <- tm_map(text.corpus, stripWhitespace)
  text.corpus <- tm_map(text.corpus, PlainTextDocument)
  tdm <- TermDocumentMatrix(text.corpus)
  
  m <- as.matrix(tdm)
  # calculate the frequency of words
  v <- sort(rowSums(m), decreasing=TRUE)
  myNames <- names(v)
  d <- data.frame(word=myNames, freq=v)
  wordcloud(d$word, d$freq, min.freq=1)
  
  ### Additional taxa suggested
  # Remove line breaks
  gap.addt <- gsub(pattern="\n", replacement=" ", x=tax.addition,
                   ignore.case=T, perl=F, fixed=F, useBytes=F)
  add.txt <- VectorSource(gap.addt)
  text.corpus <- Corpus(add.txt)
  
  text.corpus <- tm_map(text.corpus, tolower)
  text.corpus <- tm_map(text.corpus, removePunctuation)
  text.corpus <- tm_map(text.corpus, removeNumbers)
  text.corpus <- tm_map(text.corpus, removeWords, stopwords("english"))
  
  library(SnowballC)
  text.corpus <- tm_map(text.corpus, stemDocument)
  detach("package:SnowballC")
  
  text.corpus <- tm_map(text.corpus, stripWhitespace)
  text.corpus <- tm_map(text.corpus, PlainTextDocument)
  tdm <- TermDocumentMatrix(text.corpus)
  
  m <- as.matrix(tdm)
  # calculate the frequency of words
  v <- sort(rowSums(m), decreasing=TRUE)
  myNames <- names(v)
  d <- data.frame(word=myNames, freq=v)
  wordcloud(d$word, d$freq, min.freq=1)
  
  ### Notes or comments about the results
  # Remove line breaks
  gap.not <- gsub(pattern="\n", replacement=" ", x=gap.notes,
                   ignore.case=T, perl=F, fixed=F, useBytes=F)
  not.txt <- VectorSource(gap.not)
  text.corpus <- Corpus(not.txt)
  
  text.corpus <- tm_map(text.corpus, tolower)
  text.corpus <- tm_map(text.corpus, removePunctuation)
  text.corpus <- tm_map(text.corpus, removeNumbers)
  text.corpus <- tm_map(text.corpus, removeWords, stopwords("english"))
  
  library(SnowballC)
  text.corpus <- tm_map(text.corpus, stemDocument)
  detach("package:SnowballC")
  
  text.corpus <- tm_map(text.corpus, stripWhitespace)
  text.corpus <- tm_map(text.corpus, PlainTextDocument)
  tdm <- TermDocumentMatrix(text.corpus)
  
  m <- as.matrix(tdm)
  # calculate the frequency of words
  v <- sort(rowSums(m), decreasing=TRUE)
  myNames <- names(v)
  d <- data.frame(word=myNames, freq=v)
  wordcloud(d$word, d$freq, min.freq=1)
}

lapply(crops, function(x){crop_select(x)})














