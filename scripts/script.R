#załadowanie bibliotek
library(tm)
library(hunspell)
library(stringr)
library(proxy)
library(lsa)
library(corrplot)
library(dendextend)
library(flexclust)
library(fossil)
library(topicmodels)
library(wordcloud)
library(slowraker)



#zmiana kataogu roboczego
workDir <- "C:\\Users\\USER\\Documents\\GitHub\\Projekty08\\Projekt08"
setwd(workDir)

#definicja katalog�w funkcjonalnych
inputDir <- ".\\data"
outputDir <- ".\\results"
scriptsDir <- ".\\scripts"
workspacesDir <- ".\\workspaces"
dir.create(outputDir, showWarnings = FALSE)
dir.create(workspacesDir, showWarnings = FALSE)

#utworzenie korpusu dokment�w
corpusDir <- paste(
  inputDir, 
  "Teksty",
  sep = "\\"
)
corpus <- VCorpus(
  DirSource(
    corpusDir,
    pattern = "*.txt",
    encoding = "UTF-8"
  ),
  readerControl = list(
    language = "pl_PL"
  )
)
#usuni�cie z tekst�w podzia�u na akapity
pasteParagraphs <- content_transformer(function(x,char) paste(x, collapse = char))
corpus <- tm_map(corpus, pasteParagraphs, " ")

#wst�pne przetwarzanie
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, content_transformer(tolower))
stoplistFile <- paste(
  inputDir, 
  "stopwords_pl.txt",
  sep = "\\"
)
stoplist <- readLines(stoplistFile, encoding = "UTF-8")
corpus <- tm_map(corpus, removeWords, stoplist)
corpus <- tm_map(corpus, stripWhitespace)

#usuni�cie em dash i 3/4
removeChar <- content_transformer(function(x,pattern) gsub(pattern, "", x))
corpus <- tm_map(corpus, removeChar, intToUtf8(8722))
corpus <- tm_map(corpus, removeChar, intToUtf8(190))

#lematyzacja
polish <- dictionary(lang="pl_PL")
lemmatize <- function(text) {
  simpleText <- str_trim(as.character(text))
  parsedText <- strsplit(simpleText, split = " ")
  newTextVec <- hunspell_stem(parsedText[[1]], dict = polish)
  for (i in 1:length(newTextVec)) {
    if (length(newTextVec[[i]]) == 0) newTextVec[i] <- parsedText[[1]][i]
    if (length(newTextVec[[i]]) > 1) newTextVec[i] <- newTextVec[[i]][1]
  }
  newText <- paste(newTextVec, collapse = " ")
  return(newText)
}
corpus <- tm_map(corpus, content_transformer(lemmatize))

#usuni�cie rozszerze� z nazw plik�w
cutExtensions <- function(document){
  meta(document, "id") <- gsub(pattern = "\\.txt$", replacement = "", meta(document, "id"))
  return(document)
}
corpus <- tm_map(corpus, cutExtensions)

#eksport zawarto�ci korpusu do plik�w tekstowych
preprocessedDir <- paste(
  outputDir, 
  "Teksty-przetworzone",
  sep = "\\"
)
dir.create(preprocessedDir)
writeCorpus(corpus, path = preprocessedDir)

#wy�wietlenie zawarto�ci pojeedynczego dokumentu
writeLines(as.character(corpus[[1]]))
writeLines(corpus[[1]]$content)

#utworzenie macierzy cz?stosci
tdmTfAll <- TermDocumentMatrix(corpus)
dtmTfAll <- DocumentTermMatrix(corpus)
tdmBinAll <- TermDocumentMatrix(
  corpus, 
  control = list(
    weighting = weightBin
  )
)
tdmTfidfAll <- TermDocumentMatrix(
  corpus, 
  control = list(
    weighting = weightTfIdf
  )
)
tdmTfBounds <- TermDocumentMatrix(
  corpus, 
  control = list(
    bounds = list(
      global = c(2,16)
    )
  )
)
tdmTfidfBounds <- TermDocumentMatrix(
  corpus, 
  control = list(
    weighting = weightTfIdf,
    bounds = list(
      global = c(2,16)
    )
  )
)
dtmTfidfBounds <- DocumentTermMatrix(
  corpus, 
  control = list(
    weighting = weightTfIdf,
    bounds = list(
      global = c(2,16)
    )
  )
)
dtmTfBounds <- DocumentTermMatrix(
  corpus, 
  control = list(
    bounds = list(
      global = c(2,16)
    )
  )
)

#konwersja macirzy rzadkich do macierzy klasycznch
tdmTfAllMatrix <- as.matrix(tdmTfAll)
dtmTfAllMatrix <- as.matrix(dtmTfAll)
tdmBinAllMatrix <- as.matrix(tdmBinAll)
tdmTfidfAllMatrix <- as.matrix(tdmTfidfAll)
tdmTfBoundsMatrix <- as.matrix(tdmTfBounds)
tdmTfidfBoundsMatrix <- as.matrix(tdmTfidfBounds)
dtmTfidfBoundsMatrix <- as.matrix(dtmTfidfBounds)
dtmTfBoundsMAtrix <- as.matrix(dtmTfBounds)

#eksport macirzy cz?sto?ci do pliku .csv
matrixFile <- paste(
  outputDir, 
  "tdmTfidfBounds.csv",
  sep = "\\"
)
write.table(
  tdmTfidfBoundsMatrix,
  file = matrixFile,
  sep = ";",
  dec = ",",
  col.names = NA
)
##################################### mds.r

#skalowanie wielowymiarowe (MDS)
distCos <- dist(dtmTfidfBoundsMatrix, method = "cosine")
distCosMatrix <- as.matrix(distCos)
mds <- cmdscale(distCos, eig = TRUE, k=2)

#rysowanie wykresu w oknie aplikacji
legend <- paste(paste("d", 1:19, sep = ""), rownames(distCosMatrix), sep = "<-")
x <- mds$points[,1]
y <- mds$points[,2]
plot(
  x,
  y,
  col = "blue",
  xlab = "Synthetic variable 1", 
  ylab = "Syntehtic variable 2",
  main = "Multidimensional Scalling"
)
text(
  x,
  y,
  labels = paste("d", 1:19, sep = ""),
  col = "blue",
  pos = 4
)
legend("bottom", legend, cex = 0.6, text.col = "blue")

#eksport wykresu do pliku .png
plotFile <- paste(
  outputDir, 
  "mds.png",
  sep = "\\"
)
png(file = plotFile)
plot(
  x,
  y,
  xlab = "Synthetic variable 1", 
  ylab = "Syntehtic variable 2",
  main = "Multidimensional Scalling",
  col = "blue",
  xlim = c(-0.5,0.5)
)
text(
  x,
  y,
  labels = paste("d", 1:19, sep = ""),
  col = "blue"
)
legend("bottom", legend, cex = 0.6, text.col = "blue")
dev.off()

####################### pca.r

#analiza głównych składowych
pca <- prcomp(dtmTfidfBounds)

#wykres dokumentów w przestrzeni dwuwymiatowej
legend <- paste(paste("d", 1:19, sep = ""), rownames(dtmTfidfBounds), sep = "<-")
options(scipen = 5)
x <- pca$x[,1]
y <- pca$x[,2]
plot(
  x, 
  y, 
  pch = 1, 
  col = "blue"
)
text(
  x, 
  y, 
  paste("d", 1:19, sep = ""), 
  col = "blue",
  pos = 4
)
legend(0.01, 0.05, legend, text.font = 3, cex = 0.5, text.col = "blue")

#eksport wykresu do pliku .png
plotFile <- paste(
  outputDir, 
  "pca.png",
  sep = "\\"
)
png(file = plotFile)
options(scipen = 5)
plot(
  x, 
  y, 
  #xlim = c(,),
  #ylim = c(,),
  pch = 1, 
  col = "blue"
)
text(
  x, 
  y, 
  paste("d", 1:19, sep = ""), 
  col = "blue",
  pos = 3
)
legend("bottomright", legend, text.font = 3, cex = 0.5, text.col = "blue")
dev.off()

####################### lsa.r

#analiza ukrytych wymiarów semantycznych (dekompozycja wg. wartości osobliwych)
lsa <- lsa(tdmTfidfBoundsMatrix)
lsa$tk #odpowiednik macierzy U, współrzędne wyrazów
lsa$dk #odpowiednik macierzy V, współrzędne dokumentów
lsa$sk #odpowiednik macierzy D, znaczenie składowych

#przygotowanie danych do wykresu
coordTerms <- lsa$tk%*%diag(lsa$sk)
coorDocs <- lsa$dk%*%diag(lsa$sk)
terms <- c("przyjęcie", "susan", "samochód", "andrew", "chłopak", "demon", "gabrielle", "katherine")
termsImportance <- diag(lsa$tk%*%diag(lsa$sk)%*%t(diag(lsa$sk))%*%t(lsa$tk))
importantTerms <- names(tail(sort(termsImportance),8))
coordTerms <- coordTerms[terms,]
#coordTerms <- coordTerms[importantTerms,]
legend <- paste(paste("d", 1:19, sep = ""), rownames(coorDocs), sep = "<-")
x1 <- coorDocs[,1]
y1 <- coorDocs[,2]
x2 <- coordTerms[,1]
y2 <- coordTerms[,2]

#wykres dokumentów i wybranych słów w przestrzeni dwuwymiatowej
options(scipen = 5)
plot(
  x1, 
  y1, 
  xlim = c(-0.2,0.05),
  #ylim = c(,),
  pch = 1, 
  col = "orange"
)
points(
  x2, 
  y2, 
  pch = 2, 
  col = "brown"
)
text(
  x1, 
  y1, 
  paste("d", 1:19, sep = ""), 
  col = "orange",
  pos = 4
)
text(
  x2, 
  y2, 
  rownames(coordTerms), 
  col = "brown",
  pos = 4
)
legend("bottomleft", legend, cex = 0.7, text.col = "orange")

#eksport wykresu do pliku .png
plotFile <- paste(
  outputDir, 
  "lsa.png",
  sep = "\\"
)
png(file = plotFile)
options(scipen = 5)
plot(
  x1, 
  y1, 
  xlim = c(-0.2,0.05),
  #ylim = c(,),
  pch = 1, 
  col = "orange"
)
points(
  x2, 
  y2, 
  pch = 2, 
  col = "brown"
)
text(
  x1, 
  y1, 
  paste("d", 1:19, sep = ""), 
  col = "orange",
  pos = 4
)
text(
  x2, 
  y2, 
  rownames(coordTerms), 
  col = "brown",
  pos = 4
)
legend("bottomleft", legend, cex = 0.5, text.col = "orange")
dev.off()

############# clustering.r

#analiza skupie?
#hierarchiczna
#parametry metody:
#1. macierz cz?sto?ci
#a. waga (weighting)
#b. zakres amiennych (bounds)
#2. miara odleg?o?ci (euclidean, jaccard, cosine)
#3. odleg?o?? pomi?dzy skupieniami (single, complete, ward.D2)

#eksperyment 1
dist1 <- dist(dtmTfAllMatrix, method = "euclidean")
hclust1 <- hclust(dist1, method = "ward.D2")
plot(hclust1)
barplot(hclust1$height, names.arg = 19:1)

#eksperyment 2
dist2 <- dist(dtmTfidfBoundsMatrix, method = "cosine")
hclust2 <- hclust(dist2, method = "ward.D2")
plot(hclust2)
barplot(hclust2$height, names.arg = 19:1)

#eksperyment 3
dist3 <- dist(coorDocs, method = "cosine")
hclust3 <- hclust(dist3, method = "ward.D2")
plot(hclust3)
barplot(hclust3$height, names.arg = 19:1)

#podzia? obiekt?w na skupienia przy zadanej liczbie klas
#eksperyment 1
clusters1 <- cutree(hclust1, k = 4)
clustersMatrix1 <- matrix(0, 20, 4)
rownames(clustersMatrix1) <- names(clusters1)
for (i in 1:20){
  clustersMatrix1[i,clusters1[i]] <- 1
}
corrplot(clustersMatrix1)

#eksperyment 2
clusters2 <- cutree(hclust2, k = 3)
clustersMatrix2 <- matrix(0, 20, 3)
rownames(clustersMatrix2) <- names(clusters2)
for (i in 1:20){
  clustersMatrix2[i,clusters2[i]] <- 1
}
corrplot(clustersMatrix2)

#eksperyment 3
clusters3 <- cutree(hclust3, k = 3)
clustersMatrix3 <- matrix(0, 20, 3)
rownames(clustersMatrix3) <- names(clusters3)
for (i in 1:20){
  clustersMatrix3[i,clusters3[i]] <- 1
}
corrplot(clustersMatrix3)

#por?wnanie wynik?w eksperyment?w
dendrogram1 <- as.dendrogram(hclust1)
dendrogram2 <- as.dendrogram(hclust2)
dendrogram3 <- as.dendrogram(hclust3)

Bk_plot(
  dendrogram1, 
  dendrogram2, 
  add_E = FALSE,
  rejection_line_asymptotic = FALSE,
  main = "Indeks Fawlksa - Mallowsa",
  ylab = "Indeks Fawlksa - Mallowsa"
)

Bk_plot(
  dendrogram1, 
  dendrogram3, 
  add_E = FALSE,
  rejection_line_asymptotic = FALSE,
  main = "Indeks Fawlksa - Mallowsa",
  ylab = "Indeks Fawlksa - Mallowsa"
)

Bk_plot(
  dendrogram2, 
  dendrogram3, 
  add_E = FALSE,
  rejection_line_asymptotic = FALSE,
  main = "Indeks Fawlksa - Mallowsa",
  ylab = "Indeks Fawlksa - Mallowsa"
)

#niehierarchiczna (k-?rednich)
#parametry metody:
#1. macierz cz?sto?ci
#a. waga (weighting)
#b. zakres amiennych (bounds)
#2. zak?adana liczba skupie?

#eksperyment 4
kmeans1 <- kmeans(dtmTfidfBounds, centers = 3)
clustersMatrix4 <- matrix(0, 20, 3)
rownames(clustersMatrix4) <- names(kmeans1$cluster)
for (i in 1:20){
  clustersMatrix4[i,kmeans1$cluster[i]] <- 1
}
corrplot(clustersMatrix4)

#wsp??czynnik zbie?no?ci klasyfikacji przy zadanej liczbie klas
randEx1Ex4 <- rand.index(clusters1, kmeans1$cluster)
randEx2Ex4 <- rand.index(clusters2, kmeans1$cluster)
randEx3Ex4 <- rand.index(clusters3, kmeans1$cluster)
randEx2Ex3 <- rand.index(clusters2, clusters3)

######################### lda.r

#analiza ukrytej alokacji Dirichlet'a
nWords <- ncol(dtmTfAll)
nTopics <- 3
lda <- LDA(
  dtmTfAll, 
  k = nTopics, 
  method = "Gibbs", 
  control = list(
    burnin = 2000, 
    thin = 100, 
    iter = 3000
  )
)
perplaxity <- perplexity(lda, dtmTfAll)
results <- posterior(lda)

#prezentacja tematów
par(mai = c(1, 2, 1, 1))
topic1 <- head(sort(results$terms[1,], decreasing = TRUE), 20)
barplot(
  rev(topic1),
  horiz = TRUE,
  las = 1, 
  main = "Temat 1",
  xlab = "Prawdopodobieństwo",
  col = 'orange'
)
topic2 <- head(sort(results$terms[2,], decreasing = TRUE), 20)
barplot(
  rev(topic2),
  horiz = TRUE,
  las = 1, 
  main = "Temat 2",
  xlab = "Prawdopodobieństwo",
  col = 'turquoise'
)
topic3 <- head(sort(results$terms[3,], decreasing = TRUE), 20)
barplot(
  rev(topic3),
  horiz = TRUE,
  las = 1, 
  main = "Temat 3",
  xlab = "Prawdopodobieństwo",
  col = 'violet'
)

#prezentacja dokumentów
document1 <- results$topics[1,]
barplot(
  rev(document1),
  horiz = TRUE,
  las = 1, 
  main = rownames(results$topics)[1],
  xlab = "Prawdopodobieństwo",
  col = 'orange'
)
document9 <- results$topics[9,]
barplot(
  rev(document9),
  horiz = TRUE,
  las = 1, 
  main = rownames(results$topics)[9],
  xlab = "Prawdopodobieństwo",
  col = 'turquoise'
)
document16 <- results$topics[16,]
barplot(
  rev(document16),
  horiz = TRUE,
  las = 1, 
  main = rownames(results$topics)[16],
  xlab = "Prawdopodobieństwo",
  col = 'violet'
)

#udział tematów w słowach
words1 <- c("susan", "chłopak", "demon")
round(results$terms[,words1],2)

words2 <- c("danielle", "kobieta", "inkwizytorka")
round(results$terms[,words2],2)

#podział dokumentów na skupienia na podstawie dominujących tematyk

######################### keywords.r

#dla pierwszego dokumentu
##wagi tf jako miara wa?no?ci s??w
keywordsTf1 <- head(sort(dtmTfAllMatrix[3,], decreasing = TRUE))
keywordsTf1
##wagi tfidf jako miara wa?no?ci s??w
keywordsTfidf1 <- head(sort(dtmTfidfBoundsMatrix[3,], decreasing = TRUE))
keywordsTfidf1
##lda jako miara wa?no?ci s??w
importance1 <- c(results$topics[3,]%*%results$terms)
names(importance1) <- colnames(results$terms)
keywordsLda1 <- head(sort(importance1, decreasing = TRUE))
keywordsLda1
##chmura tag?w
par(mai = c(0,0,0,0))
wordcloud(corpus[3], max.words = 200,colors=brewer.pal(8, "PuOr"))
##algorytm RAKE
text1 <- as.character(corpus[3])
rake1 <- slowrake(txt = text1, stem = FALSE, stop_pos = NULL)
print(rake1[[3]])

