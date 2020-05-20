#załadowanie bibliotek
library(tm)
library(hunspell)
library(stringr)

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
  "Teksty",
  sep = "\\"
)
dir.create(preprocessedDir)
writeCorpus(corpus, path = preprocessedDir)

#wy�wietlenie zawarto�ci pojeedynczego dokumentu
writeLines(as.character(corpus[[1]]))
writeLines(corpus[[1]]$content)

