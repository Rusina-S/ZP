############################
#                          #  
#  R for beginners         #
#  Part 1                  #
#                          # 
#  Getting started with R  #
#                          #
#  SMCS                    #
#                          #  
###########################

# INTRODUCTION ----

# RStudio offers many options (console location, autowrap, custom themes, ...) and shortcuts (usual ones like CTRL+C ; CTRL+S, but also CTRL + Enter)

# Always comment your code (useful for you, your future self and other users) - useful shortcut : CTRL + MAJ + C

# In addition to comments, don't hesitate to create sections

# SECTION 1 ----

# This is a section

# SECTION 2 ----

# This is another section

# Useful shortcuts for sections : ALT + O (collapse all sections) ; Alt + Maj + O (Unfold all sections)
# (There are a list of useful shortcuts)

# To execute code, you can click on Run or press CTRL + Enter

# Check working directory (usually the folder where the script is)
getwd()

# Set the working directory with complete path
# Note that you must either use "/" or "\\"
# - man jau ir īstā dir - nelikšu citu. Šo dabūju no GitHub caur šo darbību: https://happygitwithr.com/new-github-first#new-github-first
#lai varētu strā'dat caur GitHUb un kolņēt visu izdarīto GitHubā
#setwd("C:/Users/vanbenedena/Documents/Instats/Rbase part 1 training material")

# Check the updated working directory
getwd()

# list the object in R memory
ls()

# empty the working environment (the memory of your R session)
rm(list = ls())

# Simple math
1+1+1

# Store in an object
a <- 1

# Simple math (but with objects)
a + a + a

# I can now check the working environment with ls() and remove with rm(list = ls())

# install the package car
# install.packages("car")

# load the package car that is installed on the computer
library(car)

# install several packages at once - Do not run now if already installed
# install.packages(c("car","openxlsx"))

# get help on a function
help(seq)
?seq
# get help about a package
help(package = car)


####
# 1. VECTORS ----
####

#   1.1 Vector creation ------------------------------
# create a vector with c, seq and rep

# Create numeric vectors
c(2, 4, -1)
c(1, 5/6, 2^3, -.05)

# Note: by default, R display 7 decimals. If you want to modify this behaviour, you may use options(digits=2) for two decimals

# Create a string vector
c("hd", "HSK")

# What happens if I combine numeric and strings?
c("e", 2)

# Create a sequential vector
c(1, 2, 3, 4, 5, 6, 7)
1:7
1:100
5:-5

seq(from=1, to=7, by = 1)
seq(from=1, to=7, by = 0.2)
seq(from=1, to=7, by = 2.5)
seq(from=1, to=7, length.out = 10)

# repeat 1 ten times
c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1)
rep(1,10)

# What happens ?
rep(10,1)
# Order of arguments is important

# repeat "a" five times
rep("a", 5)
# repeat group1 twice and group2 three times
rep(c("group1","group2"), c(2,3))
# repeat group1 and group2 twice
rep(c("group1","group2"), 2)
rep(c("groupe1","groupe2"), c(2,2))
rep(c("groupe1","groupe2"), each=2)

# A more "complicated" example
rep(1:3, times = 3, each = 2)

#   1.2 Assignment ---------------------------
# assign a value to a name to create an OBJECT or VARIABLE
# in this way it is stored in the (session) memory (temporary)

# different ways to do it 

# <- 
x <- c(2.1, 5, -4, 1, 5) 
x
c(2.1, 5, -4, 1, 5) -> x1
# =
x2 = c(2.1, 5, -4, 1, 5)
# Careful, c(2.1, 5, -4, 1, 5) = x2 will throw an error

# assign
assign("x3", c(2.1, 5, -4, 1, 5))

# If you want to print the result of the assignment as well, use ()
(x <- c(2.1, 5, -4, 1, 5) )

# Check the memory
ls()


#   1.3 Vectors management -------------------------------------
# display/modify one or more elements of a vector
x
x[3]
x[c(1,4,3)]

# modify an element 
x[3] <- 0
x
x[3] <- -4

# display all the elements except the third one
x[-3] 
# display all the elements except the first and second one
x[-c(1,2)]

# display using the positions
x[c(1,2,4)]
# is equivalent to a vector of TRUE and FALSE
x[c(TRUE, TRUE, FALSE, TRUE, FALSE)]

#   1.4 Types -------------------------
# type of an R object (important because the results of some functions depends on the type of the used object)

class(x)
is(x)
is.numeric(x)

#   1.5 Convert vectors ---------------------------------
# convert the type of a vector

# convert into a character vector
x_char <- as.character(x)
x_char

# convert into a numeric vector
as.numeric(x_char)
as.numeric(c(x_char, "ads"))

# convert a numeric vector into a logical one
as.logical(c(x,0))

# convert a logical vector into a numeric one
as.numeric(c(TRUE, TRUE, FALSE))

#   1.6 Length of a vector ---------------------------

length(x)

# modify the length
length(x) <- 8
length(x) <- 5
length(x) <- 4


#   1.7 Computation -------------------


x <- c(2.1, 5, -4, 1)
y <- c(0,-7, 1, .25)
x; y

x + y

# 2 is recycled
x + 2

x + c(0,1)

# !!! what happens in the following line ???
x + c(1,2,3)
# when vectors have different lengths, R attempts to replicate ("recycle") the shorter one 
# when longest length is not multiple of shortest it generates a warning

x * y

min(x)
max(x)
sum(x)
prod(x)
cumsum(x)
cumprod(x)

sqrt(x)

# You can also round numeric
round(cos(x), digits = 3)



#   1.8 Comparisons --------------------

# USE OF <,>,<=,>=,==,!= OPERATORS TO COMPARE VALUES
# THIS RESULTS IN A LOGICAL VALUE (TRUE or FALSE)
# CONSTRUCT LOGICAL EXPRESSIONS USING &,| and ! 

x
x <= c(1,6,3,4)
x <= 2
# find the position (index) of the elements for which the comparison is true
which(x <= 2)
# display the elements of x less or equal than 2 
x[x<=2]
# replace the elements less than 2 by NA (missing value)
x[x<=2] <- NA

x <- c(2.1, 5, -4, 1)

# combine comparisons
# & AND, | OR, ! NOT
0 <= x & x <=3
x==1 | x > 4
!(x==1 | x > 4)

# & has priority over | in the same way as multiplication over addition in arithmetic
# USE () WHE NEEDED TO CHANGE order of PRIORITY

# Elementwise
x
x <= 2

# all are TRUE
all(x <= 2)

# any is TRUE
any(x <=2)

#   1.9 Paste -----------------------
# paste elements to create a vector

paste(c("a","a","a"), c(1,2,3))

# "" is the empty string of characters
paste(c("a","a","a"), c(1,2,3), sep="")

paste("a", 1:3, sep="")

paste(1:2,3, sep="")

# Simplified function: paste0
paste0("a", 1:3)


#   1.10 Character vectors --------------------------------

x <- paste(c("BE", "BE", "FR", "ES", "BE"), 
           1:5, 
           sep="")
x

# positions of the elements containing BE
grep("BE", x)
# positions of the elements containing E
grep("E", x)
# display the elements containing E
x[grep("E", x)]

# extract a character string between two values
country <- substr(x, start = 1, stop = 2)
country

# replace the pattern pre by suf if present
sub("pre", "suf", c("boat", "prefix"))

# split the elements based on the space
math <- c("Niels Henrik Abel", 
          "Evariste Galois", 
          "Emmy Noether")
math
strsplit(math, split = " ")

# put the letters in lower letters
tolower(math)
# put the letters in upper letters
a <- toupper(math)

#   1.11 Sorting ---------------------
# sort a vector

x <- c(2.1, 5, -4, 1, 1)
# in ascending order
sort(x)
# in decreasing order
sort(x, decreasing = TRUE)

# find the order you need to take the elements in a way that the elements are sorted
order(x)
x[order(x)]

# ranks of the elements
rank(x)

# reverse the position of the elements
rev(x)


#   1.12 Missing values--------------------


# NA = missing value
x <- c(0,4,NA,7)
# "NA" is not a missing value
c(0,4,"NA",7)

# find missing values
is.na(x)
# find the positions of the missing values
which(is.na(x))
# count the number of missing values in x
sum(is.na(x))
# know if there is at least one missing value in x
anyNA(x)

# impact of the missing values on the mean
mean(x)
# option na.rm
help(mean)
mean(x, na.rm = TRUE)
# other way: using na.omit
mean(na.omit(x))


###
# 2. FACTORS --------------------------------------------
###

# Factors are used to represent categorical variables
# Many functions or statistical analyses required your categorical variables to be stored as factors
# Under the hood, it works as a numeric with labels

# from a character vector
f1 <- factor(c("t1","t4","t1","t3","t4","t1","t3","t3"))
f1

# from a numeric vector and adding labels on levels
v <- c(1,1,0,1,0)
v2 <- factor(v, levels=c(0,1), labels = c("bad", "good"))
# be careful with the order of the levels
v3 <- factor(v, levels=c(1,0), labels = c("good", "bad"))

# levels of f1
levels(f1)
# number of levels
nlevels(f1)
# change the reference level (useful in comparing levels after an ANOVA for example)
f2 <- relevel(f1, ref = "t3")
f2
f1

# convert a factor into a numeric vector
v
v2
as.numeric(v2)

# frequency table
table(f1)



###
# 3. DATES ----------------------------------------------
###


# Get today's date (POSIXct)
today <- Sys.time()
str(today)

# Get today's date as a Date object - does not include the time
today.date <- as.Date(today)
str(today.date)

# a specific data type to store date info
# comes with a specific format (region/language specific)

#   3.1 Date manipulation ----

# !!! be careful !!! 
# error - because the second date is a string
today-"2020-01-01"

# Need either the Date type
today.date-as.Date("2020-01-01")

# Or the POSIX type
today - as.POSIXct("2020-01-01")

# Also
difftime(today, as.POSIXct("2020-01-01"), units = "secs")
difftime(today.date, as.Date("2020-01-01"), units = "secs")

# Or in hours (possible also in minutes, days, weeks, see the help page)
difftime(today, as.POSIXct("2020-01-01"), units = "hours")
difftime(today.date, as.Date("2020-01-01"), units = "hours")

# Create a vector of date
dates <- c("02/27/92", "02/27/99", "01/14/92")

# Create a vector of time
times <- c("23:03:20", "22:29:56", "01:03:30")

# Combine both
x <- paste(dates, times)

# Convert this date-time vector into a POSIX using strptime
y <- strptime(x, "%m/%d/%y %H:%M:%S")

weekdays(y, abbreviate = FALSE)
months(y, abbreviate = TRUE)
quarters(y, abbreviate = FALSE)
# year
format(y, "%Y")
format(y, "%y")



###
# 4. LISTS ------------------------------------------------
###

# Lists are a way to store multiple objects (vectors, but also other lists or dataframes)

# create
tahiti <- list(
  factor(c("Boeing", "Airbus")),
  c(15,11,14),
  c("Brussels", "Berlin", "London")
)

tahiti <- list(
  plane = factor(c("Boeing", "Airbus")),
  duration=c(15,11,14),
  departure=c("Brussels", "Berlin", "London")
)

unlist(tahiti)



# display an element
tahiti$duration
tahiti[[2]]
tahiti[["duration"]]

# Not recommended
tahiti$du

# show the first and third element of the element departure of the list tahiti
tahiti$depart[c(1,3)]
tahiti[[3]][c(1,3)]

# select multiple elements from a list to create a new list
tahiti[c(1,3)]
str(tahiti[c(1,3)])

# no list 
str(tahiti[[2]])

#   4.1 Information -----------------
# get information

attributes(tahiti)
names(tahiti)

str(tahiti)

###
# 5. MATRICES ----------------------------------------------
###

# Matrices contain two dimensions : rows and columns
# They can only contain a single type of values (numeric, string or logical)

# Create a matrix (by default it is created in column)
matrix(1:12,nrow=3,ncol=4)

# Create a matrix in rows
matrix(1:12,nrow=3,ncol=4,byrow=TRUE)

# Create a matrix from vectors by rows
rbind(1:3,c(2,3,9),rep(7,3),6:8)

# Create a matrix from vectors by column
cbind(1:3,c(2,3,9),rep(7,3),6:8)

# Check information on the matrix
m <- matrix(1:12,nrow=3,ncol=4)

dim(m)
length(m)

# Give names to lines or columns
dimnames(m) <- list(paste("row",1:3,sep=""),
                    paste("column",1:4,sep=""))
m

# Extracting elements - this time we have two dimensions
# We must specify the row, and the column
m[3,4]

# To select all the data on the second row
m[2,]

# To select all the data in the third column
m[,3]

# You can select data from multiple rows and columns
m[c(1,2),c(1,3)]

# Operations are possible on matrices, the same way as vectors
1+m

sqrt(m)

# You can extract the diagonal of a matrix
diag(m)

# Or transpose a matric
t(m)

# You can also perform matrix multiplication
# (To get help on a special operator : ?"%*%"
t(m) %*% m

# It is also possible to solve matrices or systems
m2 <- matrix(c(1,0,1,1),ncol=2)
solve(m2)
solve(m2, c(1,-1))


###
# 6. DATA FRAMES ----------------------------------------------
###

# A dataframe is a mixture of a list and a matrix : it has the shape of a matrix but the columns can be of different classes

# Create a dataframe from scratch - columns can be seen as variables and rows as observations
mydf <- data.frame(firstCol = rep(c("A","B"),5),
                   secondCol = 1:10,
                   thirdCol = factor(rep(c(0,1), each = 5), levels = c(0,1),
                                     labels = c("Failure","Success")),
                   fourthCol = sample(1:100,size = 10))

# Most of the time, we will not create a dataframe but import it.

# importing the file Baby_health.csv which is in the Data folder (you first need to set the working directory to use the following syntax)
# help(read.table)

# It is always important to check the structure of your data before importing it
# Imported character variables are strings, this can be changed with the stringsAsFactors argument

# Complete path
baby <- read.table("C:/Users/vanbenedena/Documents/Instats/Rbase part 1 training material/Data/Baby_health.csv",
                   header = TRUE,
                   sep = ";",
                   dec = ",", stringsAsFactors = FALSE)

# If working directory is already specified, you can specify the relative path
getwd()
baby <- read.table("Data/Baby_health.csv",
                   header = TRUE,
                   sep = ";",
                   dec = ",", stringsAsFactors = FALSE)

# read.csv2 could work here since it is a special case of the read.table function: names are in the file, the comma is the decimal symbol and the semicolon is the separator symbol
# See help(read.csv) or ?read.csv
baby <- read.csv2("Data/Baby_health.csv")

# It is also possible to import Excel files, but that requires loading a package
library(openxlsx)
baby <- read.xlsx("Data/Baby_health.xlsx")


# Check information about our dataframe
str(baby)
View(baby)

# You can see your dataframe by clicking on it in the Global Environment or using View()

# column names
colnames(baby)
names(baby)
# row names
rownames(baby)
# row and column names
dimnames(baby)

# display the first 6 rows (by default)
head(baby)
# display the last 6 rows  (by default)
tail(baby)

# show the content of all the "cells" on the first row
baby[1,]
# show the content of all the "cells" on the first and third row
baby[c(1,3),]

# show the content of all the "cells" on the first and second columns
baby[,c(1,2)]
# show the content of all the "cells" on the second column
baby[,2]
baby[, "Sex"]
baby$Sex

# show the content of the "cell" on the first row and second column
baby[1,2]
baby[1,"Sex"]

# show the content of the "cells" on the first row and second and fourth columns
baby[1, c(2,4)]
baby[1,c("Sex","Weight")]

# show the content of the "cells" on the 4th, 5th and 6th rows and 2nd and 4th columns
baby[4:6, c(2,4)]


# is a vector
baby[[2]]
# !!! second column as a data.frame: baby[2]

# dimensions
dim(baby)
# number of rows
nrow(baby)
# number of columns
ncol(baby)

# create new variables/columns
baby$Weight.kg <- baby$Weight / 1000
baby <- within(baby, Weight.kg2 <- Weight/1000)


# recode a variable with the function recode of the car package
library(car)
help(recode)
table(baby$Pollution_household)
baby$Pollution_household2 <- recode(baby$Pollution_household, 
       '1 = "No pollution";
       2:3 = "Pollution" ', as.factor = TRUE)

table(baby$Pollution_household2)

# simpler solution if only two values in the new variable like here: 
baby$Pollution_household3 <- ifelse(baby$Pollution_household == 1, 
                                    yes = "No pollution", 
                                    no = "Pollution")

# take a subset according to a condition
baby.girl <- subset(baby, Sex == 1)
baby.girl2 <- baby[baby$Sex ==1, ]

# split the data frame in a list based on categories
split(baby, baby$Sex)

#   6.1 Add columns ---------------------------------------
# add new columns from another data frame

# Let's recreate our factor f1 if it is not in our global environment
# f1 <- factor(c("t1","t4","t1","t3","t4","t1","t3","t3"))


# df1
df1 <- data.frame(person=1:8,
                        treatment=f1,
                        duration=rep(c(2,1,4,7),c(2,3,2,1))
)
# df2
df2 <- read.csv2("Data/data2.csv")

df1; df2
names(df1); names(df2)
# df1$person and df2$patient contain the same information

# add columns
# help(merge)
df12 <- merge(x = df1, y = df2, by.x = "person", by.y = "patient", all = TRUE)
head(df12)
df12_bis <- merge(x = df1, y = df2, by.x = "person", by.y = "patient")
head(df12_bis)


# ADVANCED compare ID sets in df1 and df2
# x %in% y: takes each element of x and checks whether it belongs to y
# df1$person %in% df2$patient
# df2$patient %in% df1$person 

#   6.2 Add rows -------------------------
# add new rows here from another data frame
# "stacking" / "concatenation"

# df3
df3 <- read.csv2("Data/data3.csv")

names(df12); names(df3)

# add rows - error
df <- rbind(df12, df3)
# ! double-check names to avoid overwrite or confusion
# ! check your data structure before anything else 


df4<-df3[,c("person","treatment","duration")]
df5<-df12[,c("person","treatment","duration")]

dim(df4);dim(df5)
names(df4);names(df5)

mydf <- rbind(df4, df5)
mydf
dim(mydf)


###
# 7. EXPORT AND SAVE ----------------------------------------
###

# SAVE OBJECT(S) on FILE

# save in RData format
save(df1, file = "Data/myData.RData")

# load an RData file
# LOADS OBJECT IN YOUR ENVIRONMENT 
load(file = "Data/myData.RData")

# in a text format
write.table(df1, file = "Data/df.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# in csv
write.csv(df1, file = "Data/df.csv", row.names = FALSE, quote = FALSE)

# in excel
# install.packages("openxlsx")
# library(openxlsx)
write.xlsx(df1, file = "Data/df.xlsx")

