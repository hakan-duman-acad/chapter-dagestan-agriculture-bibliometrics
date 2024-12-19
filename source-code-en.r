n_packages <- c("tidyverse", "bibliometrix","clipr","ggthemes")
packages <- n_packages[!(n_packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(packages)

library(bibliometrix)
library(tidyverse)
library(ggthemes)

M <- convert2df(file = "./data/wos.bib", dbsource = "wos", format = "bibtex")
M <- M |> filter(PY <= 2023)

results <- biblioAnalysis(M, sep =";")

summary(results)

results$MostCitedPapers[1:15,]

CR <- localCitations(M, sep = ";")

print(CR$Papers[1:15,])

print(CR$Authors[1:10,])

DF <- dominance(results, k = 10)
print(DF)

authors=gsub(","," ",names(results$Authors))
indices <- Hindex(M, field = "author", elements=authors, sep = ";", years = 70)

indices$H |> arrange(desc(h_index)) |> head(10)

indices$H |> arrange(desc(g_index)) |> head(10)

indices$H |> arrange(desc(m_index)) |> head(10)

topAU <- authorProdOverTime(M, k = 10, graph = TRUE)

L <- lotka(results)
## Printing Beta Coefficient
L$Beta
## Printing Constant Coefficient
L$C
## Printing  Goodness of Fit Score
L$R2
## Printing p-value for Two-sample Kolmogorov-Smirnov Test (Beta=2)
L$p.value

## Plotting Scientific Productivity with Theoretical Lotka Law Values
Observed=L$AuthorProd[,3]
Theoretical=10^(log10(L$C)-2*log10(L$AuthorProd[,1]))
plot(L$AuthorProd[,1],Theoretical,type="l",col="red",ylim=c(0, 1), 
xlab="Articles", ylab="Freq. of Authors",main="Scientific Productivity")
lines(L$AuthorProd[,1],Observed,col="blue")
legend(x="topright",c("Theoretical (B=2)","Observed"),col=c("red","blue"),
lty = c(1,1,1),cex=0.6,bty="n")

A <- cocMatrix(M, Field = "SO", sep = ";")

a <-sort(Matrix::colSums(A), decreasing = TRUE)[1:10]
print(a)

dt <- tibble(Journal = names(a), Count = a)  |>
	mutate(Journal = str_to_title(Journal))  
dt  |>
	ggplot(aes(y=Count,x = reorder(Journal, Count)))+
	geom_bar(stat="identity")+
	coord_flip() + xlab("Journal Name") + ylab("Document Count")

b <- bradford(M)

b$table  |> filter(Zone == "Zone 1")

b$graph + ggtitle("") + xlab("")

NetMatrix <- biblioNetwork(M, analysis = "co-citation", network = "references", sep = ";")
net=networkPlot(NetMatrix, n = 75, Title = "",
             	size=TRUE, remove.multiple=FALSE,
            	labelsize=0.3)

netstat <- networkStat(NetMatrix, stat = "network", type = "degree")
summary(netstat)

fn <- "co-cite-refs.txt"
tmp <- net$cluster_res
tmp <- tmp |> group_by(cluster) |>
	arrange(cluster,desc(btw_centrality)) |>
	mutate(BC=1:n()) |>
	arrange(cluster,desc(clos_centrality)) |>
	mutate(CC=1:n()) |>
	arrange(cluster,desc(pagerank_centrality)) |>
	mutate(PR=1:n())
tmp$MD <- apply(tmp |> select(c(BC,CC,PR)) ,
                            	1, median)
tmp <- tmp |> group_by(cluster) |> arrange(cluster,MD) |>
	mutate_at(vars(4,5),\(x) round(x,5)) |>
	mutate_at(vars(3),\(x) format(round(x,3), nsmall = 3)) |>
	slice(1:10)
tmp <- tmp |>
	mutate(Paper = str_to_upper(vertex)) |>
	left_join(CR$Papers |>
            	select(Paper,DOI))
tmp |> clipr::write_clip(allow_non_interactive = TRUE)
write.table(file = fn,x = tmp)
for (i in 1:nrow(tmp)) {
	cat("\n\n\n-----------------------------------\n\n\n", file = fn, append = TRUE)
	name <- str_to_upper(str_remove_all(string = tmp$vertex[i], pattern = "\\d|-"))
	name <- str_remove(name,"\\s+$")
	year <- str_extract(string = tmp$vertex[i], pattern = "\\d{4}")
	str_tmp <- sprintf("
	Vertex name = %s; Cluster = %s;
	Betweenness Centrality (BC) = %.3f;  Closeness Centrality (CC) = %.6f; PageRank (PR) score = %.3f;
	BC rank = %s; CC Rank = %s; PR rank = %s \n
                   	", tmp[i,1], tmp[i,2], tmp[i,3], tmp[i,4], tmp[i,5],
                  	tmp[i,6], tmp[i,7], tmp[i,8])
	cat(str_tmp, file = fn, append = TRUE)
	doi <- tmp[i,11] |> unlist()
	cat("DOI :", file = fn, append = TRUE)
	cat(doi, file = fn, append = TRUE)
    cat("\n\nABSTRACT:\n", file = fn, append = TRUE)
    M |> filter(DI == tmp[i,11] |> unlist()) |> select(AB) |> 
        mutate(AB = str_to_sentence(AB)) |> unlist() |> cat(file = fn, append = TRUE)
}

NetMatrix <- biblioNetwork(M, analysis = "coupling", network = "references", sep = ";", short = TRUE,
  shortlabel = FALSE)
net=networkPlot(NetMatrix,  normalize = "association",
            	label.n = 20,#halo=TRUE,
            	weighted=T, n = 75,
            	Title = "",
            	type = "fruchterman",
            	size = TRUE, edgesize = 3,
            	labelsize = 0.3)

netstat <- networkStat(NetMatrix, stat = "network", type = "degree")
summary(netstat)

fn <- "coupling-list.txt"
tmp <- net$cluster_res
tmp <- tmp |> group_by(cluster) |>
	arrange(cluster,desc(btw_centrality)) |>
	mutate(BC=1:n()) |>
	arrange(cluster,desc(clos_centrality)) |>
	mutate(CC=1:n()) |>
	arrange(cluster,desc(pagerank_centrality)) |>
	mutate(PR=1:n())
tmp$MD <- apply(tmp |> select(c(BC,CC,PR)) ,
                            	1, median)
tmp <- tmp |> group_by(cluster) |> arrange(cluster,MD) |>
	mutate_at(vars(4,5),\(x) round(x,5)) |>
	mutate_at(vars(3),\(x) format(round(x,3), nsmall = 3)) |>
	slice(1:10)
tmp <- tmp |>
	mutate(Paper = str_to_upper(vertex)) |>
	left_join(CR$Papers |>
            	select(Paper,DOI))
tmp |> clipr::write_clip(allow_non_interactive = TRUE)
write.table(file = fn,x = tmp)
for (i in 1:nrow(tmp)) {
	cat("\n\n\n-----------------------------------\n\n\n", file = fn, append = TRUE)
	name <- str_to_upper(str_remove_all(string = tmp$vertex[i], pattern = "\\d|-"))
	name <- str_remove(name,"\\s+$")
	year <- str_extract(string = tmp$vertex[i], pattern = "\\d{4}")
	str_tmp <- sprintf("
	Vertex name = %s; Cluster = %s;
	Betweenness Centrality (BC) = %.3f;  Closeness Centrality (CC) = %.6f; PageRank (PR) score = %.3f;
	BC rank = %s; CC Rank = %s; PR rank = %s \n
                   	", tmp[i,1], tmp[i,2], tmp[i,3], tmp[i,4], tmp[i,5],
                  	tmp[i,6], tmp[i,7], tmp[i,8])
	cat(str_tmp, file = fn, append = TRUE)
	doi <- tmp[i,11] |> unlist()
	cat("DOI :", file = fn, append = TRUE)
	cat(doi, file = fn, append = TRUE)
    cat("\n\nABSTRACT:\n", file = fn, append = TRUE)
    M |> filter(DI == tmp[i,11] |> unlist()) |> select(AB) |> 
        mutate(AB = str_to_sentence(AB)) |> unlist() |> cat(file = fn, append = TRUE)
}

## Calculating and Plotting Collaboration Between Countries Network


# NetMatrix <- biblioNetwork(M, analysis = "collaboration", network = "countries", sep = ";")
# net=networkPlot(NetMatrix, n = 75, Title = "",
#              	size=TRUE, remove.multiple=TRUE,
#             	labelsize=0.3,cluster="none")

lst <- M |> count(PY) |> 
    mutate(tot = cumsum(n), 
           cut_points = cut(
               tot, breaks = seq(0,to = nrow(M)+1,length.out = 5),
               labels = 1:4))
# print(lst)

syns <- c("dagestan;daghestan","chechnya;chechnia")
r.terms <- NULL
thema <- function(dt = M) {
    res <- thematicMap(dt, field = "TI", n = 100, minfreq = 5, size = 0.5, n.labels = 5, repel = TRUE, synonyms = syns, remove.terms = r.terms)
    return (res)
}


for(i in 1:2) {
    lst2 <- lst$PY[lst$cut_points == i]
    lab <- sprintf("Time Slice = %d (Studies published between %s and %s)",i,head(lst2,1),tail(lst2,1))
    res <- thema(M|> filter(PY %in% lst2))
    plt <- res$map +
        ggtitle (subtitle = lab,label = "" )
    plot(plt)
    cat("\n\n Time Slice :",i,"\n\n")
    res$clusters |> as.data.frame() |> print()
    tmp <- res$documentToClusters |> left_join(M) |> 
    arrange(Assigned_cluster) |>
    select(Assigned_cluster,TI,DI,UT,AB) |>
    group_by(Assigned_cluster)
    tmp |> count(Assigned_cluster) |>
    print()
    tmp |>
    sample_n(10,replace = TRUE) |>
    distinct() |>
    write_csv(paste0("ts-",i,".csv"))
}


