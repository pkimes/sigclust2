#' @description a set of simple examples to check basic functionality of HSigClust package
#' 
#' @author Patrick Kimes
#' @title simple toy examples
#' @name toy_examples


#I need to get use to rbenchmark'ing my code
library(rbenchmark)
library(Rclusterpp)
# library(fastcluster) #improves hclust, the clustering not the dist calculation

N <- 100
p <- 10
delta <- .8
x <- rbind(matrix(rnorm(n=N*p, mean=c(-delta, 0), sd=1), nrow=N, ncol=p, byrow=TRUE),
           matrix(rnorm(n=N*p, mean=c(delta, 0), sd=1), nrow=N, ncol=p, byrow=TRUE),
           matrix(rnorm(n=N*p, mean=c(0, delta*sqrt(3)), sd=1), nrow=N, ncol=p, byrow=TRUE))
labels <- rep(1:3, each=N)

#looks like it at least runs and produces p-values which seem reasonable.
firstrun <- HSCtest(x, metric="euclidean", linkage="average", square=TRUE, icovest=2)


#read in dataset
n <- nrow(x)
ids <- paste("sample", 1:n) 
rownames(x) <- ids
  
#using correlation dist.
mycordist <- as.dist(1 - cor(t(x)))
hc_Joel <- hclust(mycordist, method="average")
hcd_Joel <- as.dendrogram(hc_Joel)
#using ward's method
hc_Ward <- hclust(dist(x, "euclidean")^2, "ward")
hcd_Ward <- as.dendrogram(hc_Ward)

labelColors <- c("#CDB380", "#036564", "#EB6841", "#EDC951")
labelColors <- c("#D50065", "#FF5F00", "#580EAD", "#8EEB00", "#052D6E")

# function to get color labels
colLab <- function(n) {
  a <- attributes(n)
  attr(n, "edgePar") <- c(a$edgePar, list(col="#4D4F57", lty=c(2, 2)))
  if (is.leaf(n)) {
#     a <- attributes(n)
    labCol <- labelColors[labels[which(ids == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, 
                            list(pch=15, col=labCol, cex=.3,
                                 lab.cex=1, lab.col=labCol))
    attr(n, "label") <- "--"
#     attr(n, "edgePar") <- c(a$nodePar, list(pch=15, col=labCol, cex=.3))
  }
  n
}
  
# function to get color labels
colLab <- function(n) {
  a <- attributes(n)
  attr(n, "edgePar") <- c(a$edgePar, list(col=c("#4D4F57","#4D4F57"), lty=c(2, 1), 
                                          p.lty=c(3, 3), p.lwd=c(2, 2)))
  labCol <- labelColors[labels[which(ids == a$label)]]
  attr(n, "nodePar") <- c(a$nodePar, 
                          list(pch=15, col=labCol, cex=.3,
                               lab.cex=1, lab.col=labCol))
  attr(n, "label") <- "--"
  n
}

clusMember = cutree(hc_Ward, 4)
colLab <- function(n) {
  a <- attributes(n)
  if (is.leaf(n)) {
    attr(n, "edgePar") <- c(a$edgePar, list(col="#4D4F57", lty=1, 
                                            p.lty=3, p.lwd=2))
    labCol <- labelColors[clusMember[which(names(clusMember) == a$label)]]
    attr(n, "nodePar") <- c(a$nodePar, 
                            list(pch=15, col=labCol, cex=1,
                                 lab.cex=1, lab.col=labCol))            
    attr(n, "label") <- "--"
  } else {
    if (a$height >= hc_Ward$height[26]) {
      attr(n, "edgePar") <- c(a$edgePar, list(col="#4D4F57", lty=1, 
                                              p.lty=1, p.lwd=1,
                                              t.col="#D0A7A7", t.cex=.75))      
      attr(n, "edgetext") <- "p=5.3e-4"
    }
#     attr(n, "nodePar") <- c(a$nodePar, 
#                             list(pch=13, col="#000000", cex=1,
#                                  lab.cex=.5, lab.col="#000000"))
    
  }
  n
}

# using dendrapply
clusDendro <- dendrapply(hcd_Ward, colLab)
# attr(clusDendro[[2]], "edgePar") <- list(col="red", lwd=2)
# attr(clusDendro[[1]], "edgePar") <- list(col="red", lwd=2)


# make plot
# par(bg = "#FFFFFF")
plot(clusDendro, main = "Prettier Dendrogram", dLeaf=10, root.node=FALSE) #leaflab="none", pch=NA)

ggdendro::ggdendrogram(hc, rotate=FALSE)
hcdata <- ggdendro::dendro_data.hclust(hc, type="rectangle")
ggplot() + 
  geom_segment(data=ggdendro::segment(hcdata), aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_text(data=ggdendro::label(hcdata), 
            aes(x=x, y=y-10, label=label, hjust=1, vjust=.5, angle=90), 
            size=3) +
  theme_bw() + ylim(-200, 500)#scale_y_continuous(expand=c(0.35, 0))




