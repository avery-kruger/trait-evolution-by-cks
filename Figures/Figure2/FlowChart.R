library(viridis)
library(ggplot2)



# Dummy data
set.seed(586)
a <- 1:10
b <- paste0("Obs. ", seq(1,10))
data <- expand.grid(X=a, Y=b)
data$Z <- rnorm(nrow(data))
data$Z[1:50] <- data$Z[51:100]

# Give extreme colors:
ggplot(data[1:50,], aes(X, Y, fill= Z)) +
  geom_tile(color="black") +
  scale_fill_viridis(name = "Value") +
  scale_y_discrete("Observation") +
  scale_x_discrete("Data", expand = c(0.025,0.025)) +
  geom_segment(aes(x=0.5,y=4.5,xend=10.5,yend=4.5), size = 1) +
  geom_segment(aes(x=0.5,y=5.5,xend=10.5,yend=5.5), size = 1) +
  geom_segment(aes(x=10.5,y=4.4875,xend=10.5,yend=5.5125), size = 1) +
  geom_segment(aes(x=0.5,y=4.4875,xend=0.5,yend=5.5125), size = 1) +
  geom_segment(aes(x=1,y=5.35,xend=6,yend=5.35), size = 1.5, color = "#fcfdbf") +
  geom_segment(aes(x=2,y=5.175,xend=7,yend=5.175), size = 1.5, color = "#fc8961") +
  geom_segment(aes(x=3,y=5.0,xend=8,yend=5.0), size = 1.5, color = "#b73779") +
  geom_segment(aes(x=4,y=4.825,xend=9,yend=4.825), size = 1.5, color = "#51127c") +
  geom_segment(aes(x=5,y=4.65,xend=10,yend=4.65), size = 1.5, color = "#000004") +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.y = element_text(size = 11),)
ggsave("FlowData.png")

# Dummy data
# c <- 1:20
# d <- paste0("obs", seq(1,20))
# data <- expand.grid(X=c, Y=d)
# data$Z <- runif(400, 0, 5)

mymatrix <- as.data.frame(matrix(data$Z[91:96],nrow=1))
for(i in 2:5){
  mymatrix[i,] <- data$Z[90+i:(i+5)]
}
mymatrix <- as.matrix(mymatrix)

windows <- expand.grid(Y=5:1,X=1:6)
windows$Z <- as.vector(mymatrix)
ggplot(windows, aes(X, Y, fill= Z)) +
  geom_tile(color="black") +
  scale_fill_viridis(name = "", limits=c(min(data$Z),max(data$Z))) +
  scale_y_discrete("Window") +
  scale_x_discrete("Window Element") +
  geom_segment(aes(x=1,y=5.4,xend=6,yend=5.4), size = 2, color = "#fcfdbf") +
  geom_segment(aes(x=1,y=4.4,xend=6,yend=4.4), size = 2, color = "#fc8961") +
  geom_segment(aes(x=1,y=3.4,xend=6,yend=3.4), size = 2, color = "#b73779") +
  geom_segment(aes(x=1,y=2.4,xend=6,yend=2.4), size = 2, color = "#51127c") +
  geom_segment(aes(x=1,y=1.4,xend=6,yend=1.4), size = 2, color = "#000004") +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.y = element_text(size = 11),)
ggsave("FlowWindows.png")

set.seed(1212)
mynormal <- matrix(rnorm(60), nrow=6)
normplot <- expand.grid(Y=6:1,X=1:10)
normplot$Z <- as.vector(mynormal)
ggplot(normplot, aes(X, Y, fill= Z)) +
  ggtitle("Normal Matrix") +
  geom_tile(color="black") +
  scale_fill_viridis(name = "", limits=c(min(data$Z),max(data$Z))) +
  scale_y_discrete("", expand = c(0,0)) +
  scale_x_discrete("", expand = c(0,0)) +
  theme_classic() +
  theme(axis.line = element_blank(),
        text = element_text(size=15),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5))
ggsave("FlowNormal.png")

mymult <- mymatrix %*% mynormal
temp <- expand.grid(X=5:1,Y=1:10)
temp$Z <- as.vector(mymult)

# Give extreme colors:
ggplot(temp, aes(Y, X, fill= Z)) +
  geom_tile() +
  scale_fill_viridis(name = "") +
  scale_y_discrete("") +
  scale_x_discrete("") +
  theme_classic() +
  theme(axis.line = element_blank())
ggsave("FlowWindowXNorm.png")

myclamp <- mymult
myclamp[myclamp < 0] <- 0
temp$clamp <- as.vector(myclamp)
ggplot(temp, aes(Y, X, fill= clamp)) +
  geom_tile() +
  scale_fill_viridis(name = "", limits=c(min(temp$Z),max(temp$Z))) +
  scale_y_discrete("") +
  scale_x_discrete("") +
  theme_classic() +
  theme(axis.line = element_blank())
ggsave("FlowClamp.png")

mymeans <- colMeans(myclamp)
newtemp <- expand.grid(X=1:10,Y=1:5)
newtemp$mean <- NA
newtemp$mean[21:30] <- mymeans


mydata <- mymeans
for(j in 1:4){
  mymatrix <- as.data.frame(matrix(data$Z[(91:96)-10*j],nrow=1))
  for(i in 2:5){
    mymatrix[i,] <- data$Z[90-10*j+i:(i+5)]
  }
  mymatrix <- as.matrix(mymatrix)
  mymult <- mymatrix %*% mynormal
  myclamp <- mymult
  myclamp[myclamp < 0] <- 0
  mydata <- rbind(mydata,colMeans(myclamp))
}

allfeatures <- expand.grid(Y=paste0("Obs. ",1:5),X=1:10)
allfeatures$Z <- as.vector(mydata[c(5,4,3,2,1),])
allfeatures$Z[allfeatures$Z>1.5] <- allfeatures$Z[allfeatures$Z>1.5]*0.5

ggplot(newtemp, aes(X, Y, fill= mean)) +
  geom_tile() +
  scale_fill_viridis(name = "", na.value="white",
                     limits = c(min(allfeatures$Z),max(allfeatures$Z))) +
  scale_y_discrete("") +
  scale_x_discrete("") +
  theme_classic() +
  theme(axis.line = element_blank())
ggsave("FlowOutput.png")

ggplot(allfeatures, aes(X, Y, fill= Z)) +
  geom_tile(color="black") +
  scale_fill_viridis(name = "Value") +
  scale_y_discrete("Observation") +
  scale_x_discrete("Features", expand = c(0.025,0.025)) +
  geom_segment(aes(x=0.5,y=4.5,xend=10.5,yend=4.5), size = 1) +
  geom_segment(aes(x=0.5,y=5.5,xend=10.5,yend=5.5), size = 1) +
  geom_segment(aes(x=10.5,y=4.4875,xend=10.5,yend=5.5125), size = 1) +
  geom_segment(aes(x=0.5,y=4.4875,xend=0.5,yend=5.5125), size = 1) +
  theme_classic() +
  theme(text = element_text(size = 15),
        axis.text.y = element_text(size = 11),)
ggsave("FlowAllFeatures.png")
