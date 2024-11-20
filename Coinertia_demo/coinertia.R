#please install miRNAtap, made4, plotly packages for the coinertia plot
#also may need Seurat, pracma,  dplyr package 
#setwd("~/data/CaiT/sRNA/PSCSRII/Revision2/")
library(Seurat)
library(pracma)

#process the expression infor
#load("merge_mRNA.dat")
#cells2<-cells
#Idents(cells2)<-"celltype"
#load("merge_miRNA.dat")
#cells1<-cells
#Idents(cells1)<-"celltype"
#data_mir<-as.matrix(AverageExpression(cells1)$RNA)
#data_rna<-as.matrix(AverageExpression(cells2)$RNA)

# or directly used the processed data
mirna<-read.csv( "data_mir.csv")
data_mir<-as.matrix(mirna[,2:dim(mirna)[2]])
rownames(data_mir)<-mirna[,1]
rna<-read.csv( "data_rna.csv")
data_rna<-as.matrix(rna[,2:dim(rna)[2]])
rownames(data_rna)<-rna[,1]

#filter data for coinertia analysis, too many miRNAs/mRNAs lead to very slow plot
summary(as.vector(data_mir))
inx<-apply(data_mir,1, function(x){any(x>2)})
sum(inx)
#transform data by reducing the weight of potential outliers
data.m<-sqrt(data_mir[inx,])

summary(as.vector(data_rna))
inx<-apply(data_rna,1, function(x){any(x>1.5)})
sum(inx)
data.g<-sqrt(data_rna[inx,])

library(made4)
#coinertia analysis
coin<- cia(data.m, data.g)
infor<-coin$coinertia
infor

# preprare the coinertia plot
library(dplyr)  
library(plotly)

#check if it works for your GUI
p1<-plot_ly(infor$co, x=~Comp1, y=~Comp2, color = (log(rowSums(data.m))), text=~rownames(infor$co), type = "scatter", mode="markers", hoverinfo = "text")
p1

#miRNA target prediction using miRNAtap package
library(miRNAtap)
library(org.Hs.eg.db)

mir= rownames(data.m)
mir = gsub("hsa-", "", mir)

idmaps = c()
for(i in 1:length(mir)){
  t=c()
  predictions = getPredictedTargets(mir[i], species = 'hsa', method = 'geom', min_src = 3)
  if(!is.null(predictions)) {
  tidmaps<-select(org.Hs.eg.db, rownames(predictions), c("ENTREZID", "SYMBOL"))
#  print(dim(tidmaps))
  if(sum(rownames(data.g) %in% tidmaps$SYMBOL)>0) {
    t= cbind(rownames(data.g)[rownames(data.g) %in% tidmaps$SYMBOL],paste("hsa.",  gsub("-", "." , mir[i]), sep="") )
    idmaps=rbind(idmaps, t)
  }
  }  
}

#separate data for plot
#miRNA targets
data2<-cbind(infor$li[idmaps[,1],], idmaps)
colnames(data2)<-c("x", "y", "name", "miRNA")

#other mRNAs
inx<-rownames(infor$li)[!rownames(infor$li) %in% data2$name]
data3<-cbind(infor$li[inx,], inx, "NA") 
colnames(data3)<-c("x", "y", "name", "miRNA")
data3<- data3[!grep("^A[A-Z]\\d+\\.\\d", rownames(data3)),]

#miRNA expression for plot
data1<-cbind(infor$co, rownames(infor$co), rownames(infor$co))
colnames(data1)<-c("x", "y", "name", "miRNA")
data<- rbind(data1, data2, data3)
datag <- data %>%  highlight_key(~miRNA) 
class(datag)

mytext1 = apply(round(data.m,2), 1, paste, collapse=",")
mytext2 = rownames(data.m)
mytext3 = paste(substr(colnames(data.m), 1,4), sep=",", collapse = ",")
mytext =paste("<br>", mytext2,  "<br>", mytext3, "<br>", mytext1, sep = ""  )
fig1 <- plot_ly(data=datag) %>% filter(grepl("hsa", name)) %>%
  add_markers(
    x = ~x, 
    y = ~y, 
    type  = "scatter",
    mode  = "markers",
    text = mytext,
    color = log(rowSums(data.m)),
    hoverinfo = "text"
  ) 

#fig1

#mRNA expression for plot
d<-data.g[c(as.character(data2$name),gsub("\\.([A-Z])", "-\\1", as.character(data3$name)) ), ]
mytext1 = apply(round(d,2), 1, paste, collapse=",")
mytext2 = rownames(d)
mytext3 = paste(substr(colnames(d), 1,4), sep=",", collapse = ",")
mytext =paste("<br>", mytext2,  "<br>", mytext3, "<br>", mytext1, sep = ""  )
fig2 <- plot_ly(data = datag) %>%  # initiate plot with same data frame
  filter(!grepl("hsa", name)) %>%                     # use dplyr verb on plotly object
  add_markers(
    x     = ~x,
    y     = ~y,
    text = mytext,
    color = log(rowSums(d) ),
    hoverinfo = "text"
    ) 
#fig2

#show global correlation
data3<-data.frame(x=infor$mX[,1], y=infor$mX[,2], x2 = infor$mY[,1], y2= infor$mY[,2], names = rownames(infor$mX))
fig3<- plot_ly(data3, text=~names) %>% add_markers(~x, ~y) %>% add_annotations( 
                   x = ~x,
                   y = ~y,
                   xref = "x", yref = "y",
                   axref = "x", ayref = "y",
                   text = ~names,
                   xanchor = 'left',
                   opacity= 0.5,
                   showarrow = T,
                   ax = ~x2,
                   ay = ~y2
)
#fig3

s1 <- subplot(fig1, fig2)%>%  
  layout(annotations = list( 
    list(x = 0.15 , y = 1.2, text = "miRNA Subplot", showarrow = F, xref='paper', yref='paper'), 
    list(x = 0.85 , y = 1.2, text = "mRNA Subplot", showarrow = F, xref='paper', yref='paper')) 
  )
#s1

fig <- subplot(fig3, s1,  nrows = 2, margin = 0.07) %>%  
  layout(title = "miRNA mRNA CIA analysis result", 
         plot_bgcolor='#e5ecf6', 
         showlegend=FALSE,showlegend2=FALSE, 
         margin = 0.01) 

fig %>%  highlight(on = "plotly_click",selectize = TRUE)

#save the plot 
library(htmlwidgets)
ciaplot<- fig %>%  highlight(on = "plotly_click",selectize = TRUE)
saveWidget(ciaplot, "coinertia_demo.html", selfcontained = T)
