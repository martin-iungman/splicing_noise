library(tidyverse)
library(furrr)

list.files(pattern = "tsv$")
list_gates<-map(list.files(path="Amplicon-Seq",pattern = "tsv$"),~read_tsv(.x)%>%rowwise()%>%mutate(tot=sum(c_across(where(is.double)))))
names(list_gates)<seq(1,5)

map_dfr(list_gates, ~sum(.x$tot))%>%t()%>%tibble(.name_repair = "universal")%>%rownames_to_column(var="Gates")%>%rename(counts=".")%>%
  ggplot(aes(Gates,counts))+geom_col()+ggtitle("Number of counts per Gate")

map2_dfr(list_gates, names(list_gates),~.x%>%as_tibble()%>%summarise(across(starts_with("rep"), sum))%>%mutate(Gate=.y))%>%pivot_longer(where(is.double), names_to="rep", values_to="counts")%>%
  ggplot(aes(Gate, counts))+geom_violin()+geom_point(aes(col=rep))+ggtitle("Number of counts per Gate")

seqnames<-map(list_gates, ~.x%>%select(seqname, type))%>%list_rbind()%>%unique()

list_gates[[1]]%>%ggplot(aes(tot))+geom_density()+scale_x_log10()
list_gates[[1]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 1")
list_gates[[2]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 2")
list_gates[[3]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 3")
list_gates[[4]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 4")
list_gates[[5]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 5")

long_df<-map(seq(1,5),~list_gates[[.x]]%>%select(seqname, tot)%>%mutate(gate=.x))%>%list_rbind()

long_df%>%
  ggplot(aes(tot,fill=Gates))+geom_density(alpha=0.7)+scale_x_log10()

long_df%>%group_by(seqname)%>%summarise(suma=sum(tot))%>%ggplot(aes(suma))+geom_density()+scale_x_log10()+
  ggtitle("Total counts per sequence")
  

get_stats<-function(df){
  k<-nrow(df)
  if(k!=5){
    df<-bind_rows(df, data.frame(seqname=unique(df$seqname),tot= 1, gate=seq(1,5)[!seq(1,5)%in%df$gate]))
  }
  a<-sum(df$tot)
  vctr<-map(1:k, ~rep(df$gate[.x],times=df$tot[.x]))%>%list_c()
  entropy=-1/log(5)*sum(df$tot/a*log(df$tot/a))
  extropy=(-1/(4*log(5/4)))*sum((1-df$tot/a)*log(1-df$tot/a))
  return(data.frame(seqname=unique(df$seqname),positive_gates=k,sum=a,mean=mean(vctr),median=median(vctr), sd=sd(vctr),mad=mad(vctr, constant=1), entropy=entropy, extropy=extropy))
}


plan(multisession, workers = 8)
list_stats<-future_map(unique(long_df$seqname), ~long_df[long_df$seqname==.x,]%>%get_stats(), .progress = T)%>%list_rbind()

list_stats%>%ggplot(aes(extropy))+geom_density()
list_stats%>%ggplot(aes(entropy))+geom_density()
list_stats%>%ggplot(aes(median))+geom_histogram()
list_stats%>%ggplot(aes(mean))+geom_density()
list_stats%>%ggplot(aes(sd))+geom_density()
list_stats%>%ggplot(aes(extropy, fill=as.factor(median)))+geom_density(alpha=0.7)

list_stats<-left_join(list_stats, seqnames, by="seqname")
list_stats%>%ggplot(aes(median, fill=type))+geom_histogram(position = "dodge")
list_stats%>%ggplot(aes(extropy, fill=type))+geom_density(alpha=0.7)
list_stats%>%filter(type%in%c("cancer-promoter-wt","cancer-promoter-mut"))%>%ggplot(aes(extropy, fill=type))+geom_density(alpha=0.7)
list_stats%>%ggplot(aes(sum,extropy))+geom_point(size=0.7)+scale_x_log10()
list_stats%>%ggplot(aes(sum,entropy))+geom_point(size=0.7)+scale_x_log10()


get_stats(long_df%>%filter(seqname==list_stats$seqname[15198]))
long_df%>%filter(seqname==list_stats$seqname[12331])%>%ggplot(aes(gate,tot))+geom_col()


