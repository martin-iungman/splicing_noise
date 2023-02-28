library(tidyverse)
library(furrr)

list.files(path="Amplicon-seq", pattern = "tsv$", full.names=T)
list_gates<-map(list.files(path="Amplicon-Seq",pattern = "tsv$",full.names=T),
                ~read_tsv(.x)%>%rowwise()%>%mutate(tot=sum(c_across(where(is.double)))))
names(list_gates)<-seq(1,5)

map(list_gates, ~sum(.x$tot))%>%list_c()%>%tibble(.name_repair = "universal")%>%rownames_to_column(var="Gates")%>%rename(counts=".")%>%
  ggplot(aes(Gates,counts))+geom_col(fill="#3FA0FF")+ggtitle("Number of counts per Gate")+theme_bw()

counts_per_rep<-map2(list_gates, names(list_gates),~.x%>%as_tibble()%>%summarise(across(starts_with("rep"), sum))%>%mutate(Gate=as.factor(.y)))%>%list_rbind()%>%
  pivot_longer(where(is.double), names_to="rep", values_to="counts")
counts_per_rep%>%
  ggplot(aes(Gate, counts))+geom_violin()+
  geom_point(aes(col=rep), size=2)+
  ggtitle("Number of counts per Gate")+
  scale_color_brewer(palette = "Set1", type="qual")+theme_bw()

counts_per_rep%>%ggplot(aes(rep, counts, fill=Gate))+
  geom_col(position = "fill")+
  scale_fill_brewer(palette ="Set1",type="qual")+
  theme_bw()+xlab("Replicate")


# Assesing the consistency of replicates distributions
min_range_raw=c()
max_range_raw=c()
min_range_norm=c()
max_range_norm=c()
raw_plots<-list()
norm_plots<-list()
for(i in 1:5){
  tots<-counts_per_rep%>%filter(Gate==i)
  a<-list_gates[[i]]%>%pivot_longer(starts_with("rep"), values_to = "counts_rep",names_to = "Replicates")%>%left_join(tots, by=c("Replicates"="rep"))%>%
    mutate(counts_rep_norm=counts_rep/counts, counts_tot_norm=tot/counts)
  min_range_raw<-min(min_range_raw,a$counts_rep)
  max_range_raw<-max(max_range_raw, a$counts_rep)
  min_range_norm<-min(min_range_norm,a$counts_rep_norm)
  max_range_norm<-max(max_range_norm, a$counts_rep_norm)
  plot_norm<-a%>%ggplot(aes(counts_rep_norm, fill=Replicates))+geom_density(alpha=0.7)+ggtitle(paste0("Gate ",i, ": Normalized"))+scale_color_brewer(palette = "Set1", type="qual")+theme_bw()+theme(legend.position = "none", plot.title = element_text(size = 10), axis.title = element_text(size=7))+xlab("Normalized Counts")
  #assign(paste0("density_norm_gate",i), plot_norm, envir = parent.frame())
  norm_plots[[i]]<-plot_norm
  plot<-a%>%ggplot(aes(counts_rep, fill=Replicates))+geom_density(alpha=0.7)+ggtitle(paste0("Gate ",i, ": Raw"))+scale_color_brewer(palette = "Set1", type="qual")+theme_bw()+theme(legend.position = "none", plot.title = element_text(size = 10), axis.title = element_text(size=7))+xlab("Raw Counts")
  #assign(paste0("density_raw_gate",i), plot, envir = parent.frame())
  raw_plots[[i]]<-plot
}
raw_plots<-map(raw_plots, ~.x+xlim(min_range_raw, max_range_raw)+scale_x_log10())
norm_plots<-map(norm_plots, ~.x+xlim(min_range_norm, max_range_norm)+scale_x_log10())

legend_plot<-plot_norm+theme(legend.position="bottom")
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}
shared_legend <- extract_legend(legend_plot)
grid.arrange(grobs=c(raw_plots, norm_plots, list(shared_legend)),layout_matrix=matrix(c(1:5,11,6:10,11), ncol=2, byrow = F))

# grid.arrange(
#   arrangeGrob(arrangeGrob(density_raw_gate1, density_raw_gate2, density_raw_gate3, density_raw_gate4, density_raw_gate5, nrow=5),
#   arrangeGrob(density_norm_gate1, density_norm_gate2, density_norm_gate3, density_norm_gate4, density_norm_gate5, nrow=5),ncol=2), 
#   shared_legend, ncol=1, heights=c(1, 0.1)
# )

grid.arrange(grobs=arrangeGrob(raw_plots, norm_plots), ncol=2)


# list_gates[[1]]%>%ggplot(aes(tot))+geom_density()+scale_x_log10()
# list_gates[[1]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 1")
# list_gates[[2]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 2")
# list_gates[[3]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 3")
# list_gates[[4]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 4")
# list_gates[[5]]%>%pivot_longer(starts_with("rep"), values_to = "Counts",names_to = "Replicates")%>%ggplot(aes(Counts, fill=Replicates))+geom_density(alpha=0.7)+scale_x_log10()+ggtitle("Gate 5")

long_df<-map(seq(1,5),~list_gates[[.x]]%>%select(seqname, tot)%>%mutate(gate=as.factor(.x)))%>%list_rbind()

long_df%>%
  ggplot(aes(tot,fill=gate))+geom_density(alpha=0.5)+scale_x_log10()+scale_fill_brewer(palette="Set1", type="qual")+theme_bw()

long_df%>%group_by(seqname)%>%summarise(suma=sum(tot))%>%ggplot(aes(suma))+geom_density()+scale_x_log10()+
  ggtitle("Total counts per sequence")+ theme_bw()
  

get_stats<-function(df){
  k<-nrow(df)
  df$gate<-as.numeric(df$gate)
  if(k!=5){
    df<-bind_rows(df, data.frame(seqname=unique(df$seqname),tot= 1, gate=seq(1,5)[!seq(1,5)%in%df$gate]))
  }
  a<-sum(df$tot)
  vctr<-map(1:5, ~rep(df$gate[.x],times=df$tot[.x]))%>%list_c()
  entropy=-1/log(5)*sum(df$tot/a*log(df$tot/a))
  extropy=(-1/(4*log(5/4)))*sum((1-df$tot/a)*log(1-df$tot/a))
  return(data.frame(seqname=unique(df$seqname),positive_gates=k,sum=a,mean=mean(vctr),median=floor(median(vctr)), sd=sd(vctr),mad=mad(vctr, constant=1), entropy=entropy, extropy=extropy))
}


plan(multisession, workers = 2)
list_stats<-future_map(unique(long_df$seqname), ~long_df[long_df$seqname==.x,]%>%get_stats(), .progress = T)%>%list_rbind()

list_stats%>%ggplot(aes(extropy))+geom_density()
list_stats%>%ggplot(aes(entropy))+geom_density()
list_stats%>%ggplot(aes(median))+geom_histogram()
list_stats%>%ggplot(aes(mean))+geom_density()
list_stats%>%ggplot(aes(sd))+geom_density()
list_stats%>%ggplot(aes(extropy, fill=as.factor(floor(median))))+geom_density(alpha=0.7)

seqnames<-map(list_gates, ~.x%>%select(seqname, type))%>%list_rbind()%>%unique()
list_stats<-left_join(list_stats, seqnames, by="seqname")
list_stats%>%ggplot(aes(floor(median), fill=type))+geom_histogram(position = "fill")
list_stats%>%ggplot(aes(extropy, fill=type))+geom_density(alpha=0.7)
list_stats%>%filter(type%in%c("cancer-promoter-wt","cancer-promoter-mut"))%>%ggplot(aes(extropy, fill=type))+geom_density(alpha=0.7)
wilcox.test(list_stats%>%filter(type%in%c("cancer-promoter-wt"))%>%select(extropy)%>%as_vector(), list_stats%>%filter(type%in%c("cancer-promoter-mut"))%>%select(extropy)%>%as_vector())
list_stats%>%ggplot(aes(sum,extropy))+geom_point(size=0.7)+scale_x_log10()
list_stats%>%ggplot(aes(sum,entropy))+geom_point(size=0.7)+scale_x_log10()


#get_stats(long_df%>%filter(seqname==list_stats$seqname[15198]))
long_df%>%filter(seqname==list_stats$seqname[12331])%>%ggplot(aes(gate,tot))+geom_col()


