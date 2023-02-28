library(shiny)
library(tidyverse)
library(furrr)
library(data.table)
list.files(path="../Amplicon-seq", pattern = "tsv$", full.names=T)
plan(multisession, workers = 5)
list_gates<-future_map(list.files(path="../Amplicon-Seq",pattern = "tsv$",full.names=T),
                       ~read_tsv(.x)%>%rowwise()%>%mutate(tot=sum(c_across(where(is.double)))))
plan(sequential)
names(list_gates)<-seq(1,5)
long_df<-map(seq(1,5),~list_gates[[.x]]%>%select(seqname, tot)%>%mutate(gate=as.factor(.x)))%>%list_rbind()
print("files loaded")
get_stats<-function(df){
  k<-length(df$tot[df$tot!=0])
  df$gate<-as.numeric(df$gate)
  a<-sum(df$tot)
  vctr<-map(1:length(unique(df$gate)), ~rep(df$gate[.x],times=df$tot[.x]))%>%list_c()
  if(k==1){extropy=0}else{
    extropy=(-1/((k-1)*log(k/(k-1))))*sum((1-df$tot/a)*log(1-df$tot/a))
  }
  mean=mean(vctr)
  median=median(vctr)
  mode<-df$gate[which.max(df$tot)]
  cv2<-(sd(vctr)/mean(vctr))^2
  fano_factor<-var(vctr)/mean(vctr)
  mad=mad(vctr, constant=1)
  sum=a
  if(k!=5){
    vctr<-c(vctr,seq(1,5)[!seq(1,5)%in%df$gate])
    a=sum(a, 5-length(unique(df$gate)))}
  entropy=(-1/log(5))*sum(df$tot/a*log(df$tot/a))
  return(data.frame(seqname=unique(df$seqname),positive_gates=k,sum,mean,median, mode,cv2,fano_factor,mad, entropy, extropy))
}
plan(multisession, workers = 5)
list_stats<-future_map(unique(long_df$seqname), ~long_df[long_df$seqname==.x,]%>%get_stats(), .progress = T)%>%list_rbind()
print("stats calculated")
list_stats<-list_stats%>%mutate(across(where(is.numeric),round,3))


ui <- fluidPage(
  tags$head(
    tags$link(href = "https://fonts.googleapis.com/css?family=Roboto+Mono", rel = "stylesheet"),
    tags$style(HTML('
      * {
        font-family: Roboto Mono;
        font-size: 100%;
      }
      #sidebar {
         background-color: #a3c1ad;
         border: 0px;
      }
      .rt-th {
        display: none;
      }
      .rt-noData {
        display: none;
      }
      .rt-pagination-nav {
        float: left;
        width: 1%;
      }
    '))
  ),

      sidebarLayout(
        sidebarPanel(
          id = "sidebar",
          tags$h1("Splicing Outcome Assay"),
          tags$br(),
          textInput("prom",label="Promoter ID to plot:",value="FP010289_CENPW_2"),
          tags$br(),
          
          tabsetPanel(
            tabPanel(
              "Histogram", value = "hist",
              plotOutput("hist")
            ),
            tabPanel(
              "Mean", value = "mean",
              plotOutput("density_mean")
            ),
            tabPanel(
              "CV²", value = "cv2",
              plotOutput("density_cv2")
            ),
            tabPanel(
              "Fano factor", value = "fano",
              plotOutput("density_fano")
            ),
            tabPanel(
              "Entropy", value = "entropy",
              plotOutput("density_entropy")
            ),
            tabPanel(
              "Extropy", value = "extropy",
              plotOutput("density_extropy")
            ),
            tabPanel(
              "Total counts", value = "sum",
              plotOutput("density_sum")
            )
          ),
          width=4
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tags$br(),
          tags$br(),
          dataTableOutput("df")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  hist_plot<-reactive({
    req(input$prom)
    df<-long_df%>%filter(seqname==input$prom)
    if(nrow(df)!=5){
      df<-bind_rows(df, data.frame(seqname=input$prom,gate=as_factor(which(!seq(1,5)%in%df$gate)),tot=0))
    }
    ggplot(df,aes(gate,tot))+geom_col(fill="#5f9ea0")+ggtitle(input$prom)+xlab("Gates")+ylab("Counts")
  })
  density_mean<-reactive({
    value=NULL
    if(input$prom%in%list_stats$seqname){
      value=list_stats%>%filter(seqname==input$prom)%>%select(mean)%>%as.numeric()
    }
    list_stats%>%ggplot(aes(mean))+geom_density(fill="#317873")+geom_vline(xintercept=value,linetype="dashed")
  })
  density_cv2<-reactive({
    value=NULL
    if(input$prom%in%list_stats$seqname){
      value=list_stats%>%filter(seqname==input$prom)%>%select(cv2)%>%as.numeric()
    }
    list_stats%>%ggplot(aes(cv2))+geom_density(fill="#49796b")+theme_bw()+geom_vline(xintercept=value,linetype="dashed")
  })
  density_fano<-reactive({
    value=NULL
    if(input$prom%in%list_stats$seqname){
      value=list_stats%>%filter(seqname==input$prom)%>%select(fano_factor)%>%as.numeric()
    }
    list_stats%>%ggplot(aes(fano_factor))+geom_density(fill="#a3c1ad")+geom_vline(xintercept=value,linetype="dashed")
  })
  density_entropy<-reactive({
    value=NULL
    if(input$prom%in%list_stats$seqname){
      value=list_stats%>%filter(seqname==input$prom)%>%select(entropy)%>%as.numeric()
    }
    list_stats%>%ggplot(aes(entropy))+geom_density(fill="#49796b")+geom_vline(xintercept=value,linetype="dashed")
  })
  density_extropy<-reactive({
    value=NULL
    if(input$prom%in%list_stats$seqname){
      value=list_stats%>%filter(seqname==input$prom)%>%select(extropy)%>%as.numeric()
    }
    list_stats%>%ggplot(aes(extropy))+geom_density(fill="#a3c1ad")+geom_vline(xintercept=value,linetype="dashed")
  })
  density_sum<-reactive({
    value=NULL
    if(input$prom%in%list_stats$seqname){
      value=list_stats%>%filter(seqname==input$prom)%>%select(sum)%>%as.numeric()
    }
    list_stats%>%ggplot(aes(sum))+geom_density(fill="#a3c1ad")+geom_vline(xintercept=value,linetype="dashed")+xlab("Total counts")+scale_x_log10()
  })
  
      output$df <- renderDataTable(as.data.table(list_stats,  options = list(autoWidth = TRUE), filter = list(
        position = 'bottom', clear = FALSE, plain = TRUE
      )))
      output$hist<-renderPlot(hist_plot())
      output$density_mean<-renderPlot(density_mean())
      output$density_cv2<-renderPlot(density_cv2())
      output$density_fano<-renderPlot(density_fano())
      output$density_entropy<-renderPlot(density_entropy())
      output$density_extropy<-renderPlot(density_extropy())
      output$density_sum<-renderPlot(density_sum())
      
}

# Run the application 
shinyApp(ui = ui, server = server)
