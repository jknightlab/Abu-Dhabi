list.of.packages <- c("ggplot2", "shiny")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

suppressMessages(library(shiny))
suppressMessages(library(ggplot2))
suppressMessages(library(data.table))

# eQTL results
load("Peak_eQTL.RData")
all.eqtl.tested <- readRDS("eQTL_results_final_lm_cis.rds")

# Interactions
comb.int <- readRDS("Comb.int.rds")

comb.int <- comb.int[comb.int$Gene %in% peak.eQTL$Gene[peak.eQTL$qvalue < 0.05], ]
comb.int$DiseaseFDR <- p.adjust(comb.int$Disease, method="fdr")
comb.int$GlucoseFDR <- p.adjust(comb.int$Glucose, method="fdr")
comb.int$HbA1cFDR <- p.adjust(comb.int$HbA1c, method="fdr")

# trans.eqtl <- read.delim("Significant_trans_eQTL_Bonferroni.txt", stringsAsFactors=F)
# trans.eqtl <- trans.eqtl[order(trans.eqtl$eQTL_pval), ]
# trans.eqtl <- trans.eqtl[!duplicated(trans.eqtl$Gene), ]
# trans.geno <- readRDS("trans.geno.rds")

# Differential Expression
load("DE.Rdata")

# GWAS
# GWAS <- read.table("AD_no_parents_no_outliers_T2D.assoc", header=TRUE)
snp.pos <- data.frame(fread("SNP_pos.txt"))

genelist <- read.delim("Gene_info.txt", header=T, stringsAsFactors=F, sep="")

ld <- read.delim("pairwise_LD.ld.gz", sep="", stringsAsFactors=F)

# Define UI ----
ui <- fluidPage(

  # App title
  titlePanel("eQTL in the Abu Dhabi cohort"),

  # Selectors as first row
  fluidRow(
    column(4,
           align="center",
           selectizeInput("geneSelector", choices=NULL, selected=NULL, label="Please select a gene name", multiple=FALSE)),

    column(8,
           align="center",
           radioButtons("factor",
                        label = "Please select an environmental factor",
                        choices = list("Disease Status" = "Disease",
                                       "HbA1c" = "HbA1c",
                                       "Glucose" = "Glucose"),
                        selected=character(0),
                        inline=TRUE))),

  #Text corresponding to eQTL and interaction p value as next row
  fluidRow(
    column(4,
           align="center",
           textOutput("eQTL_pval")),

    column(8,
           offset = 4,
           align="center",
           textOutput("Int_pval"))),

  # Links to type 2 diabetes GWAS results

  # uiOutput("tab"),

  # eQTl figure, interaction figures
  fluidRow(
    column(4,
           plotOutput("fig1")),

    conditionalPanel(
      condition = "input$geneSelector != 'NULL'",

      column(4,
             plotOutput("fig2")),

      column(4,
             plotOutput("fig3")))),

  #Add a download button under each plot
  fluidRow(
    column(4,
           align="center",
           downloadButton("down_eQTL", label="Download the plot")),

    column(4,
           align="center",
           downloadButton("down_int1", label="Download the plot")),

    column(4,
           align="center",
           downloadButton("down_int2", label="Download the plot"))),

  # Local association plot
  fluidRow(
    column(12, plotOutput("LAP"))
  ),

  # Differential expression and eQTL table
  fluidRow(
    column(4,
           align="center",
           textOutput("DE_pval"))),

  fluidRow(
    column(4,
           plotOutput("fig4")),

    column(4,
           dataTableOutput("table1"))

    ),


  fluidRow(
    column(4,
           align="center",
           downloadButton("down_DE", label="Download the plot"))

  )

)



# Define server logic ----
server <- function(input, output, session) {

  updateSelectizeInput(session, "geneSelector",
                       choices = peak.eQTL,
                       selected="CUTALP",
                       options = list(labelField='Gene', searchField='Gene',
                                      valueField='Gene',
                                      render = I("{ option: function(item, escape) {return '<div><strong>' + escape(item.Gene) + '</span></div>';} }")),
                       server = TRUE)


  ###############################################################################
  # Text for displaying p values
  ###############################################################################

  #Text for eQTL p value
  output$eQTL_pval <- renderText({
    paste("eQTL p value ", signif(peak.eQTL[peak.eQTL$Gene==input$geneSelector, "eQTL_pval"], 2),
          "\n q value ", signif(peak.eQTL[peak.eQTL$Gene == input$geneSelector, "qvalue"], 2), sep="")
  })


  #Function for plotting eQTL figure
  eqtl.plot.side <- function(gene, snp){

    eqtl <- data.frame("Expression"=exprs[gene, ],
                       "Genotype"=as.factor(geno.peak[, snp]))
    eqtl <- eqtl[complete.cases(eqtl), ]

    ggplot(eqtl, aes(Genotype, Expression, colour=Genotype)) +
      geom_boxplot() +
      geom_point(position=position_jitter(width=0.2)) +
      xlab(label=snp) +
      ylab(label=paste("Adjusted", gene, "Expression")) +
      theme_bw() +
      ggtitle("eQTL (Linear Model)")
  }

  #Text for eQTL int value
  output$Int_pval <- renderText({
    paste("eQTL interaction p value ", signif(comb.int[comb.int$Gene==input$geneSelector, input$factor], 2), sep="")
  })

  # url <- a("Type 2 Diabetes GWAS associations",
           # href=paste0("http://www.type2diabetesgenetics.org/gene/geneInfo/", input$geneSelector))
  # output$tab <- renderUI({
    # tags$a("URL link:", url)
  # })

  # Functions for plotting eQTL interactions
  int.boxplot <- function(gene, snp, int){

    eqtl <- data.frame("Expression"=exprs[gene, ],
                       "Genotype"=as.factor(geno.peak[, snp]),
                       "Modulator"=info[, int])
    eqtl <- eqtl[complete.cases(eqtl), ]

    if(class(eqtl$Modulator) == "factor"){
      eqtl <- eqtl[eqtl$Modulator %in% c("Control", "T2DM"),]
      ggplot(eqtl, aes(Genotype, Expression, colour=Modulator)) +
        facet_grid(~ Modulator) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(position=position_jitter(width=0.2)) +
        xlab(label=snp) +
        ylab(label=paste("Adjusted", gene, "Expression")) +
        theme_bw() +
        ggtitle("eQTL Interaction (Linear Model)")


    # ggplot(eqtl, aes(Genotype, Expression, colour=Modulator)) +
    #   geom_boxplot(outlier.shape=NA) +
    #   geom_point(aes(group=Modulator), position=position_jitterdodge()) +
    #   xlab(label=snp) +
    #   ylab(label=paste("Adjusted", gene, "Expression")) +
    #   theme_bw() +
    #   ggtitle("eQTL Interaction (Linear Model)")
    } else {

    # if variable is continuous
      eqtl <- eqtl[order(eqtl$Modulator), ]
      eqtl$Rank <- 1:nrow(eqtl)
      ggplot(eqtl, aes(Genotype, Expression, colour=Modulator)) +
        scale_color_gradient2(midpoint=mean(range(eqtl$Modulator, na.rm=T)),
                                       low="blue", mid="gold", high="red",
                              limits=range(eqtl$Modulator, na.rm=T)) +
        geom_boxplot(outlier.shape=NA) +
        geom_point(aes(group=Modulator), position=position_jitter(width=0.2)) +
        xlab(label=int) +
        ylab(label=paste("Adjusted", gene, "Expression")) +
        theme_bw() +
        ggtitle("eQTL Interaction (Linear Model)")
    }

  }

  int.dotplot <- function(gene, snp, int){

    eqtl <- data.frame("Expression"=exprs[gene, ],
                       "Genotype"=as.factor(geno.peak[, snp]),
                       "Modulator"=info[, int])
    eqtl <- eqtl[complete.cases(eqtl), ]
    if(int == "Disease"){
      eqtl <- eqtl[eqtl$Modulator %in% c("Control", "T2DM"), ]
    }

    ggplot(eqtl, aes(Modulator, Expression, colour=Genotype, group=Genotype)) +
      geom_point(position=position_jitter(width=0.2)) +
      geom_smooth(method="lm", se=FALSE, na.rm=TRUE) +
      xlab(label=int) +
      ylab(label=paste("Adjusted", gene, "Expression")) +
      theme_bw() +
      ggtitle("eQTL Interaction (Linear Model)")
  }


  #Text for DE p value
  output$DE_pval <- renderText({
    paste("Differential expression FDR ",
          signif(comb.res[comb.res$Gene==input$geneSelector, input$factor], 2), sep="")
  })


  # Function for plotting DE figure
  de.plot <- function(gene, variable){

    if(variable == "Disease"){
      my.factor <- dds.disease
      exprn <- dis.exprn
      disease <- dds.disease
    } else if(variable == "HbA1c") {
      exprn <- hba1c.exprn
      my.factor <- dds.hb1ac
      disease <- dds.disease[!is.na(my.factor)]
      my.factor <- my.factor[!is.na(my.factor)]
    } else if(variable == "Glucose") {
      exprn <- glu.exprn
      my.factor <- dds.glucose
      disease <- dds.disease[!is.na(my.factor)]
      my.factor <- my.factor[!is.na(my.factor)]
    }

    exprn <- data.frame("Factor"=my.factor,
                        "Expression"=exprn[gene, ],
                        "Disease"=disease)

    if(class(exprn$Factor) == "factor"){
      # if variable is a factor
      ggplot(exprn, aes(Factor, Expression, colour=Disease)) +
        geom_boxplot() +
        geom_point(position=position_jitter(width=0.2)) +
        xlab(label=variable) +
        ylab(label=paste("Adjusted", gene, "Expression")) +
        theme_bw() +
        # scale_y_log10() +
        ggtitle("Differential Gene Expression (DESeq2)")
    } else {

      # if variable is continuous
      ggplot(exprn, aes(Factor, Expression)) +
        geom_smooth(method="lm") +
        geom_point(aes(colour=Disease)) +
        xlab(label=variable) +
        ylab(label=paste("Adjusted", gene, "Expression (DESeq2)")) +
        theme_bw()
    }
  }

  make.fancy.locus.plot <- function(gene){
    results <- subset(all.eqtl.tested, Gene == gene)
    if(nrow(results) == 0){
      print("There is no significant association")
    } else {
    locus <- data.frame(SNP=results$SNP,
                        CHR=snp.pos$Chr[match(results$SNP, snp.pos$Edit)],
                        POS=snp.pos$BP[match(results$SNP, snp.pos$Edit)],
                        PVAL=results$eQTL_pval)
    locus$SNP <- snp.pos$Full[match(locus$SNP, snp.pos$Edit)]
    snp <- as.character(locus$SNP[which.min(results$eQTL_pval)])

    snp.ld <- subset(ld, SNP_A == snp)
    locus$RSQR <- snp.ld$R2[match(locus$SNP, snp.ld$SNP_B)]
    locus$RSQR[which.min(locus$P)] <- 1
    hit <- locus[which.min(locus$PVAL), ]
    chr <- hit$CHR

    locus <- locus[complete.cases(locus), ]

    #
    # genes in the region
    #
    min.pos <- min(locus$POS, na.rm=T)
    max.pos <- max(locus$POS, na.rm=T)
    genes.in.locus <- subset(genelist,
                             (genelist$START > min.pos & genelist$START < max.pos & genelist$CHR == chr) |
                               (genelist$STOP > min.pos & genelist$STOP < max.pos & genelist$CHR == chr) )
    genes.in.locus <- rbind(genes.in.locus, genelist[genelist$GENE == gene, ])

    #
    # size of the region
    #
    min.pos <- min(c(locus$POS,
                     genelist$START[genelist$GENE == gene],
                     genelist$STOP[genelist$GENE == gene]), na.rm=T) - 1000
    max.pos <- max(c(locus$POS, genelist$START[genelist$GENE == gene],
                     genelist$STOP[genelist$GENE == gene]), na.rm=T) + 1000
    size.pos <- max.pos - min.pos
    center.pos <- min.pos + ( size.pos / 2 )
    center.100kb.pos <- round(center.pos / 100000) * 100000
    offset.100kb.pos <- round((size.pos/3) / 100000) * 100000

    #
    # range of y-axis
    #
    # this dedicates 33% of the yaxis to the genes, labels, recomb rate
    range <- ceiling(-log10(hit$PVAL)) + 1
    offset <- ( range * 4 / 3 ) - range
    big.range <- range + offset

    ystart.gene <- - offset

    #
    # SNPS
    #
    markers.in.strong.ld <- subset(locus, (row.names(locus) != snp &
                                             locus$RSQR >= 0.8))
    markers.in.moderate.ld <- subset(locus, (row.names(locus) != snp &
                                               locus$RSQR >= 0.5 &
                                               locus$RSQR < 0.8))
    markers.in.weak.ld <- subset(locus, (row.names(locus) != snp &
                                           locus$RSQR >= 0.2 &
                                           locus$RSQR < 0.5))
    markers.not.in.ld <- subset(locus, (row.names(locus) != snp &
                                          locus$RSQR<0.2))


    par(mar=c(4, 4, 3, 4))

    #
    # start plot
    #
    plot(1, type="n", xlim=c(min.pos, max.pos), ylim=c(-offset, range),
         xlab="", ylab="", main=gene, axes=F)


    #
    # axes, titles and legends
    #
    mtext(paste("Chromosome", chr, "position (kb)", sep=" "), side=1, line=2.5)
    axis(1, at=c(center.100kb.pos - offset.100kb.pos,
                 center.100kb.pos,
                 center.100kb.pos + offset.100kb.pos),
         labels=c((center.100kb.pos - offset.100kb.pos) / 1000,
                  center.100kb.pos / 1000,
                  (center.100kb.pos + offset.100kb.pos) / 1000), las=1)

    axis(2, at=seq(0, range, 5), labels=seq(0, range, 5), las=1)
    mtext("Observed (-logP)", side=2, at=(range/2), line=2)


    box()
    lines(c(min.pos, max.pos), c(0,0), lty="dotted", lwd=1, col="black")

    #
    # plot the genotyped markers
    #
    points(markers.not.in.ld$POS, -(log10(markers.not.in.ld$PVAL)), pch=23,
           cex=1.0, bg="white")
    points(markers.in.weak.ld$POS, -(log10(markers.in.weak.ld$PVAL)), pch=23,
           cex=1.25, bg="yellow")
    points(markers.in.moderate.ld$POS, -(log10(markers.in.moderate.ld$PVAL)),
           pch=23, cex=1.25, bg="orange")
    points(markers.in.strong.ld$POS, -(log10(markers.in.strong.ld$PVAL)), pch=23,
           cex=1.25, bg="red")


    #
    # this is the peak SNP
    #
    points(hit$POS, -(log10(hit$PVAL)), pch=23, cex=2.5, bg="red")
    text(hit$POS, -(log10(hit$PVAL)),
           labels=c(paste(hit$SNP, "\n", "P=", signif(hit$PVAL, 3), sep="")),
           pos=4, offset=2)

    #
    # plot the genes
    #
    for ( i in 1:nrow(genes.in.locus) ) {
      if ( genes.in.locus[i, ]$STRAND == "+" ) {
        arrows(max(genes.in.locus[i, ]$START, min.pos), -offset + (big.range / 8),
               min(genes.in.locus[i, ]$STOP, max.pos), -offset + (big.range / 8), length=0.05,
               lwd=2, code=2, lty="solid", col="darkgreen")
        if ( ! is.na(genes.in.locus[i, ]$GENE) ) {
          text(genes.in.locus[i, ]$START + (genes.in.locus[i, ]$SIZE / 2),
               -offset + ( big.range / 20 )+ (big.range / 8),
               labels=genes.in.locus[i,]$GENE,
               cex=0.8, font=3, col=ifelse(genes.in.locus[i, ]$GENE == gene, "red", "black"))
        }
      }
      else if ( genes.in.locus[i, ]$STRAND == "-" ) {
        arrows(max(genes.in.locus[i, ]$START, min.pos), -offset,
               min(genes.in.locus[i, ]$STOP, max.pos), -offset, length=0.05,
               lwd=2, code=1, lty="solid", col="darkgreen")
        if ( ! is.na(genes.in.locus[i, ]$GENE) ) {
          text(genes.in.locus[i, ]$START + (genes.in.locus[i, ]$SIZE / 2),
               -offset + ( big.range / 20 ), labels=genes.in.locus[i,]$GENE,
              cex=0.8, font=3, col=ifelse(genes.in.locus[i, ]$GENE == gene, "red", "black"))
        }
      }

    }
    }
  }


  ###############################################################################
  # Plotting and downloading figures
  ###############################################################################

  # Plot eQTL figure
  output$fig1 <- renderPlot({
    req(input$geneSelector)
    eqtl.plot.side(input$geneSelector, as.character(peak.eQTL[peak.eQTL$Gene==input$geneSelector, "SNP"]))
  })

  # Download the eQTL figure as a pdf
  output$down_eQTL <- downloadHandler(
    filename = function(){
      paste(input$geneSelector, "_eQTL", ".pdf", sep="")
      },
    content = function(file){
      pdf(file, useDingbats=FALSE)
      print(eqtl.plot.side(input$geneSelector, as.character(peak.eQTL[peak.eQTL$Gene==input$geneSelector, "SNP"])))
      dev.off()
    }
)

  # Plot interaction figures
  output$fig2 <- renderPlot({
    req(input$geneSelector, input$factor)
    int.boxplot(input$geneSelector, as.character(peak.eQTL[peak.eQTL$Gene==input$geneSelector, "SNP"]), input$factor)
  })

  output$fig3 <- renderPlot({
    req(input$geneSelector, input$factor)
    int.dotplot(input$geneSelector, as.character(peak.eQTL[peak.eQTL$Gene==input$geneSelector, "SNP"]), input$factor)
  })

  # Download first interaction figure
  output$down_int1 <- downloadHandler(
    filename = paste(input$geneSelector, input$factor, "_int1_eQTL", ".pdf", sep=""),
    content = function(file){
      pdf(file, useDingbats=FALSE)
      print(int.boxplot(input$geneSelector, as.character(peak.eQTL[peak.eQTL$Gene==input$geneSelector, "SNP"]), input$factor))
      dev.off()
    })

  #Plot DE figure
  output$fig4 <- renderPlot({
    req(input$geneSelector, input$factor)
    de.plot(input$geneSelector, input$factor)
  })

  #Download the DE figure as a pdf
  output$down_DE <- downloadHandler(
    filename = paste(input$geneSelector, input$factor, "_DE", ".pdf", sep=""),
    content = function(file){
      pdf(file, useDingbats=FALSE)
      print(de.plot(input$geneSelector, input$factor))
      dev.off()
    })


  # eQTL results table
  output$table1 <- renderDataTable(subset(all.eqtl, Gene == input$geneSelector),
                                   options=list(pageLength=10))

  output$LAP <- renderPlot({
    req(input$geneSelector)
    make.fancy.locus.plot(input$geneSelector)
  })

}


# Run the app ----
shinyApp(ui = ui, server = server)

