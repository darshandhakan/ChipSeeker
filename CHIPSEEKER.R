#!/usr/bin/env Rscript

options(java.parameters = c("-XX:+IgnoreUnrecognizedVMOptions", "-Xmx8192m"))

suppressMessages(library(ggplot2))
suppressMessages(library(ChIPseeker))
suppressMessages(library(optparse))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(annotables))
suppressMessages(library(xlsx))
suppressMessages(library(stringr))
suppressMessages(library(R.utils))

print (version)

option_list <- list(make_option(c("-n", "--projectName"),
                                default=basename(getwd()),
                                dest="projectName",
                                help="name of the project used for storing images and tables [default: name of the current directory]."),
                    make_option(c("-b", "--Beddir"),
                                dest="bed",
                                help="Path of the Bed Files Folder"),
                    make_option(c("-s", "--sample_data"),
                                dest = "sample_data",
                                help="Name of the Sample Data Folder"),
                    make_option(c("-t", "--tx"),
                                dest = "txdb",
                                help="Enter the TXDB database required for annotations"),
                    make_option(c("-u", "--human"),
                                 default=TRUE,
                                 dest="Human",
                                 help = "If the database is for humans")
)

                    

parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list,
                       description="Perform ChipSeeker Analysis",
                       epilogue="Analysis of Chip Bed Files using ChipSeeker")

opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

##Check mandetory inputs:

if ( is.null(opt$bed) ) {
  stop("--Path to the Bed Files containing folder must be provided. See script usage (--help)")
}

if ( is.null(opt$sample_data) ) {
  stop("--Name of the Sample data Table containing samples and the corresponding BED files needed shoud be provided. See script usage (--help)")
}

BedDir=opt$bed
Sample_table=opt$sample_data
TXDB=opt$txdb
HUMAN=opt$Human

Run_ChiPseeker <- function(Bed_Dir, SAMPLE_TABLE, TXDB, Human, Custom=NULL)
{
  # Load libraries
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(clusterProfiler)
  library(annotables)
  library(org.Hs.eg.db)
  library(xlsx)
  library(stringr)
  library(R.utils)
  
  print (getOption("clusterProfiler.download.method"))
  #options(download.file.method="wget")
  #R.utils::setOption("clusterProfiler.download.method", "wget")
  #print (getOption("clusterProfiler.download.method"))
  
  ## Get the Files as List
  FILES <- list.files(path = Bed_Dir, pattern = ".bed", all.files = TRUE, full.names = TRUE, recursive = TRUE)
  Sample_table <- read.csv(file = paste0(Bed_Dir,"/",SAMPLE_TABLE), sep = "\t", header = TRUE)
  FILES2 <- str_replace(string = FILES, pattern = ".*/", replacement = "")
  Sample_table <- Sample_table[match(FILES2, Sample_table$Bed_File),]
  FILES <- as.list(FILES)
  SAMPLES <- Sample_table$Sample_Name
  names(FILES) <- SAMPLES
  
  print (FILES)
  print (Sample_table)
  
  ## Load the Txdb ##
  
  if (!is.null(Custom))
  {
    txdb <- loadDb(TXDB)
  }
  else {
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  assign("txdb", value = txdb, envir = .GlobalEnv)
  print ("PEAK Annotate Step")
  
  ## Generate the Peak Annot List
  peakAnnoList <- lapply(FILES, annotatePeak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
  print (peakAnnoList)
  
  print ("Creating Directories")
  ## Create Output Directores
  
  PLOT_OUTPUT <- paste0(Bed_Dir, "/PLOTS")
  if (!dir.exists(PLOT_OUTPUT))
  {
    dir.create(PLOT_OUTPUT)
  }
  else {
    print("dir exists")
  }
  if (!dir.exists(paste0(PLOT_OUTPUT,"/INDIVIDUAL_SAMPLE_PLOTS")))
  {
    dir.create(paste0(PLOT_OUTPUT,"/INDIVIDUAL_SAMPLE_PLOTS"))
  }
  else {
    print("dir exists")
  }
  
  TABLE_OUTPUT <- paste0(Bed_Dir, "/TABLES")
  if (!dir.exists(TABLE_OUTPUT))
  {
    dir.create(TABLE_OUTPUT)
  }
  else {
    print("dir exists")
  }
  
  ## Plot the distribution of the Annotations as pie plot:
  print ("Generating Pie plots and Bar plots")
  
  COLORS <- c("dodgerblue", "indianred", "pink1", "orange", "purple", "green", "darkgreen", "dodgerblue4", "cyan")
  FTS <- c("Promoter", "5' UTR", "3' UTR", "1st Exon", "Other Exon", "1st Intron", "Other Intron", "Downstream (<=300)", "Distal Intergenic")
  names(COLORS) <- FTS
  assign("COLORS", value = COLORS, envir = .GlobalEnv)
  assign("peakAnnoList", value = peakAnnoList, envir = .GlobalEnv)
  assign("SAMPLES", value = SAMPLES, envir = .GlobalEnv)
  for (i in seq_along(SAMPLES))
  {
    # Pie Plot
    png(paste0(PLOT_OUTPUT, "/INDIVIDUAL_SAMPLE_PLOTS/", SAMPLES[i], "_Peak_Annotation_PIE_PLOT.png"), height = 12, width = 12, units = 'in', res = 300)
    pie(peakAnnoList[[SAMPLES[i]]]@annoStat[,2], labels = peakAnnoList[[SAMPLES[i]]]@annoStat[,1], col = COLORS[as.character(peakAnnoList[[SAMPLES[i]]]@annoStat[,1])])
    legend("topright", legend = as.character(peakAnnoList[[SAMPLES[i]]]@annoStat[,1]), fill = COLORS[peakAnnoList[[SAMPLES[i]]]@annoStat[,1]])
    dev.off()
    
    # Bar Plot
    PLOT <- ggplot(peakAnnoList[[SAMPLES[i]]]@annoStat, aes(x=Feature, y=Frequency, fill = Feature)) + geom_bar(stat = "identity") + scale_fill_manual(values = COLORS[as.character(peakAnnoList[[SAMPLES[i]]]@annoStat[,1])]) + theme(axis.text.x = element_text(size = 12, face = "bold", angle = 45), axis.text.y = element_text(size = 14, face = "bold"), axis.title = element_text(size = 14, face = "bold"), legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12, face = "bold"), panel.background = element_blank())
    png(paste0(PLOT_OUTPUT, "/INDIVIDUAL_SAMPLE_PLOTS/", SAMPLES[i], "_Peak_Annotation_BAR_PLOT.png"), height = 12, width = 12, units = 'in', res = 300)
    plot(PLOT)
    dev.off()
  }
  print ("Now generating the Combined Bar plot of Region wise annotation of peaks")
  
  png(paste0(PLOT_OUTPUT, "/COMBINED_REGION_DISTRIBUTION_PLOT.png"), height = 14, width = 12, units = 'in', res = 300)
  plot(plotAnnoBar(peakAnnoList))
  dev.off()
  
  print ("Plot the Distribution of TF binding relative to TSS")
  png(paste0(PLOT_OUTPUT, "/COMBINED_TF_DISTRIBUTION_PLOT.png"), height = 14, width = 12, units = 'in', res = 300)
  plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci \n relative to TSS")
  dev.off()
  
  
  ## Generate tables:
  
  #+++++++++++++++++++++++++++
  # xlsx.writeMultipleData
  #+++++++++++++++++++++++++++++
  # file : the path to the output file
  # ... : a list of data to write to the workbook
  xlsx.writeMultipleData <- function (file, ...)
  {
    require(xlsx, quietly = TRUE)
    objects <- list(...)
    fargs <- as.list(match.call(expand.dots = TRUE))
    objnames <- as.character(fargs)[-c(1, 2)]
    nobjects <- length(objects)
    for (i in 1:nobjects) {
      if (i == 1)
        write.xlsx(objects[[i]], file, sheetName = objnames[i])
      else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                      append = TRUE)
    }
  }
  
  # - titleStyle : style object to use for title
  xlsx.addTitle<-function(sheet, rowIndex, title, titleStyle){
    rows <-createRow(sheet,rowIndex=rowIndex)
    sheetTitle <-createCell(rows, colIndex=1)
    setCellValue(sheetTitle[[1,1]], title)
    setCellStyle(sheetTitle[[1,1]], titleStyle)
  }
  
  
  #############################################
  print ("Also perform GO Enrichment Analysis")
  
  for (i in seq_along(SAMPLES))
  {
    annot <- as.data.frame(peakAnnoList[[SAMPLES[i]]]@anno)
    ## Get the Gene Symbols
    if (Human == TRUE)
    {
      entrezids <- unique(annot$geneId)
      entrez2gene <- grch38 %>% filter(entrez %in% entrezids) %>% dplyr::select(entrez, symbol)
      m <- match(annot$geneId, entrez2gene$entrez)
      annot <- cbind(annot[,1:14], geneSymbol=entrez2gene$symbol[m])
    }
    ## Now add tor the excel sheet
    wb<-createWorkbook(type="xlsx")
    
    # Title and sub title styles
    TITLE_STYLE <- CellStyle(wb)+ Font(wb,  heightInPoints=16,
                                       color="dodgerblue4", isBold=TRUE, underline=1)
    SUB_TITLE_STYLE <- CellStyle(wb) +
      Font(wb,  heightInPoints=14,
           isItalic=TRUE, isBold=FALSE)
    # Styles for the data table row/column names
    TABLE_ROWNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE)
    TABLE_COLNAMES_STYLE <- CellStyle(wb) + Font(wb, isBold=TRUE) +
      Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
      Border(color="black", position=c("TOP", "BOTTOM"),
             pen=c("BORDER_THIN", "BORDER_THICK"))
    
    sheet <- createSheet(wb, sheetName = paste0(SAMPLES[i], "_PEAK_ANNOTAITON"))
    xlsx.addTitle(sheet, rowIndex=1, title=paste0(SAMPLES[i], "Annotations to each of the Peak Calls Detected (TSS Sites)"), titleStyle = TITLE_STYLE)
    addDataFrame(annot, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE, rownamesStyle = TABLE_ROWNAMES_STYLE)
    saveWorkbook(wb, paste0(TABLE_OUTPUT, "/", SAMPLES[i], "_TSS_PEAK_ANNOTATION_TABLE.xlsx"))
    
    ## Generate the Peaks Coverage plots:
    
    PEAK_FILE <- readPeakFile(peakfile = FILES[[i]])
    png(paste0(PLOT_OUTPUT, "/INDIVIDUAL_SAMPLE_PLOTS/", SAMPLES[i], "_PEAK_COVERAGE_PLOT.png"), height = 12, width = 12, units = 'in', res = 300)
    plot(covplot(PEAK_FILE, weightCol="V5", chrs = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")))
    dev.off()
    
    ## Run GO Enrichment Analysis:
    if (Human == TRUE)
    {
      entrezids <- annot$geneId %>% as.character() %>% unique()
      ego <- enrichGO(gene = entrezids, keyType = "ENTREZID", OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
      cluster_summary <- data.frame(ego)
      if (nrow(cluster_summary) > 0)
      {
        wb2<-createWorkbook(type="xlsx")
        TITLE_STYLE2 <- CellStyle(wb2)+ Font(wb2,  heightInPoints=16,
                                             color="dodgerblue4", isBold=TRUE, underline=1)
        SUB_TITLE_STYLE2 <- CellStyle(wb2) +
          Font(wb2,  heightInPoints=14,
               isItalic=TRUE, isBold=FALSE)
        # Styles for the data table row/column names
        TABLE_ROWNAMES_STYLE2 <- CellStyle(wb2) + Font(wb2, isBold=TRUE)
        TABLE_COLNAMES_STYLE2 <- CellStyle(wb2) + Font(wb2, isBold=TRUE) +
          Alignment(wrapText=TRUE, horizontal="ALIGN_CENTER") +
          Border(color="black", position=c("TOP", "BOTTOM"),
                 pen=c("BORDER_THIN", "BORDER_THICK"))
        
        sheet <- createSheet(wb2, sheetName = paste0(SAMPLES[i], "_GO_ENRICHMENT"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0(SAMPLES[i], " GO Enrichment of the Peaks Identified"), titleStyle = TITLE_STYLE2)
        addDataFrame(cluster_summary, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
        
        
        # Plot GO Enrichment:
        png(paste0(PLOT_OUTPUT, "/INDIVIDUAL_SAMPLE_PLOTS/", SAMPLES[i],"_GO_ENRICHMENT_PLOT.png"), height = 12, width = 14, units = 'in', res = 300)
        plot(dotplot(ego, showCategory=30))
        dev.off()
      }
      
      print ("Now generating the KEGG enrichment")
      
      # Plot the KEGG Enrichment:
      ekegg <- enrichKEGG(gene = entrezids, organism = 'hsa', pvalueCutoff = 1, keyType="kegg")
      if (nrow(ekegg@result) > 0)
      {
        png(paste0(PLOT_OUTPUT, "/INDIVIDUAL_SAMPLE_PLOTS/", SAMPLES[i],"_KEGG_ENRICHMENT_PLOT.png"), height = 12, width = 14, units = 'in', res = 300)
        plot(dotplot(ekegg))
        dev.off()
        sheet <- createSheet(wb2, sheetName = paste0(SAMPLES[i], "_KEGG_ENRICHMENT"))
        xlsx.addTitle(sheet, rowIndex=1, title=paste0(SAMPLES[i], " KEGG Enrichment of the Peaks Identified"), titleStyle = TITLE_STYLE2)
        addDataFrame(ekegg@result, sheet, startRow=3, startColumn=1, colnamesStyle = TABLE_COLNAMES_STYLE2, rownamesStyle = TABLE_ROWNAMES_STYLE2)
        saveWorkbook(wb2, paste0(TABLE_OUTPUT, "/",  SAMPLES[i], "_FUNCTIONAL_ENRICHMENT_TABLE.xlsx"))
      }
    }
  }
  
  ## Generate the Promoiter TSS plot for all samples:
  
  promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
  print ("Promoters fetched")
  tagMatrixList <- lapply(as.list(FILES), getTagMatrix, windows=promoter)
  print ("tagmatrix generated")
  assign("tagMatrixList", value = tagMatrixList, envir = .GlobalEnv)
  png(paste0(PLOT_OUTPUT, "/COMBINED_TSS_PROMOTER_PEAK_PLOTS.png"), height = 12, width = 14, units = 'in', res = 300)
  plot(plotAvgProf(tagMatrixList, xlim=c(-1000, 1000), conf=0.95,resample=500, facet="row"))
  dev.off()
}
Run_ChiPseeker(Bed_Dir = BedDir, SAMPLE_TABLE = Sample_table, TXDB = TXDB, Human = HUMAN)

