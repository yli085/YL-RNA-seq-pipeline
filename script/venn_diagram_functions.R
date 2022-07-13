## venn.diagram customized function
## Input parameters
## s1, s2: the first and the second input gene lists
## s1_name, s2_name: the names to show on the Venn diagram
## file_dir: the file directory for the output figure

library(writexl)
library(VennDiagram)

venn.diagram.two <- function(s1, s2, s1_name, s2_name, file_dir) {
  
  vennFilename <- paste0(file_dir, s1_name, "&", s2_name, ".png")
  vennDiagramColors <- c("dodgerblue", "darkorange1")
  dataForVennDiagram <- list(s1, s2)
  
  s1_rename <- paste0(s1_name, " (", length(s1), ")")
  s2_rename <- paste0(s2_name, " (", length(s2), ")")
  category.names <- c(s1_rename, s2_rename)
  
  venn.diagram(dataForVennDiagram,
               category.names = category.names,
               filename = vennFilename,
               imagetype = 'png',
               disable.logging = T, 
               scaled = T,
               col = vennDiagramColors,
               fill = vennDiagramColors,
               lwd = 2,
               alpha = 0.25,
               fontfamily = 'sans',
               margin = 0.2,
               cat.dist = 0.05,
               cat.default.pos = "outer",
               cat.pos = c(-25, 25),
               cat.col = vennDiagramColors,
               cat.fontface = 'bold',
               cat.fontfamily = 'sans',
               cat.cex = 1,
               cex = 1)
  
  # the intersection 
  intersect <- data.frame(Gene = intersect(s1, s2))
  excelFilename <- paste0(file_dir, s1_name, "&", s2_name, ".xlsx")
  write_xlsx(intersect, excelFilename)
  
}



# the input is 3 sets
venn.diagram.three <- function(s1, s2, s3, s1_name, s2_name, s3_name, file_dir){
  
  vennFilename <- paste0(file_dir, s1_name, "&", s2_name, "&", s3_name, ".png")
  vennDiagramColors <- c('#EA4335', '#FBBC05', '#4285F4')
  dataForVennDiagram <- list(s1, s2, s3)
  
  s1_rename <- paste0(s1_name, " (", length(s1), ")")
  s2_rename <- paste0(s2_name, " (", length(s2), ")")
  s3_rename <- paste0(s3_name, " (", length(s3), ")")
  category.names <- c(s1_rename, s2_rename, s3_rename)
  
  venn.diagram(
    dataForVennDiagram,
    category.names = category.names,
      
    filename = vennFilename,
    disable.logging = T, 
    
    # Output features
    imagetype="png" ,
    margin=0.2,
    
    # Circles
    lwd = 2,
    col = vennDiagramColors,
    fill = vennDiagramColors,
    alpha=0.25,
    
    # Numbers
    cex = 1,
    fontfamily = "sans",
    
    # Set names
    cat.col = vennDiagramColors,
    cat.cex = 1,
    cat.fontface = "bold",
    cat.fontfamily = "sans",
    # cat.pos = c(-35, 35, 120),
    cat.dist = 0.07,

  )
  
  # the intersection 
  intersect <- data.frame(Gene = intersect(intersect(s1, s2), s3))
  excelFilename <- paste0(file_dir, s1_name, "&", s2_name, "&", s3_name, ".xlsx")
  write_xlsx(intersect, excelFilename)
  
}
