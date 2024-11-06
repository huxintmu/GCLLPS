#' Calculate Gastric Cancer Lifespan LncRNA Prognostic Signature (GCLLPS)
#'
#' This function calculates the GCLLPS based on a given gene expression matrix.
#' The expression matrix should have rows as samples and columns as genes.
#' If the matrix is missing required genes, an error will be thrown.
#'
#' @param my_expr A numeric matrix with rows as samples and columns as genes.
#' @return A data frame with sample names and their GCLLPS scores.
#' @export
predict_GCLLPS <- function(my_expr) {
  required_genes <- c("AC007364.1","AC007557.2","AC092620.2",
                      "CDC42-IT1","COLCA1","CTA-29F11.1",
                      "CTC-338M12.4","CTD-2619J13.13","DLG3-AS1",
                      "FOXP1-IT1","GS1-166A23.2","HHIP-AS1",
                      "HOTAIR","LINC00094","LINC00152",
                      "LINC00520","LINC00847","LINC00973",
                      "LINC01089","LINC01094","MIR4435-1HG",
                      "NFE4","RP11-190A12.8","RP11-212P7.2",
                      "RP11-250B2.6","RP11-255P5.2","RP11-379F4.6",
                      "RP11-440I14.2","RP11-465B22.8","RP11-524D16__A.3",
                      "RP11-549J18.1","RP11-792A8.4","RP3-428L16.2",
                      "SERPINB9P1","SNHG22","TUG1","UBAC2-AS1","ZEB1-AS1")
  missing_genes <- required_genes[!required_genes %in% colnames(my_expr)]
  if (length(missing_genes) > 0) {
    stop("Error: The following genes are missing in the expression matrix and are required for GCLLPS calculation: ",
         paste(missing_genes, collapse = ", "), ". Please provide a complete expression matrix.")
  }
  load(system.file("data/fit.rdata", package = "GCLLPS"))
  GCLLPS <- cbind(rownames(my_expr), my_GCLLPS=predict(fit, newdata = my_expr)$predicted)
  return(GCLLPS)
}
