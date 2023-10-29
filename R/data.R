

#' @title Whole blood reference matrix of 200 miRNAs and 6 blood cell subtypes
#' @description Reference-based cell-type fraction estimation methods rely on a
#' prior defined reference matrix. We used miRNA expression profiling data from
#' purified blood cell types from Simonas et al (2012) to construct this whole
#' blood reference matrix by using ANOVA method. It contains 200 miRNAs of 6
#' blood cell subtypes:
#' \itemize{
#'   \item CD4+ T-cells
#'   \item CD8+ T-cells
#'   \item B-cells
#'   \item NK-cells
#'   \item Monocytes
#'   \item Neutrophils
#' }
#' @format A matrix with 200 rows (miRNAs) and 6 variables (Cell Types)
#' @usage data(refMatrix_Blood0.m)
#' @docType data
#' @keywords datasets
#' @name refMatrix_Blood0.m
"refMatrix_Blood0.m"

#' @title Whole blood reference matrix of 145 miRNAs and 6 blood cell subtypes
#' @description Reference-based cell-type fraction estimation methods rely on a
#' prior defined reference matrix. We used miRNA expression profiling data from
#' purified blood cell types from Simonas et al (2012) to construct this whole
#' blood reference matrix by using DEM methods. It contains 145 miRNAs of 6
#' blood cell subtypes:
#' \itemize{
#'   \item CD4+ T-cells
#'   \item CD8+ T-cells
#'   \item B-cells
#'   \item NK-cells
#'   \item Monocytes
#'   \item Neutrophils
#' }
#' @format A matrix with 145 rows (miRNAs) and 6 variables (Cell Types)
#' @usage data(refMatrix_Blood1.m)
#' @docType data
#' @keywords datasets
#' @name refMatrix_Blood1.m
"refMatrix_Blood1.m"

#' @title Whole blood reference matrix of 140 miRNAs and 6 blood cell subtypes
#' @description Reference-based cell-type fraction estimation methods rely on a
#' prior defined reference matrix. We used miRNA expression profiling data from
#' purified blood cell types from FANTOM5 project to construct this whole
#' blood reference matrix by using DEM method. It contains 140 miRNAs of 6
#' blood cell subtypes:
#' \itemize{
#'   \item CD4+ T-cells
#'   \item CD8+ T-cells
#'   \item B-cells
#'   \item NK-cells
#'   \item Monocytes
#'   \item Neutrophils
#' }
#' @format A matrix with 140 rows (miRNAs) and 6 variables (Cell Types)
#' @usage data(refMatrix_Blood2.m)
#' @docType data
#' @keywords datasets
#' @name refMatrix_Blood2.m
"refMatrix_Blood2.m"


#' @title miRNA expression profiles (TCGA LAML)
#' @description An example data table that containing the miRNA read count from
#' bulk miRNA-seq of 12 Whole Blood samples form TCGA LAML (Acute Myeloid Leukemia)
#' cohort. Only used for demonstration purpose.
#' @format A matrix with 1807 rows (miRNAs) and 12 variables (samples)
#' @usage data(miRExpCount_LAML.m)
#' @docType data
#' @keywords datasets
#' @name miRExpCount_LAML.m
"miRExpCount_LAML.m"


#' @title Breast tissue reference matrix of 271 miRNAs and 10 cell subtypes
#' @description Reference-based cell-type fraction estimation methods rely on a
#' prior defined reference matrix. We used miRNA expression profiling data from
#' purified blood cell types from FANTOM5 to construct this whole blood reference
#' matrix. It contains 271 miRNAs of 6 blood cell subtypes:
#' \itemize{
#'   \item Epithelial
#'   \item Fibroblast
#'   \item Adipocyte
#'   \item CD4+ T-cells
#'   \item CD8+ T-cells
#'   \item B-cells
#'   \item NK-cells
#'   \item Monocytes
#'   \item Macrophage
#'   \item Neutrophils
#' }
#' @format A matrix with 271 rows (miRNAs) and 10 variables (Cell Types)
#' @usage data(refMatrix_Breast.m)
#' @docType data
#' @keywords datasets
#' @name refMatrix_Breast.m
"refMatrix_Breast.m"


#' @title miRNA expression profiles (TCGA BRCA)
#' @description An example data table that containing the miRNA read count from
#' bulk miRNA-seq of 12 Breast tissue samples form TCGA BRCA (Breast Invasive Carcinoma)
#' cohort. Only used for demonstration purpose.
#' @format A matrix with 2235 rows (miRNAs) and 12 variables (samples)
#' @usage data(miRExpCount_BRCA.m)
#' @docType data
#' @keywords datasets
#' @name miRExpCount_BRCA.m
"miRExpCount_BRCA.m"

