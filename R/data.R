

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


#' @title Reference matrices for 9 different tissues
#' @description We used miRNA expression profiling data from purified cells
#' obtained from FANTOM5 to construct reference for solid tissues including the
#' Breast, Lung, Liver, Kidney, Prostate, Esophagus, Respiratory_tract, Vessel
#' and Skin. For breast tissue, it contains 10 cell subtypes dominated by the
#' epithelial, fibroblast, adipocyte and immune-cell component of seven cell types.
#' whereas for other tissues, the reference mainly contains 9 cell types that
#' dominated by the epithelial, fibroblast and immune-cell components.
#' \itemize{
#'   \item Breast
#'   \item Lung
#'   \item Liver
#'   \item Kidney
#'   \item Prostate
#'   \item Esophagus
#'   \item Respiratory_tract
#'   \item Vessel
#'   \item Skin
#' }
#' @format A list with 9 different reference matrices
#' @usage data(SolidTissueMatrix)
#' @docType data
#' @keywords datasets
#' @name SolidTissueMatrix
"SolidTissueMatrix"


#' @title Generic reference matrices for solid tissues
#' @description Generic reference matrices consists of two references used for the
#' Hierarchical mode of DeconmiR (hDeconmiR), which is used for cell type deconvolution
#' of any solid tissue. The first is the reference centroids for epithelial cells,
#' fibroblasts and total immune cells, and the second is the reference for 7 subtypes of
#' immune cells, such as B-cells, NK cells, monocytes and etc.
#' \itemize{
#'   \item ref1
#'   \item ref2
#' }
#' @format A list with 2 different reference matrices
#' @usage data(GenericMatrix)
#' @docType data
#' @keywords datasets
#' @name GenericMatrix
"GenericMatrix"

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


