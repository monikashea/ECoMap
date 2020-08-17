# ECoMap
Code and data demonstrating ECoMap approach for mapping the Wisconsin Tension Zone. The data and code accompany the article: Shea, M.E., M.K. Clayton, P.A. Townsend, S. Berg, H. Elza, and D.J. Mladenoff. 2020. Identifying ecotone location using the co-occurrence property. Journal of Vegetation Science. https://doi.org/10.1111/jvs.12929

Files include:  
1. TZ_ECoMap.R: The main R script implementing the ECoMap approach for mapping the Wisconsin Tension Zone.
2. TZ_ECoMap_movezones.R: The R script exploring what happens when alternate initial boundary locations are used. Includes code for making Fig. 4. 
3. TZ_ECoMap_TablesFigs.R: The R script used to produce tables and figures, minus Fig. 4 and the supplemental figures. Includes code fully implementing ECoMap approach.
4. TZ_ECoMap_ISAcellsize: The R script exploring using different number of samples and different cell sizes for indicator species analysis. This is from the supplement Appendix S2.
5. TZ_ECoMap_bandwidth: The R script exploring the effect of using different bandwidths for density mapping. This is from the supplement Appendix S2.
6. TZ_ECoMap_RatioCV_loess: The R script exploring the effect of using different histogram cutoff values, bin sizes, and loess model span values when determining the critical point relative ratio value. This is from the supplement Appendix S2.
7. all_trees_pred_12Aug2019: Witness tree dataset, differentiated. See Appendix S1 for information on differentiation. The raw dataset may be requested from the authors. Metadata for dataset is in "metadata.docx" file.
8. 6mi_Grid_ecoregions_area.csv: Area of each 6mi grid cell within each Forest Service ecoregion. Metadata for dataset is in "metadata.docx" file.
9. 6mi_Grid: Folder containing shapefile and associated files for the 6mi_Grid used in the analysis. Metadata for dataset is in "metadata.docx" file.
10. Output: Empty folder that us referred to in code.
11. Supplement: Folder containing subfolders (2km_Grid, 4km_Grid, 12mi_Grid, 18mi_Grid), each containing shapefiles and associated files for differently-sized grids, used in the analysis described in Appendix S2. Metadata for datasets is in "metadata.docx" file.
12. WI_HARN_mask: Folder containing shapfile and associated files for the Wisconsin boundary, using the NAD83/Wisconsin Transverse Mercator map projection (EPSG:3070). Metadata for dataset is in "metadata.docx" file.
13. GLOUserGuide_v3.pdf: Document providing background on the raw dataset used to create all_trees_pred_12Aug2019.csv. The dataset has been undated since this documentation was published, but it still mostly applies. For more detailed information on updates or anything else about the dataset, please contact the authors.
14. metadata.docx: Metadata document for the datasets in this repository.

All datasets are licensed under the Creative Commons Attribution 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/ or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA. 

If there are any questions or requests for data/code that is not available here, please contact the authors.
