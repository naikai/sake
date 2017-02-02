# Gene List Filtering
Yu-Jui Ho and Toby Aicher, M Hammell Lab  
`r Sys.Date()`  

Data clustering analysis can be strongly affected by the input gene lists used for calculating distances between samples. Ideally, one would prefer to use a gene list with large variations across samples, and where these variations reflect biological differences that separate cell types. Several options are provided in SAKE to allow for filtering of the input gene list before proceeding to data clustering.  We also provide several visualization tools to evalutae the effects of gene list filtering on the resulting data including: preliminary correlation analysis between all samples, gene network visualization tools, and simple pairwise scatter plots between samples. 


## Filtering

We provide three different methods to filter genes from your dataset before downstream analysis. 

- **Whole transcriptome:** select whole transcriptome if you want to use all expressed genes for downstream analysis. 

- **Top Ranks:** you can rank and filter genes based on simple statistics calculated for each gene across all samples. There are currently four ways of ranking genes: mean, median, median absolute deviation (MAD), and variance. Using the top ranked genes by "mean" or "median" equates to using the most abundant genes for downstream analysis. Using the top ranked genes by "variance" equates to using the genes that change the most across samples, but this gene list can often be dominated by noise from low abundance genes.  In our experience, we find using gene lists with the top ranked `Median Absolute Deviation (MAD)` often gives the most reproducible results, as these genes have relatively high abundance and also show differential abundance between samples.  We suggest using the Top **5000 - 8000** MAD genes for `single-cell RNA-Seq` data. 

- **Upload gene list:** You also have the option to upload a text file which contains your pre-determined genes of interest. This is typically done when you have a list of marker genes used to test a specific hypothesis and want to see how samples cluster accordingly. 

An example gene list file should look like [this](../inst/extdata/genefile/genelist.EMT.txt):

Gene        | 
------------| -------------
AHNAK       |
BMP1        |
CALD1       |
CAMK2N1     |
CDH2        |
COL1A2      |
COL3A1      |
COL5A2      |
FN1         |

                     
The first row should be a character string **Gene**. The following rows should be the names/IDs of your gene of interest. 
                
## Sample Correlation and Gene Networks

We provide two tools to visualize correlations in your data and observe how filtering affects these correlations.  

### Sample Correlation 

Sample correlation creates a matrix with your samples as the rows and columns. The tool calculates the `Pearson correlation` between each of your samples using the genes you selected above. More similar samples are red, while dissimilar samples are blue. Large groupings of red indicate tight clusters of samples^[Wei et al, https://github.com/taiyun/corrplot]

<img src="Figures/Ting/Sample_Correlation.png" width="400px" height="420px" />

In the case shown above, we might estimate there are likely 5 major subclusters in this data set^[Ting *et al*, Single-Cell RNA Sequencing Identifies Extracellular Matrix Gene Expression by Pancreatic Circulating Tumor Cells, Cell Reports, 2014]. This is just to provide an overview of how the samples correlate with each other using `Pearson correlation`. However, `Pearson correlation` can easily be altered by a few genes with high expression. Therefore, we tend to be more confident about the clustering results identified from using the NMF method which will be described in the next section.  
                    
                    
### Gene Network 
This tool identifies genes whose expression is correlated with one another across samples and creates a visual network to view these associations. Currently we use the absolute value of `Pearson correlation` to assess the strength of the connectivity^[Gandrud et al, https://christophergandrud.github.io/networkD3/]. Hub genes indicate genes that are highly connected to other genes within the same network -- these genes can often provide good candidate marker genes for further investigation.

In this example, we used the **Top MAD 1500** genes from the same dataset shown in the heatmap above for running this gene expression network module. 

<img src="Figures/Ting/SelTopMAD1500.png" width="400px" height="50px" />

<img src="Figures/Ting/Gene_network.png" width="400px" height="420px" />

## Sample scatter plot

You can compare the expression levels of whole transcriptome between two samples with a scatter plot. Genes that lie above the reference line are more highly expressed in the sample on the Y-axis, while genes that lie below the reference line are more highly expressed in the sample in the X-axis. Users can choose to log transform their gene counts, display *R^2^* coefficent of correlation, or show a two-fold change line to assess genes that have large changes in expression between the two samples. 

<img src="Figures/Ting/Scatter_plot.png" width="700px" height="450px" />

***

Continue on the next section [NMF](NMF.Rmd)

