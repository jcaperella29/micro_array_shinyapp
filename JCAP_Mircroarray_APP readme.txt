JCAP_Mircroarray_APP readme. 

Summary of the Main Functions of this App 
This app can take your expression and phenotype data and perform several tasks when you click on the respective buttons. 
The first task that the app can perform is differential expression and feature selection. Differential expression is performed using Limma. The top 200 hits are then subjected to feature selection, which is performed using a random forest classifier. 20 genes that are deemed to be most important in terms of how great of an impact randomizing/permutating their expression values have on the performance of the classifier are output on the page. The resulting table contains the Ensembl IDs, HGNC symbols, p-values, respective log fold changes , and  Probe Ids of your top 20 genes. The table can be exported as a .csv file. This process also removes both genes that are lowly expressed among all samples  and outlier samples. 
The second task that this app can perform is to prepare a PCA plot of your data in which the points are colored with regard to the phenotype of interest (see the legend for details). Hovering your cursor over points on this plot will tell you which sample the points relate to.
The  third task that this app can perform is to prepare a UMAP  of your data colored in a similar fashion to the PCA plot. Hovering your cursor over points on this plot will tell you what sample the points relate to.
It  is noteworthy that the first, second and third tasks make use of parallel processing for parts of their respective functions making the application quite fast.
The fourth task that this application performs is using EnrichR to map  your top 20 genes to pathways. This can be performed for all of them, the upregulated genes, and the downregulated genes.
The fifth task that this app can perform is to prepare a volcano plot showing the relationship  between the adjusted p-value and log fold change with regard to the top 20 genes. Hovering your cursor over points on this plot will inform you regarding what gene (HGNC symbol) a given point relates to. 
All plots have sliders and points on the plot that can be used for saving the plot (camera icon) or expanding it (a box icon).
Prerequisites
1.You need to have your expression data in as a file on your device.
2. You need to have an accessible file containing your phenotype information. 
Please be aware that the program can only look at one phenotype at a time with this app. However, after examining your data with regard to one phenotype, you can re-examine your data with regard to another phenotype. 
The expression data and phenotype data can be in either .csv or .txt format. 


Mechanics/How to Use 
Inputs and Main Functions 
1.On the top of the purple sidebar on the left of the interface, you will see a file browse button with the label  “Input expression data” above it. Use this button to input your  expression . This can in be .txt or .csv format. When the file is uploaded, a blue ribbon will appear beneath the button that says ”Upload complete.”
2.Beneath that button, you will see another file browse button. This will have the label “Input Phenotype Data” above it. Use this button to input your phenotype data. This can also be either a .txt or .csv file. When the file is uploaded, a blue ribbon will appear beneath the button that says ”Upload complete.” 
3.Beneath that button, you will see a drop-down menu allowing you to tell the application what chip was used to collect your data. 
4.Below that you will see a button labeled “Perform Differential Expression Analysis.” Clicking this  button after uploading your expression data  and phenotype data will cause your data to undergo quality control (removing outlier samples, scaling and lowly expressed genes), differential expression, followed by feature selection  and annotation of the top twenty hits. The results of this process will be displayed in a table containing the Ensembl IDs, log fold changes,  HGNC symbols, full names of the protein product and the adjusted p-value. The table can be  viewed  by clicking on the “ Differential Expression Results” tab. FYI, the tab is considered “clicked “ whenever it is bright blue. A notification will also appear on the screen when the process starts and is completed.
5.Beneath that  button, there is a button labeled “Display PCA plot.” Clicking this button after uploading your expression and phenotype data will cause a PCA plot of your data to  appear in the “PCA plot” tab. Samples will be colored based on the phenotype. A legend is shown on the right of the plot. A notification will also appear on the screen when the process starts and is completed.
6. Beneath that  button, there is a button labeled “Display UMAP plot.” Clicking this button after uploading your expression  and phenotype data will cause a UMPA plot of your data  to appear in the “UMAP” tab. Both of the plots can be saved via right clicking the plot or using the camera icon on the plot.  A notification will also appear on the screen when the process starts and is completed.
7. Beneath that button,  there are 3 buttons that will allow you to use EnrichR to map your top genes to pathways if you click them after the table of the top 20 hits has been prepared. 
The first of those buttons runs EnrichR on all 20 genes. It is labeled as "Enrich Pathways(All Genes)." The button beneath  that button performs the analysis only on your upregulated genes. It is labeled as "Enrich Pathways(Upregulated Genes).” The last of those 3 buttons runs the analysis only on your downregulated genes. It is  labeled as "Enrich Pathways(Downregulated  Genes).”  A notification will also appear on the screen when the pathway analysis  process starts and is completed.
8. The next button going towards the bottom of the page is labeled "Display Volcano Plot." If this button is clicked after the table of top 20 genes has been produced, it creates a volcano plot displaying the relationship between the  adjusted p-value and log fold change of your top 20 genes. The plot will be shown in the “Volcano Plot” tab.  A notification will also appear on the screen when the process starts and  is completed. 

Sliders
Below the button in #8,  you will see sliders that will allow you to adjust parameters on the plots. 
1.The first slider  allows you to adjust how many PCA  components are shown in the PCA plot.
2.The second slider allows you to adjust the number of UMAP neighbors in the UMAP plot.
3.The third slider moves the volcano plot along the adjusted p-value axis.
4.The fourth slider moves the volcano plot along the log fold change axis.
Download Buttons
1.The next  set of buttons are download buttons. The first of the download buttons is labeled as  “Export Results as CSV” and will allow you to download the “Differential Expression Results” table  as a .csv file. Note: do not click this button before the “Differential Expression Results” has been prepared.
2.The second (going towards the bottom of the page) download button is labeled as "Export EnrichR results for all 20 genes as a CSV." Clicking this button after performing pathway enrichment on all 20 of your top genes will allow you to download the results of enrichment as a .csv file.
3.The third  (going towards the bottom of the page) download button is labeled as "Export EnrichR results for upregulated genes as a CSV." Clicking this button after performing pathway enrichment on your  upregulated  genes will allow you to download the results of enrichment as a .csv file.
4. The fourth and final download button is labeled as "Export EnrichR results for downregulated genes  as a CSV." Clicking this button after performing pathway enrichment on your  downregulated  genes will allow you to download the results of enrichment as a .csv file.
Notes
The “Perform Differential Expression Analysis,” Display PCA plot,” and “Display UMAP plot” will not perform their respective tasks unless the expression and phenotype files are uploaded. In the same way, the “Display Volcano plot” button will not run if the Differential Results table is not prepared. 
This readme file can be found by clicking the “Readme” tab.

