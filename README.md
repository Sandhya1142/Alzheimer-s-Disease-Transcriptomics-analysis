# Alzheimer's Disease Transcriptomics-analysis

**1. OVERVIEW**

RNAseq transcriptomics analyis has been implemented to identify the differentially expressed genes and perform Gene Set Enrichment Analysis in the chosen dataset consisting of 120 individuals consisting of healthy individuals and individuals across 4 different stages of AD- having mild cognitive impairment, cognitive impairment, Alzheimer's disease and both Alzheimer's disease & cognitive impairment.


**2. RESEARCH QUESTION ADDRESSED**

How does the gene expression change across different stages of cognitive impairement and Alzheimer's Disease  as compared to healthy individuals and what biological pathways are associated with these gene expression changes ?

**3. DATASET**

The following dataset has been taken from the **ENCODE** dataset (The Encyclopedia of DNA Elements) from 120 patients spanning healthy/control and disease-associated cognitive states.

**1** : No cognitive impairemnt (Healthy): 43 individuals


**2** : Mild Cognitive impairement: 31 individuals


**3** : Cognitive impairement: 4 individuals


**4** : Alzheimer's Disease: 35 individuals


**5** : Alzheimer's Disease & Cognitive impairement: 7 individuals




**4. ANALYSIS WORKFLOW**

<img width="1040" height="720" alt="Slide1" src="https://github.com/user-attachments/assets/560b6cee-70a8-437e-998e-ff25163f607f" />




**5. RESULTS**
Since significant number of DEGs were not observed for the groups under study, GSEA was performed to identify the set of genes that are activated or suppressed between the disease and control groups. The following were the groups under study:

<img width="1040" height="720" alt="Slide2" src="https://github.com/user-attachments/assets/a53ccbb8-e8f4-4515-a5d1-67f49dc776c6" />


Group 1: No cognitive impairment VS mild cognitive impairment
- **2 DEGs** were observed with padj < 0.05 & |Log2FC| > 1
- GSEA identified enrichment related to phosphatase activity and ubiquitin/ligase-associated processes.



Group 2: No cognitive impairment VS cognitive impairment
- **5 DEGs** were observed with padj < 0.05 & |Log2FC| > 1 
- GSEA indicated changes involving synaptic vesicle and exocytic vesicle processes and ATP synthesis   



Group 3: No cognitive impairment VS Alzheimer’s disease
- **0 DEGs** were observed  with padj < 0.05 & |Log2FC| > 1
- GSEA identified enrichment involving cytosolic ribosome and translation-related processes.


Group 4: No cognitive impairment VS Alzheimer’s disease & Cognitive impairment
- **1 DEG** were observed with padj < 0.05 & |Log2FC| > 1
- From GSEA, we observe that the top set of activated genes belong to oligodendrocyte differentialtion and regulation of RNA splicing.

**Gap Addressed using GSEA**: A lack of significant DEGs does not necessarily mean that there are no biologically meaningful transcriptomic changes.

 

**6. TOOLS AND SKILLSET ACQUIRED**
RNA-Seq analysis

Transcriptomics

R

Bash Scripting

Bioconductor

DESeq2

Differential Expression Analysis 

Gene Set Enrichment Analysis (GSEA)

Data Visualization

GitHub
