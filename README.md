# Alzheimer's Disease Transcriptomics-analysis

**1. OVERVIEW**

RNAseq transcriptomics analyis has been implemented to identify the differentially expressed genes and perform Gene Set Enrichment Analysis in the chosen dataset consisting of 120 individuals consisting of healthy individuals and individuals across 4 different stages of AD- having mild cognitive impairment, cognitive impairment, Alzheimer's disease and both Alzheimer's disease & cognitive impairment.


**2. RESEARCH QUESTION ADDRESSED**

How does the gene expression change across different stages of cognitive impairement and Alzheimer's Disease  as compared to healthy individuals and what biological pathways are associated with these gene expression changes ?

**3. DATASET**

The following dataset has been taken from the **ENCODE** dataset (The Encyclopedia of DNA Elements) from 120 patients spanning healthy/control and disease-associated cognitive states.

**1** : No cognitive impairemnt (Healthy): 


**2** : Mild Cognitive impairement:


**3** : Cognitive impairement:


**4** : Alzheimer's Disease:


**5** : Alzheimer's Disease & Cognitive impairement: 




**4. ANALYSIS WORKFLOW**


<img width="1040" height="720" alt="Analysis_Workflow" src="https://github.com/user-attachments/assets/79c1fa9c-051c-456b-86f4-52e771e2dd06" />



**5. RESULTS**

The followig were the groups under study:

Group 1: No cognitive impairment VS mild cognitive impairment
- 2 DEGs were observed with padj < 0.05 & Log2FC > 1
- From GSEA, we observe that the top set of activated genes belong to phosphatase activity, ubiquitin and ligase complex.   

Group 2: No cognitive impairment VS cognitive impairment
- 5 DEGs were observed with padj < 0.05 & Log2FC > 1 and Log2FC < 1
- From GSEA, we observe that the top set of suppressed genes belong to synaptic vesicle membrane, exocytic vesicle membrane and ATP synthesis pathway.   

Group 3: No cognitive impairment VS Alzheimer’s disease
- 0 DEGs were observed  with padj < 0.05 & Log2FC > 1 / Log2FC < 1
- From GSEA, we observe that the top set of suppressed genes belong to cytosolic ribosome component, cytosolic large ribosomal subunit and cytosolic translation.    

Group 4: No cognitive impairment VS Alzheimer’s disease & Cognitive impairment
- 1 DEGs were observed with padj < 0.05 & Log2FC < 1
- From GSEA, we observe that the top set of activated genes belong to oligodendrocyte differentialtion and regulation of RNA splicing.

All the downstream analysis results have been added to the /Figures folder. 

