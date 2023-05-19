# Guide: Inferring individual-specific transcription factor activity   

![Framework](https://github.com/xl27/GTEx_aQTLs/blob/main/images/FigS1A.png)  

## Example Datasets:
We use GTEx skeletal muscle RNA-seq data and NR4A1 perturbation data  as examples:   
GTEx v8 data: [open access datasets](https://gtexportal.org/home/datasets)   
ENCODE data: [CRIPSRi experiments](https://www.encodeproject.org/experiments/ENCSR357LVC/)     
           
Results of TF activity are available here: [data](http://bussemakerlab.org/papers/Li-aQTL/TF_activity/)   

## Steps:
1) Preparation of RNA-Seq data from the GTEx project: `Process_GTEx_Count_Matrix.R`	
2) Selection of gene pairs: `Select_Gene_Pairs.Rmd`
3) Generating TF perturbation response signature using ENCODE CRISPRi data: `Perturbation_Signature.Rmd`	
4) Inferring TF activity based on gene-pair model:  
`Pair_Level_Inference_Model.py NR4A1 Muscle_Skeletal`	


## Dependencies:
- R v4.0.1   
- featureCounts v2.0.0    
- Python v3   
- TensorFlow v2.9.1    
- pandas/numpy/sys/sklearn   
