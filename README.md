This R pipeline performs stability-informed cutoff analysis linking gene expression to patient survival, exemplified on TCGA datasets.  
It systematically evaluates candidate survival cutoffs across five key criteria:

1. **Cutoff region** – cutoff falls within the main expression density peak.  
2. **Group balance** – sample split is not overly skewed, ensuring statistical power.  
3. **Tumor–normal consistency** – direction of tumor vs. normal expression matches survival risk.  
4. **Bootstrap stability** – cutoff distribution shows low variability across resamples.  
5. **Spline support** – survival trend supported by spline-based modeling.  

<img width="3140" height="1250" alt="Stabicut_logo" src="https://github.com/user-attachments/assets/4f5c211c-31d9-4665-a5f6-aece0e6f13ea" />
