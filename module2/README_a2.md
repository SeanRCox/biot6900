Sean Cox

I chose to cover colon adenocarcinoma. Colon cancer is the third most common cancer diagnosed in both men and women in the United States and the second most common cause of cancer deaths.

I was able to get everything working.

Colon Cancer 

All datasets are open
CPTAC Colon Cancer  https://www.linkedomics.org/data_download/CPTAC-COAD/
USCS TCGA Colon Cancer https://xenabrowser.net/datapages/?cohort=TCGA%20Colon%20Cancer

For the analysis, I integrated colon cancer RNA, protein, and mutation evidence from CPTAC with TCGA tumor-versus-normal RNA data, harmonizing the layers to shared gene identifiers. I compared RNA and protein patterns, combined the gene-level layer scores, and exported the top-ranked candidates to `targets_ad.csv` for the next module. To run it, open `BIOT6900_Module2_Starter.ipynb` from the `module2` folder and run the notebook cells in order, with the required source datasets available in `module2/data`. The assignment's purpose was to practice multi-omics target prioritization when the evidence layers do not all come from matched samples.
