# BIOT 6900 · Module 7

Sean Cox

In this module, I built a colon cancer knowledge graph in Neo4j using genes, phenotypes, pathways, drugs, and variants. I queried the graph with Cypher and used NetworkX PageRank to prioritize genes. The graph organizes evidence but does not establish causation or clinical effectiveness.

To run the assignment, generate `data/targets_ad_updated_data.csv` with `python module7/update_targets_ad_data.py`, then open `BIOT6900_Module7_Starter.ipynb` and run the cells in order. Set Neo4j credentials in the git-ignored `.env` file first. The script combines Module 2/3 data with GWAS Catalog results; blank GWAS values mean no matching association was found. Check gene symbols and pathway memberships before interpreting pathway connections.
