# OptMAP
 Evolution and design of metabolic networks using quality-diversity optimization algorithms

Metabolic engineering in the wet lab is limited by the number of possible combinations of knock-ins, knockouts, media compositions etc. that can be tested. 
The search or design space created by all the possible combinations of manipulations that are possible in metabolic engineering is so large that navigating it in the wet lab is a long and resource consuming journey. Therefore, we introduce a novel computational tool, called OptMAP, which is a novelty search based framework that combines quality diversity evolutionary algorithms with genome scale modeling for intelligent in silico metabolic engineering.
The main algorithm is inspired by the quality-diversity algorithm called MAP-Elites and aims to predict strain design strategies for maximal production of specific compounds.
We use an evolutionary algorithm in order to model the evolution of metabolic networks in silico in the hope to discover novel high performing strains.
Metabolic engineering manipulations suggested by OptMAP include gene knockout strategies for microbial strain optimization and medium composition suggestions, with the possibility to expand to up- and down-regulation of genes, the addition of heterologous pathways and much more.
The OptMAP framework was validated by identifying metabolic manipulations for the production of four industrially relevant biochemicals using the genome-scale metabolic models for Escherichia coli str. K-12 substr. MG1655}
OptMAP was applied to the overproduction of succinate, acetate and ethanol and the production of flavanones.

![alt text](https://github.com/shvhoye/OptMAP/blob/main/optmap_workflow.jpg?raw=true)