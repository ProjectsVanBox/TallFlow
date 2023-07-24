# TallFlow
Code for all analyses performed in Poort et al. 2023?, "Transient differentiation-state plasticity during acute lymphoblastic leukemia initiation" 

Author: Vera Poort
Contact1: v.m.poort@prinsesmaximacentrum.nl
Contact2: r.hagelaar@prinsesmaximacentrum.nl

Flow cytometry data is deposited here: [FlowRepository.org](http://flowrepository.org/) 
Thymocyte subtyping: Repository ID: FR-FCM-Z5TV. T-cell leukemia immunophenotyping: Repository ID: FR-FCM453
Z5UY.
WGS data is deposited here: https://ega-archive.org/. Archive ID: [EGAD00001011126](https://ega-archive.org/search-results.php?query=EGAD00001011126)


##Table of content
###Flow cytometry analysis
Figures 1 & 2
[QC raw .fcs files](https://github.com/ProjectsVanBox/TallFlow/blob/main/PythonGating/QC_flowfiles.R)
[Semi-automatic gating](https://github.com/ProjectsVanBox/TallFlow/tree/main/PythonGating)
[Visualisations](https://github.com/ProjectsVanBox/TallFlow/tree/main/Barcharts)
[Comparing Flow panel with EuroFlow and CITE data]()

###Variant filtering
All single cell WGS was first mapped and QC'ed using [NF-IAP](https://github.com/ToolsVanBox/NF-IAP)
For single cell WGS samples we made use of PTa Analysis TOolkit [PTATO](https://github.com/ToolsVanBox/PTATO)
[pt2229](https://github.com/ProjectsVanBox/TallFlow/tree/main/PTA_Dev_2229)
[pt2322](https://github.com/ProjectsVanBox/TallFlow/tree/main/PTA_Dev_2322)
[pt2283D](https://github.com/ProjectsVanBox/TallFlow/tree/main/PTA_pt2283D)
For bulk WGS we made use of [SMuRF](https://github.com/ToolsVanBox/SMuRF)
[bulkSeq SMuRF](https://github.com/ProjectsVanBox/TallFlow/tree/main/BulkSeq_SMuRF)
[CallableLoci](https://github.com/ProjectsVanBox/TallFlow/tree/main/CallableLociBulk)
[Drivers](https://github.com/ProjectsVanBox/TallFlow/tree/main/Drivers)
[True Positive Rate](https://github.com/ProjectsVanBox/TallFlow/tree/main/TruePosRate)
[Structural variants](https://github.com/ProjectsVanBox/TallFlow/tree/main/SVs)

###[Ageline](https://github.com/ProjectsVanBox/TallFlow/tree/main/AgeLine)
Figure 3

###Mutational Patterns
Figure 3
[SBS signatures](https://github.com/ProjectsVanBox/TallFlow/tree/main/Mutpatterns) & [plotting](https://github.com/ProjectsVanBox/TallFlow/blob/main/MutationalPatterns_Addon.R)
[Indel signatures](https://github.com/ProjectsVanBox/TallFlow/tree/main/IndelSignatures)

###Treebuilding
Figure 4
[Trees](https://github.com/ProjectsVanBox/TallFlow/tree/main/TreeBuilding)

###V(D)J analysis
Figure 5
[MIXCR](https://github.com/ProjectsVanBox/TallFlow/tree/main/MIXCR)
[RSS motive Search](https://github.com/ProjectsVanBox/TallFlow/tree/main/RSSmotifSearch)


###[Ageline](https://github.com/ProjectsVanBox/TallFlow/tree/main/AgeLine)
Figure 3







