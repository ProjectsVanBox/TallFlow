# TallFlow
Code for all analyses performed in Poort et al. 2023?, **"Transient differentiation-state plasticity during acute lymphoblastic leukemia initiation"**

Author: Vera Poort <br>
Contact1: v.m.poort@prinsesmaximacentrum.nl <br>
Contact2: r.hagelaar@prinsesmaximacentrum.nl

Flow cytometry data is deposited here: [FlowRepository.org](http://flowrepository.org/) 
+ Thymocyte subtyping: Repository ID: FR-FCM-Z5TV.
+ T-cell leukemia immunophenotyping: Repository ID: FR-FCM453Z5UY.

WGS data is deposited here: https://ega-archive.org/. Archive ID: [EGAD00001011126](https://ega-archive.org/search-results.php?query=EGAD00001011126)


## Table of content
### Flow cytometry analysis
*Figures 1 & 2* <br>
[QC raw .fcs files](https://github.com/ProjectsVanBox/TallFlow/blob/main/PythonGating/QC_flowfiles.R) <br>
[Semi-automatic gating](https://github.com/ProjectsVanBox/TallFlow/tree/main/PythonGating) <br>
[Visualisations](https://github.com/ProjectsVanBox/TallFlow/tree/main/Barcharts) <br>
[Comparing Flow panel with EuroFlow and CITE data](https://github.com/ProjectsVanBox/TallFlow/blob/main/CompFlow.R)

### Variant filtering
All single cell WGS was first mapped and QC'ed using [NF-IAP](https://github.com/ToolsVanBox/NF-IAP) <br>
For single cell WGS samples we made use of PTa Analysis TOolkit [PTATO](https://github.com/ToolsVanBox/PTATO) <br>
+ [pt2229](https://github.com/ProjectsVanBox/TallFlow/tree/main/PTA_Dev_2229) <br>
+ [pt2322](https://github.com/ProjectsVanBox/TallFlow/tree/main/PTA_Dev_2322) <br>
+ [pt2283D](https://github.com/ProjectsVanBox/TallFlow/tree/main/PTA_pt2283D) <br>

For bulk WGS we made use of [SMuRF](https://github.com/ToolsVanBox/SMuRF) <br>
[bulkSeq SMuRF](https://github.com/ProjectsVanBox/TallFlow/tree/main/BulkSeq_SMuRF) <br>
[CallableLoci](https://github.com/ProjectsVanBox/TallFlow/tree/main/CallableLociBulk) <br>


[Drivers](https://github.com/ProjectsVanBox/TallFlow/tree/main/Drivers) <br>
[True Positive Rate](https://github.com/ProjectsVanBox/TallFlow/tree/main/TruePosRate) <br>
[Structural variants](https://github.com/ProjectsVanBox/TallFlow/tree/main/SVs) <br>

### Ageline
*Figure 3* <br>
[Ageline](https://github.com/ProjectsVanBox/TallFlow/tree/main/AgeLine)

### Mutational Patterns
*Figure 3* <br>
  [SBS signatures](https://github.com/ProjectsVanBox/TallFlow/tree/main/Mutpatterns) & [plotting](https://github.com/ProjectsVanBox/TallFlow/blob/main/MutationalPatterns_Addon.R) <br>
[Indel signatures](https://github.com/ProjectsVanBox/TallFlow/tree/main/IndelSignatures)

### Treebuilding
*Figure 4* <br>
[Trees](https://github.com/ProjectsVanBox/TallFlow/tree/main/TreeBuilding)

### V(D)J analysis
*Figure 5* <br>
  [MIXCR](https://github.com/ProjectsVanBox/TallFlow/tree/main/MIXCR) <br>
  [RSS motive Search](https://github.com/ProjectsVanBox/TallFlow/tree/main/RSSmotifSearch)








