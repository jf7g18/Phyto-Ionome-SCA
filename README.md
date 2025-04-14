# Phytoplankton Single-Cell Ionome Analysis Using Microscopy and LA-TOF-ICP-MS

**Title**: Data and Code for *"A Workflow for Deriving Single-Cell Phytoplankton Ionomes via Coupled Microscopy and LA-TOF-ICP-MS"*

**Authors**: Joe S. K. Furby  
**Corresponding Author Contact**: [jf7g18@soton.ac.uk]  
**Associated Manuscript**: Submitted to *Journal of Analytical Atomic Spectroscopy* (JAAS – Fast Transient Signals)

---

## Overview

This repository contains the dataset and R script for deriving single-cell ionomes of natural phytoplankton assemblages. The workflow integrates **correlative multimodal imaging (CMI)** combining:

- Optical light microscopy for cell imaging  
- U-Net model semantic segmentation from optical imagery for single cell isolation 
- LA-TOF-ICP-MS for multi-elemental quantification

Single cell ionomes are extracted from the stacked LA-TOF-ICP-MS element concentration maps using segmentation masks (manual and U-Net model predicted). **NIST SRM 610 and 612** glasses serve as external standards, and **phosphorus (P)** is used as an internal standard to normalize cellular ionomes for semi-quantitative multi-element analysis at the single-cell level.

The goal is to enable **reproducible and high-throughput ionome analysis** of natural phytoplankton communities using an open, script-driven pipeline.

---

## Repository Contents

```plaintext
📁 data/
 ├── LA TOF-ICP-MS/                                
 │   ├── Calibration/
 │   │   ├── NIST_610_scan/                        # 3 ablation line .csv files (435 shots each) of NIST610 SRM glass    
 │   │   ├── NIST_612_scan/                        # 3 ablation line .csv files (459 shots each) of NIST612 SRM glass 
 │   │   ├── NIST Calibration Results.csv          # Model-I regression slopes and R² values from NIST calibration
 │   │   ├── NIST Concentrations.csv               # Element concentrations (µg/g) used for calibration
 │   │   └── NIST values Jochum et al. 2011.csv    # Published SRM NIST glass element concentrations (Jochum et al., 2011)
 │
 │   ├── STATION1/                                 # 90 ablation line .csv files (512 shots each)
 │   ├── STATION2/                                 # 100 ablation line .csv files (512 shots each)
 │   ├── STATION3/
 │   ├── STATION4/
 │   ├── STATION5/
 │   └── STATION6/                                 # 100 ablation line .csv files (503 shots each)

 ├── Optical Images/                               # Light microscope images of phytoplankton on sample filters (.jpg)
 ├── Output/                                       # Output figures (.png) and final single cell ionome dataset (Total Cell Ionomes.csv)
 ├── Segmentation Masks/                           # Manual and U-Net model-generated segmentation masks (.tiff)
 ├── Species and Volume/                           # Annotated cell species and volume estimates from manual and U-Net masks (.csv)
 └── U-Net Model Validation/                       # Segmentation masks of validation images (.tiff)

📁 scripts/
 └── Phyto Ionome SCA.R                            # Main R script for processing and analysis pipeline

📄 README.md                                        # This file
📄 metadata.yaml                                    # Optional: Structured metadata for repository contents
```