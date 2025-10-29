# Animal mitochondrial eDNA phylogenetic placement for the Tjörnin dataset
This repo contains the workflow and scripts for mitochondrial phylogenetic placement of animal ancient eDNA from Lake Tjörnin, Reykjavík, including reference panel assembling, tree construction, target DNA extraction, and PathPhynder-based lineage assignment.
Codes/commends for running each step of the analysis are contained in the given scrips. One would need to edit the parameters on top of each script for modifying specific file names, CPU threads, and memory limits before running the codes/ commends.

## ⚙️ Workflow Overview

### **1️⃣ Reference Panel Assembling**
`1_Reference_panel_construction.sh`  
Builds the mitochondrial reference panel for phylogenetic inference.

### **2️⃣ Phylogenetic Tree Construction**
`2_Phylogenetic_tree_construction.sh`  
Generates phylogenetic tree.

### **3️⃣ Target DNA Extraction**
`3_Target_reads_extraction.sh`  
Identifies reads from target taxa.

### **4️⃣ PathPhynder Placement**
`4_PathPhynder.sh`  
Performs PathPhynder to assign query reads to reference nodes on the phylogenetic tree for lineage placement.

---

## 🧰 Supporting Scripts

| File | Description |
|:--|:--|
| **make_beast_xml.py** | Generates BEAST input XML files from the provided template with parameters. |
| **nex_to_nwk.py** | Converts tree files from NEXUS (`.nex`) to PathPhynder-required Newick (`.nwk`) format. |
| **run_ngsLCA.py** | Performs taxonomic assignment using the ngsLCA pipeline. |
| **ultrametric_tree_template.xml** | Template for BEAST XML input specifying priors, models, and clock parameters. |

---

## 📁 Intermediate Files

The folder **`Intermediate_files/`** contains:
- Mitochondrial genome sequences used in reference panel.  
- Ultrametric tree files infered with BEAST.  

---
