# Phylogenetic-placement-for-animal-aeDNA-in-Lake-Tjornin-of-Reykjavik
Computational workflow and scripts for phylogenetic placement of animal ancient environmental DNA (aeDNA) from Lake Tjörnin, Reykjavík. Includes reference panel construction, tree inference, target species enrichment, and PathPhynder-based lineage assignment.

## ⚙️ Workflow Overview

### **1️⃣ Reference Panel Construction**
`1_Reference_panel_construction.sh`  
Builds the mitochondrial reference panel for downstream phylogenetic inference.

### **2️⃣ Phylogenetic Tree Construction**
`2_Phylogenetic_tree_construction.sh`  
Generates phylogenetic tree for PathPhynder.

### **3️⃣ Target Species Enrichment**
`3_Target_species_enrichment.sh`  
Identifies sequences belonging to target taxa for focused phylogenetic placement.

### **4️⃣ PathPhynder Analysis**
`4_PathPhynder.sh`  
Runs PathPhynder to assign query sequences to reference nodes on the phylogenetic tree for lineage placement and taxonomic validation.

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
- Mitochondrial genome sequences used as reference data.  
- Ultrametric tree files infered with BEAST.  

These files are retained for **verification**, **reproducibility**, and **inspection** of each stage in the workflow.

---
