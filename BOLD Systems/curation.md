# 🧬 Pipeline Guidelines

## 📖 Introduction
This pipeline processes retrieved sequence data, curates it into marker-specific FASTA files with taxonomy, optionally concatenates ITS subregions, and performs alignment and taxonomic validation.

---

## ⚙️ Requirements

### 🔧 Software
- **R** (≥ 4.0)
- **Python 3** (for `sativa.py`)
- **MAFFT** (sequence alignment tool)

### 📦 R Packages
- `data.table`
- `stringr`
- `seqinr`
- `R.utils`

### 📂 Input
- A tab-delimited file (named `download` in the example) with the following fields:
  - `nuc`: nucleotide sequence (can include `-` gaps)
  - `marker_code`: marker label (in our case: `ITS`, `ITS1`, `ITS2`, `trnL`)
  - `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`
  - `processid`: unique identifier

---

## 🗺️ Pipeline Overview

1. **Data curation in R**  
   - Load raw data and filter sequences  
   - Generate marker-specific FASTA files (`ITS`, `trnL`) and taxonomy tables  
   - Compress outputs  

2. **ITS concatenation (optional)**  
   - Merge ITS1, 5.8S, and ITS2 FASTA files into a single concatenated `full.ITS` file  

3. **Alignment & validation**  
   - Align curated FASTA with MAFFT  
   - Validate taxonomy assignments with SATIVA  

---

## 🚀 Step-by-Step Instructions

### Step 1 – Data Curation in R 🧹
Run the provided R script (`curation.R`) to filter raw data, remove missing/invalid entries, and generate curated FASTA + taxonomy files for ITS0 and trnL markers.

```bash
Rscript curation.R
```

**Output files:**
- `r-curation/ITS0.fasta.gz`  
- `r-curation/trnL.fasta.gz`  
- `r-curation/ITS0.tax`  
- `r-curation/trnL.tax`  

---

### Step 2 – ITS Concatenation (Optional) 🔗
Run the provided R script (`concatITS.R`) to concatenate ITS1, 5.8S, and ITS2 into a single sequence per sample.

```bash
Rscript concatITS.R
```

**Output file:**
- `*.full.ITS.fasta` (concatenated ITS sequences)

---

### Step 3 – Alignment with MAFFT 🧩
Align the curated FASTA (example for `trnL`):

```bash
mafft --auto r-curation/trnL.fasta > r-curation/trnL_aligned.fasta
```

**Output file:**
- `r-curation/trnL_aligned.fasta`

---

### Step 4 – Taxonomic Validation with SATIVA 🔍
Run [SATIVA](https://github.com/amkozlov/sativa) on the aligned FASTA + taxonomy table:

```bash
../sativa.py -s r-curation/trnL_aligned.fasta              -t r-curation/trnL.tax              -T 2 -x bot -n syntest
```

**Output files:**
- SATIVA result files (mislabels, reports, etc.)

---

## 📦 Outputs
Final outputs include:
- Curated FASTA + taxonomy tables (`r-curation/*.fasta.gz`, `r-curation/*.tax`)
- Optional concatenated ITS sequences (`*.full.ITS.fasta`)
- Aligned FASTA (`*_aligned.fasta`)
- SATIVA mislabel detection reports

---

## 💡 Notes & Tips
- File names are flexible; update script paths accordingly.  
- Ensure MAFFT and SATIVA are available in `$PATH`.  
- Large datasets may require additional memory/CPU.  
- SATIVA `-T` option controls number of threads.  
