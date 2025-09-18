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
   - Merge ITS1, 5.8S, and ITS2 FASTA files into a single concatenated `ITS` file  

3. **Alignment & validation**  
   - Align curated FASTA with MAFFT  
   - Validate taxonomy assignments with SATIVA  

---

## 🚀 Step-by-Step Instructions

### Step 1 – Data Curation in R 🧹
Run the provided R script (`tsv2fasta.R`) to use the raw download file, remove missing/invalid entries, and generate curated FASTA + taxonomy files for ITS and trnL markers.

```bash
Rscript tsv2fasta.R
```

**Output files:**
- `ITS.fasta.gz`  
- `trnL.fasta.gz`  
- `ITS.tax`  
- `trnL.tax`  

---

### Step 2 – Running ITSx

Run **ITSx** on the input FASTA file to extract and save all ITS subregions using desired number of CPU cores:

```bash
ITSx -i {path/to/}ITS.fasta -o {output/folder/path} --cpu {number of threads} --save_regions all
ITSx -i {path/to/}trnL.fasta -o {output/folder/path} --cpu {number of threads} --save_regions all
```

### Step 3 – ITS Concatenation (Optional) 🔗
Run the provided R script (`construct_ITS_region.R`) to concatenate ITS1, 5.8S, and ITS2 into a single sequence per sample.

```bash
Rscript construct_ITS_region.R
```

**Output file:**
- `full_ITS.fasta` (concatenated ITS sequences)

---

### Step 4 – Alignment with MAFFT 🧩
Align the curated FASTA (example for `trnL`):

```bash
mafft --auto full_ITS.fasta > full.ITS_aligned.fasta
mafft --auto trnL_fasta > trnL_aligned.fasta
```

**Output file:**
- `trnL_aligned.fasta`

---

### Step 5 – Taxonomic Validation with SATIVA 🔍
Run [SATIVA](https://github.com/amkozlov/sativa) on the aligned FASTA + taxonomy table:

```bash
../sativa.py -s r-curation/full_ITS_aligned.fasta -t ITS.tax -T 2 -x bot -n output
../sativa.py -s r-curation/trnL_aligned.fasta              -t trnL.tax -T 2 -x bot -n output
```

**Output files:**
- SATIVA result files (mislabels, reports, etc.)

## Apply Sativa insights

Run the provided R script (`apply_sativa.R`) to correct misclassifications, using sativa .mis output and original .tax file to update it.

---

## 📦 Outputs
Final outputs include:
- Curated FASTA + taxonomy tables (`*.fasta.gz`, `*.tax`)
- Optional concatenated ITS sequences (`*.full.ITS.fasta`)
- SATIVA mislabel detection reports

---

## 💡 Notes & Tips
- File names are flexible; update script paths accordingly.  
- Ensure MAFFT and SATIVA are available in `$PATH`.  
- Large datasets may require additional memory/CPU.  
- SATIVA `-T` option controls number of threads.  
