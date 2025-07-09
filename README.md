# Compound Data Format Conversion Tool (SMILES, One-line Graph, gSpan Format)

This tool is an integrated CLI for converting compound data (SMILES format, one-line graph CSV, gSpan format) between different formats. The converted data can be directly used in [pmopt][1] and other tools.

[1]:https://github.com/takeuchi-lab/pmopt

---

## Required Python Packages

- `pysmiles`: SMILES notation molecular reading
- `pandas`: DataFrame operations
- `networkx`: Graph operations
- `rdkit`: Molecular structure conversion

Installation example:
```sh
pip install pysmiles pandas networkx rdkit
```

---

## Usage (converter_cli.py)

### Basic Commands

```sh
python converter_cli.py --mode smiles_to_graph --input input.csv --output output.csv --smiles-column SMILES --objective-column y --hydrogen 0
python converter_cli.py --mode graph_to_gspan --input input.csv --output output.graph
python converter_cli.py --mode gspan_to_graph --input input.graph --output output.csv
python converter_cli.py --mode graph_to_smiles --input input.csv --output output.csv
```

### Mode Descriptions

#### 1. SMILES → One-line Graph CSV
- `--mode smiles_to_graph`
- Required: `--input` input CSV/Excel, `--output` output CSV
- Optional: `--smiles-column` SMILES column name, `--objective-column` objective variable column name, `--hydrogen` include hydrogen (1/0)

#### 2. One-line Graph CSV → gSpan Format
- `--mode graph_to_gspan`
- Required: `--input` input CSV, `--output` output .graph
- Optional: `--graph-column` graph column name, `--spacing` blank lines between graphs, `--target-column` objective variable column name (output in header)

#### 3. gSpan Format → One-line Graph CSV
- `--mode gspan_to_graph`
- Required: `--input` input .graph, `--output` output CSV
- Optional: `--gspan-target-column` objective variable column name (default: objective)

#### 4. One-line Graph CSV → SMILES
- `--mode graph_to_smiles`
- Required: `--input` input CSV, `--output` output CSV
- Optional: `--graph-column` graph column name, `--other-columns` other column names to output together

---

## Input/Output Examples

### SMILES→One-line Graph CSV
Input CSV:
```csv
objective,SMILES
-5.89,O=C1c3ccccc3[Se]N1c2ccccc2
```
Command:
```sh
python converter_cli.py --mode smiles_to_graph --input input.csv --output output.csv --smiles-column SMILES --objective-column objective --hydrogen 1
```
Output CSV:
```csv
objective,graph
-5.89,v 0 O v 1 C ... e 0 1 = ...
```

### One-line Graph CSV→gSpan
Input CSV:
```csv
graph,objective
v 0 O v 1 C ... e 0 1 = ...,100
```
Command:
```sh
python converter_cli.py --mode graph_to_gspan --input input.csv --output output.graph --target-column objective
```
Output .graph:
```
t # 0 objective 100
v 0 8
v 1 6
...
e 0 1 4
...
```

### gSpan→One-line Graph CSV
Input .graph:
```
t # 0 objective 100
v 0 8
v 1 6
...
e 0 1 4
...
```
Command:
```sh
python converter_cli.py --mode gspan_to_graph --input input.graph --output output.csv
```
Output CSV:
```csv
objective,graph
100,v 0 O v 1 C ... e 0 1 = ...
```

### One-line Graph CSV→SMILES
Input CSV:
```csv
graph
v 0 O v 1 C ... e 0 1 = ...
```
Command:
```sh
python converter_cli.py --mode graph_to_smiles --input input.csv --output output.csv
```
Output CSV:
```csv
smiles
O=C1c3ccccc3[Se]N1c2ccccc2
```

---

## Option Details
- `--smiles-column` : Column name containing SMILES (default: smiles)
- `--objective-column` : Objective variable column name (default: objective)
- `--hydrogen` : Whether to include hydrogen (1: include, 0: exclude, default: 0)
- `--graph-column` : Column name containing one-line graph (default: graph)
- `--spacing` : Insert blank lines between graphs in gSpan output
- `--target-column` : Column name for objective variable (target value) to output in gSpan header (e.g., objective, y, etc., optional)
- `--gspan-target-column` : Objective variable column name for gSpan input (default: objective)
- `--other-columns` : Other column names to output together when generating SMILES

---

## Important Notes
- Please pay attention to the column names and format of input files
- Errors will occur if there are unknown element symbols or bond symbols
- If there are multiple input SMILES, the longest one will be selected
- Converted data can be directly used in [pmopt][1]

### 【Important】Information Loss and Limitations in Reverse Conversion
- **When converting SMILES→One-line Graph, the following information is lost:**
  - Stereochemical information ([C@H], E/Z, etc.)
  - Explicit aromaticity flags
  - Number of hydrogens (when hydrogen=0)
  - Atomic charge and radical information
- Therefore, **reverse conversion from One-line Graph→SMILES may not completely match the original SMILES.**
- Particularly for molecules containing aromatic rings or stereochemistry, **RDKit's Kekulize processing may cause errors (Can't kekulize mol.), making reverse conversion to SMILES impossible.**
- This is due to the specifications of this tool and RDKit/pysmiles, and depends on the structure and information content of the input molecule.
- If complete reverse conversion is required, **it is recommended to also save the original SMILES.**

---

## SMILES to SVG Image Generation (draw_smiles_to_svg.py)

A tool for visualizing molecular structures in SMILES notation as SVG images.

### Basic Commands

```sh
python draw_smiles_to_svg.py input.csv output_dir
python draw_smiles_to_svg.py input.csv output_dir --smiles_col SMILES
```

### Argument Descriptions

- `csv_path`: Input CSV file (file containing SMILES)
- `out_dir`: Output directory (where SVG files will be saved)
- `--smiles_col`: Column name containing SMILES (default: smiles)

### Usage Example

Input CSV (input.csv):
```csv
smiles,name
O=C1c3ccccc3[Se]N1c2ccccc2,compound1
CC(C)CC1=CC=C(C=C1)C(C)C(=O)O,compound2
```

Command:
```sh
python draw_smiles_to_svg.py input.csv svg_output
```

Output:
- `svg_output/0.svg`: SVG image of molecular structure from row 1
- `svg_output/1.svg`: SVG image of molecular structure from row 2

### Important Notes

- If there are invalid SMILES notations, those rows will be skipped and error messages will be displayed
- Output SVG images are 300x300 pixels in size
- If the output directory doesn't exist, it will be created automatically
