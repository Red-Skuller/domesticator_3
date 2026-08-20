# Domestica

Domestica is a Python command-line tool. It automates DNA sequence optimization. The software reads protein sequences, reverse-translates them into DNA, and inserts them into a GenBank vector template. It applies constraints to optimize the sequence and can send the sequence to vendor APIs to check synthesis complexity.

## Features

*   **Reverse Translation**: Converts protein sequences to DNA using specified codon tables.
*   **Template Integration**: Replaces a target region in a GenBank template with the new DNA sequence.
*   **Sequence Optimization**: Uses `dnachisel` to resolve constraints (e.g., GC content, hairpins) and minimize k-mers.
*   **Sequence Padding**: Adds padding bases (A) to sequences that are shorter than a specified minimum length.
*   **Vendor Validation**: Connects to IDT or ThermoFisher APIs to check if the sequence is acceptable for synthesis.
*   **Parallel Processing**: Processes multiple sequences concurrently.

## Installation

Domestica requires Python. Install the package and its dependencies:

```bash
pip install -r requirements.txt
```

*(Note: Ensure you have `typer`, `dnachisel`, `biopython`, `pandas`, `openpyxl`and `httpx` installed).*

## Configuration

Configure the tool with environment variables. You can write these variables in a `.env` file in your working directory.

| Variable Name | Description | Required |
| :--- | :--- | :--- |
| `DOMESTICA_VENDOR_TARGET` | Vendor to use for evaluation (`idt` or `thermofisher`). | No |
| `DOMESTICA_PRODUCT` | Synthesis product type (e.g., `eblocks`, `genes`). | No |
| `DOMESTICA_MAX_WORKERS` | Maximum number of parallel processes (1-64). | No |
| `DOMESTICA_MIN_LENGTH` | Minimum base pair length for the final sequence. | No |
| `DOMESTICA_IDT_CLIENT_ID` | IDT API Client ID. | Yes (if vendor=idt) |
| `DOMESTICA_IDT_CLIENT_SECRET` | IDT API Client Secret. | Yes (if vendor=idt) |
| `DOMESTICA_IDT_USERNAME` | IDT Account Username. | Yes (if vendor=idt) |
| `DOMESTICA_IDT_PASSWORD` | IDT Account Password. | Yes (if vendor=idt) |
| `DOMESTICA_THERMOFISHER_CLIENT_ID` | ThermoFisher API Client ID. | Yes (if vendor=thermofisher) |
| `DOMESTICA_THERMOFISHER_CLIENT_SECRET`| ThermoFisher API Client Secret. | Yes (if vendor=thermofisher) |

## Usage

Use the `optimize` command to run the pipeline.

```bash
python -m domestica optimize [INPUT_PATH] [OUTPUT_PATH] --template [TEMPLATE_PATH] [OPTIONS]
```

### Arguments and Options

*   `INPUT_PATH`: Path to the input file. If you do not provide this, the tool optimizes the empty template.
*   `OUTPUT_PATH`: Path to save the output file(s).
*   `-t, --template`: **(Required)** Path to the `.gb` or `.genbank` template file.
*   `-v, --vendor`: Set the vendor for complexity evaluation (`idt` or `thermofisher`).
*   `-p, --product`: Set the product type (Default: `eblocks`).
*   `-w, --workers`: Set the number of parallel workers.
*   `-m, --min-length`: Set the minimum sequence length in base pairs (Default: `300`).
*   `--verbose`: Enable debug logging.

### Examples

Optimize sequences from a FASTA file and write GenBank outputs:

```bash
python -m domestica optimize input.fasta output.gb -t vector_template.gb -v idt -p eblocks --verbose
```

Optimize sequences from an Excel file and output a CSV report:

```bash
python -m domestica optimize sequences.xlsx report.csv -t vector_template.gb -w 4
```

## File Formats

### Inputs
Domestica accepts these input file types:
*   **FASTA** (`.fasta`, `.fa`): Extracts record IDs and protein sequences.
*   **Excel** (`.xlsx`, `.xls`): Reads the first column as the record ID and the second column as the protein sequence.

### Outputs
Domestica writes to these output file types:
*   **GenBank** (`.gb`, `.genbank`): Writes one GenBank file for each successful sequence. The file name includes the record ID.
*   **Excel/CSV** (`.xlsx`, `.xls`, `.csv`): Writes a flat table containing optimization status, scores, and sequences.

## Template Preparation

The input GenBank template must obey two structural rules:

1.  **Target Region**: The template must contain at least one feature with the text `!INSERT` or `INSERT` in its `label`, `note`, or `locus_tag`. Domestica replaces this region with the reverse-translated DNA.
2.  **Codon Table (Optional)**: To enforce a specific genetic translation table, add an `@EnforceTranslation` tag to a feature that overlaps the `INSERT` region.
    *   Format: `@EnforceTranslation(genetic_table="TableName")`
    *   If you do not specify a table, Domestica uses `Standard`.