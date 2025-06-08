# NCBI-Searcher
Get all records from https://www.ncbi.nlm.nih.gov with filters, flags, and special functions to make scraping the NCBI NLM databases easier.

![image_2025-06-08_064340342](https://github.com/user-attachments/assets/f2fa5801-7a9d-4bd3-bf97-f29d91a7ed49)
![image_2025-06-08_114124324](https://github.com/user-attachments/assets/5321164e-4e85-4592-a2ad-d48560e7fcf2)

## Features
- Query NCBI all databases using PubMed-style syntax.
- Output in plain text or JSON.
- Supports batching, field-based output splitting, and verbose mode.
- Credential cycling to distribute requests across multiple API keys to optimize rate limits.
- Search Query can support Regular Expressions for more specific and advanced database searching.
- GUI Version to build queries if you are not familiar with the syntax.

Search queries are based from NCBI's e-utilities or PubMed syntax (AND, NOT, OR). E.g:  
`Human papillomavirus AND (E6[Title] OR E7[Title]) NOT partial[Title] -v -o HPVRecords.txt -s strain`  
`Gammapapillomavirus NOT (E1[Title] OR E7[Title]) -o GHPV_Records.txt`  

| Operator | Function                                   |
| -------- | ------------------------------------------ |
| `AND`    | Returns results containing **both** terms  |
| `OR`     | Returns results containing **either** term |
| `NOT`    | Excludes results containing the term       |
| `()`     | Groups terms for priority                  |

| Tag                | Searches...              |
| ------------------ | ------------------------ |
| `[Title]`          | Article title            |
| `[Abstract]`       | Abstract only            |
| `[Title/Abstract]` | Title and abstract       |
| `[Author]`         | Author name              |
| `[Journal]`        | Journal name             |
| `[MeSH Terms]`     | Medical Subject Headings |
| `[PDAT]`           | Publication date         |
| `[Affiliation]`    | Author's institution     |
| `[Language]`       | Language of the article  |

| Modifier | Function                                         |
| -------- | ------------------------------------------------ |
| `*`      |  Wildcard (e.g. vaccin* → vaccine, vaccination)  |
| `"`      | Exact phrase (e.g. "gene therapy")               |
| `:`      | Date range (e.g. 2020:2023[PDAT])                |

# Available settings
- `-v` Verbose mode.
- `-s` Split file between descriptor type. E.g., `accession` or `date`. Various flags exist that are not included in the record output but can be seen in the source, like `dateM`, `dateD`, `dateY`
- `-o` Ouput file, can be any file type.
- `-j` JSON Ouput mode.
- `-c` Credentials to use, they should be written in `apikey:email` format. If you have multiple credentials they can be used with -c `apikey:email`,`apikey2:email2`,`apikey3:email3`
- `-b` Batch size, default is `50` but recommended by NCBI to be no further than `300`. This mostly controls how much records it will gather before switching to the next credential. And to make data gathering faster and organized.



Example record output:

```
(accession) Accession: XQZ12201
(date) Date: 13-APR-2025
(strain) Strain: 16
(organism) Organism: Human papillomavirus 16
(source) Isolation Source: cervical intraepithelial neoplasia CIN 3
(location) Geo Location: Russia: Saint Petersburg
(gene) Gene(s): 
(description) Description: protein E6 [Human papillomavirus 16]
(sequence) Sequence:
MHQKRTAMFQDPQERPRKLPQLCTELQTTIHDIILECVYCKQQLLRREVY
(url) GenBank URL: https://www.ncbi.nlm.nih.gov/protein/XQZ12201.1
```

The descriptors beside the actual names are used for the -s flag, to split file out based on them.
