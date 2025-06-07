# NCBI-Searcher
Get all records from https://www.ncbi.nlm.nih.gov with filters, flags, and special functions to make scraping the NCBI NLM databases easier.

[SQ] Search Query Flags:
-v                 Verbose mode
-o [filename]      Set output file (default: records.txt)
-j                 Output as JSON instead of plain text
-s [field]         Split output files by descriptor field (e.g., strain)
-c [cred_string]   Credentials (email:api_key), comma-separated for multiple credentials for faster scraping (due to throttling by the NCBI)
-b [num]           Sets record fetching batch size, default: 50

[DB] Database Query Commands:
help               Show this help message
exit               Exit the program

Search queries are based from NCBI's e-utilities or PubMed syntax (AND, NOT, OR). E.g:
`Human papillomavirus AND (E6[Title] OR E7[Title]) NOT partial[Title]`
`Gammapapillomavirus NOT (E1[Title] OR E7[Title])`

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

Example output:

```
(accession) Accession: VE6_HPV94
(strain) Strain: 94
(organism) Organism: Human papillomavirus 94
(source) Isolation Source: Unknown
(location) Geo Location: Unknown
(gene) Gene(s): E6
(description) Description: RecName: Full=Protein E6
(sequence) Sequence:
MSMGAQEPRNIFLLCRNCGISFEDLRLCCVFCTKQLTVAELTAFALRELNLVWKAGVPYGACARCLLLQGIARRLKYWQYSYYVEGVEEETKESINTQQIRCYTCHKPLVKEEKDRHRNERRRLHKISGYWRGCCAYCWTRCTVRIPQ
(url) GenBank URL: https://www.ncbi.nlm.nih.gov/protein/Q705D2.1
```

The descriptors beside the actual names are used for the -s flag, to split file out based on them.
