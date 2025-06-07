"""
MIT License

Copyright (c) [2025] [Xie Yuzho]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from Bio import Entrez, SeqIO
import re
from datetime import datetime
from urllib.error import HTTPError
import json
from contextlib import contextmanager
from itertools import cycle

@contextmanager
def dummy_context_manager():
    """
    A no-op context manager for conditional file handling. 
    In simple words: Dummy placeholder.
    """
    yield None

def parse_credentials(cred_string):
    """
    Parse credential string in format 'email:api_key,email2:api_key2,...' into list of tuples.
    """
    pairs = [p.strip() for p in cred_string.split(",")]
    return [tuple(p.split(":")) for p in pairs if ":" in p] # Looks complicated..
    # Chronological order of credentials in this function:
    # email:api_key,email2:api_key2
    # [(email:api_key) , (email2:api_key2)]
    # [(email , api_key) , (email2 , api_key2)]

def clean_query_flags(query, flag):
    """
    Remove a flag (e.g., '-v') and its argument (if any) from the query string.
    Returns cleaned query and/or the flag's value if applicable.
    """
    pattern = None
    flag_value = None
    if flag in ['-o', '-s', '-c', '-b']:
        # Flags that take an argument
        pattern = re.compile(rf"{flag}\s+([^\s]+)")
        match = pattern.search(query)
        if match:
            flag_value = match.group(1)
            query = pattern.sub("", query).strip()
    else:
        # Flags without argument like -v or -j
        if flag in query:
            flag_value = True
            query = query.replace(flag, "").strip()
    return query, flag_value

def main():
    credentials = []
    credential_cycle = cycle([])

    while True:
        try:
            total_count = 0
            output_file = "records.txt"

            db = input("[DB] Enter command (e.g., help, exit) or NCBI database (e.g., protein, nuccore): ").strip().lower()
            if db in {"exit", "x"}:
                print("Exiting...")
                break

            if db == "help":
                print("""
                [SQ] Search Query Flags:
                -v                 Verbose mode
                -o [filename]      Set output file (default: records.txt)
                -j                 Output as JSON instead of plain text
                -s [field]         Split output files by descriptor field (e.g., strain)
                -c [cred_string]   Credentials (email:api_key), comma-separated for multiple

                [DB] Database Query Commands [You are here]:
                help               Show this help message
                exit               Exit the program
                """)
                continue

            if not db:
                print("Database input cannot be empty.")
                continue

            query = input("[SQ] Enter your NCBI search query: ").strip()
            if not query:
                print("Search query cannot be empty.")
                continue

            # Parse flags and options
            query, output_file_opt = clean_query_flags(query, "-o")
            if output_file_opt:
                output_file = output_file_opt

            query, split_field = clean_query_flags(query, "-s")
            query, verbose_mode = clean_query_flags(query, "-v")
            verbose_mode = bool(verbose_mode)
            query, json_mode = clean_query_flags(query, "-j")
            json_mode = bool(json_mode)
            query, batch_size_str = clean_query_flags(query, "-b")
            batch_size = int(batch_size_str) if batch_size_str else 50

            query, cred_string = clean_query_flags(query, "-c")
            if cred_string:
                credentials = parse_credentials(cred_string)
                credential_cycle = cycle(credentials)
            else:
                credentials = []
                credential_cycle = cycle([])

            if verbose_mode:
                print(f"[VERBOSE] Searching '{db}' for: '{query}'")
                if credentials:
                    print(f"[VERBOSE] Using {len(credentials)} credential(s) for API access.")
                print(f"[VERBOSE] Output file: {output_file}")
                if split_field:
                    print(f"[VERBOSE] Splitting output files by field: {split_field}")
                print(f"[VERBOSE] Output format: {'JSON' if json_mode else 'Plain text'}")

            # Set Entrez credentials from credential list (given from -c flag) if available
            try:
                api_key, email = next(credential_cycle)
                Entrez.email = email
                Entrez.api_key = api_key
            except StopIteration:
                Entrez.email = None
                Entrez.api_key = None
                if verbose_mode:
                    print("[WARN] No credentials provided. Get API key at https://www.ncbi.nlm.nih.gov/account")

            # Search for IDs
            handle = Entrez.esearch(db=db, term=query, retmax=50000)
            record = Entrez.read(handle)
            ids = record.get("IdList", [])
            print(f"Found {len(ids)} record(s).")

            if not ids:
                print("No results found.")
                continue

            # Open main output file if not splitting
            split_files = {}
            with open(output_file, "w") if not split_field else dummy_context_manager() as main_fh:
                for batch_start in range(0, len(ids), batch_size):
                    try:
                        api_key, email = next(credential_cycle)
                        Entrez.email = email
                        Entrez.api_key = api_key
                    except StopIteration:
                        Entrez.email = None
                        Entrez.api_key = None

                    batch_ids = ids[batch_start:batch_start + batch_size] # Batches ids in sizes of batch_size

                    if verbose_mode:
                        print(f"[VERBOSE] Fetching batch {batch_start // batch_size + 1}: {len(batch_ids)} records using credential: {email if Entrez.email else 'None'}")

                    fetch_handle = Entrez.efetch(db=db, id=",".join(batch_ids), rettype="gb", retmode="text")

                    for record in SeqIO.parse(fetch_handle, "genbank"):
                        """
                        Gets all data of the record.
                        """
                        source_feature = next((f for f in record.features if f.type == "source"), None)
                        gb_url = f"https://www.ncbi.nlm.nih.gov/{db}/{record.id}"

                        organism = record.annotations.get("organism", "Unknown")
                        nameStrainMatch = re.match(r"(.+?) (\d+\w*)", organism)
                        name = nameStrainMatch.group(1) if nameStrainMatch else "Unknown"
                        strain = nameStrainMatch.group(2) if nameStrainMatch else "Unknown"

                        genes = [f.qualifiers['gene'][0] for f in record.features if 'gene' in f.qualifiers]

                        metadata = {
                            "strain": strain,
                            "accession": record.name,
                            "gene": genes,
                            "organism": organism,
                            "description": record.description,
                            "location": source_feature.qualifiers.get("geo_loc_name", ["Unknown"])[0] if source_feature else "Unknown",
                            "isolation_source": source_feature.qualifiers.get("isolation_source", ["Unknown"])[0] if source_feature else "Unknown",
                            "sequence": str(record.seq),
                            "genbank_url": gb_url,
                        }

                        # Handle file splitting.
                        if split_field:
                            key = metadata.get(split_field, "Unknown")
                            if isinstance(key, list):
                                key = "_".join(key)
                            safe_key = re.sub(r"[^\w\-\.]", "_", str(key))
                            split_filename = f"{output_file.rsplit('.', 1)[0]}_{safe_key}" 
                            if key not in split_files:
                                split_files[key] = open(split_filename, "w")
                            out_fh = split_files[key]
                        else:
                            out_fh = main_fh

                        # Write output.
                        if json_mode:
                            json.dump(metadata, out_fh, indent=4) # Dump 'metadata' (list) as json into out_fh.
                            out_fh.write("\n")
                        else:
                            out_fh.write(f"(accession) Accession: {metadata['accession']}\n")
                            out_fh.write(f"(strain) Strain: {metadata['strain']}\n")
                            out_fh.write(f"(organism) Organism: {metadata['organism']}\n")
                            out_fh.write(f"(source) Isolation Source: {metadata['isolation_source']}\n")
                            out_fh.write(f"(location) Geo Location: {metadata['location']}\n")
                            out_fh.write(f"(gene) Gene(s): {', '.join(metadata['gene'])}\n")
                            out_fh.write(f"(description) Description: {metadata['description']}\n")
                            out_fh.write(f"(sequence) Sequence:\n{metadata['sequence']}\n")
                            out_fh.write(f"(url) GenBank URL: {metadata['genbank_url']}\n")
                            out_fh.write("\n" + "="*40 + "\n\n")

                        total_count += 1

                        if verbose_mode:
                            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                            print(f"[VERBOSE] Wrote record {total_count} of {len(ids)}: {metadata['accession']} ({name} - strain {metadata['strain']}) | {timestamp}")
                        else:
                            print(f"\rTotal fetched: {total_count}", end="", flush=True)

            # Close any open split files
            # e.g records_E6.txt, records_E7.txt
            for fh in split_files.values():
                fh.close()

            print(f"\nFinished writing {total_count} records\n")

        except KeyboardInterrupt:
            print(f"\nExiting..")
            break

        except HTTPError as e:
            print(f"""
HTTP Error encountered: {e}
This usually indicates malformed or missing API key/email.
Get your credentials at https://www.ncbi.nlm.nih.gov/account
""")
            continue

if __name__ == "__main__":
    main()
