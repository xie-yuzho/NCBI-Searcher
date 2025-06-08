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


# This is the GUI version, if you are not familiar with tkinter i don't recommend modifying this.



import tkinter as tk
from tkinter import scrolledtext, messagebox
from Bio import Entrez, SeqIO
import re
from datetime import datetime
from urllib.error import HTTPError
import json
from contextlib import contextmanager
from itertools import cycle
import threading


@contextmanager
def dummy_context_manager():
    yield None

def parse_credentials(cred_string):
    pairs = [p.strip() for p in cred_string.split(",")]
    return [tuple(p.split(":")) for p in pairs if ":" in p]

def clean_query_flags(query, flag):
    pattern = None
    flag_value = None
    if flag in ['-o', '-s', '-c', '-b']:
        pattern = re.compile(rf"{flag}\s+([^\s]+)")
        match = pattern.search(query)
        if match:
            flag_value = match.group(1)
            query = pattern.sub("", query).strip()
    else:
        if flag in query:
            flag_value = True
            query = query.replace(flag, "").strip()
    return query, flag_value

class NCBIFetcherApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("NCBI Fetcher GUI")
        self.geometry("700x500")

        # Input fields
        tk.Label(self, text="NCBI Database (e.g., protein, nuccore):").pack(pady=5)
        self.db_entry = tk.Entry(self, width=50)
        self.db_entry.pack()

        tk.Label(self, text="Search Query (use flags as needed):").pack(pady=5)
        self.query_entry = tk.Entry(self, width=50)
        self.query_entry.pack()

        self.pause_event = threading.Event()  # When set, fetching is paused
        self.stop_event = threading.Event()   # When set, fetching should stop

        frame = tk.Frame(self)
        frame.pack(pady=5)

        self.query_builder_button = tk.Button(self, text="Open Query Builder", command=self.open_query_builder)
        self.query_builder_button.pack(pady=10)

        self.fetch_button = tk.Button(frame, text="Fetch", command=self.start_fetch_thread)
        self.fetch_button.pack(side=tk.LEFT, padx=5)

        self.pause_button = tk.Button(frame, text="Pause", command=self.toggle_pause)
        self.pause_button.pack(side=tk.LEFT, padx=5)

        self.stop_button = tk.Button(frame, text="Stop", command=self.stop_fetch)
        self.stop_button.pack(side=tk.LEFT, padx=5)

        # Output display box
        self.output_box = scrolledtext.ScrolledText(self, width=80, height=20)
        self.output_box.pack(padx=10, pady=10)

    def open_query_builder(self):
        qb = tk.Toplevel(self)
        qb.title("Query Builder")
        qb.geometry("500x700")

        # Canvas + scrollbar
        canvas = tk.Canvas(qb, borderwidth=0)
        scrollbar = tk.Scrollbar(qb, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        content_frame = tk.Frame(canvas)
        canvas.create_window((0, 0), window=content_frame, anchor="nw")
        content_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))

        # Main term
        tk.Label(content_frame, text="Main term(s):").pack(anchor="w", padx=10, pady=(10,0))
        main_term_entry = tk.Entry(content_frame, width=40)
        main_term_entry.pack(padx=10)

        # Logic selectors
        inner_logic_var = tk.StringVar(value="OR")
        outer_logic_var = tk.StringVar(value="OR")
        tk.Label(content_frame, text="Inner logic in each block:").pack(anchor="w", padx=15, pady=(10,0))
        for val in ("OR", "AND"):
            tk.Radiobutton(content_frame, text=val, variable=inner_logic_var, value=val).pack(anchor="w", padx=20)
        tk.Label(content_frame, text="Logic between blocks:").pack(anchor="w", padx=10, pady=(10,0))
        for val in ("OR", "AND"):
            tk.Radiobutton(content_frame, text=val, variable=outer_logic_var, value=val).pack(anchor="w", padx=20)

        # Inclusion Block 1
        tk.Label(content_frame, text="Must include block 1:").pack(anchor="w", padx=10, pady=(10,0))
        include_frame1 = tk.Frame(content_frame)
        include_frame1.pack(padx=10, fill="x")
        include_entries1 = []
        def add_include1(text=""):
            row = tk.Frame(include_frame1)
            row.pack(fill="x", pady=2)
            e = tk.Entry(row, width=20)
            e.pack(side="left", padx=5)
            te = tk.Entry(row, width=20)
            te.pack(side="right", padx=5)
            include_entries1.append((e,te))
        add_include1()
        tk.Button(content_frame, text="Add to block 1", command=add_include1).pack(pady=5)

        # Inclusion Block 2
        tk.Label(content_frame, text="Must include block 2:").pack(anchor="w", padx=10, pady=(10,0))
        include_frame2 = tk.Frame(content_frame)
        include_frame2.pack(padx=10, fill="x")
        include_entries2 = []
        def add_include2(text=""):
            row = tk.Frame(include_frame2)
            row.pack(fill="x", pady=2)
            e = tk.Entry(row, width=20)
            e.pack(side="left", padx=5)
            te = tk.Entry(row, width=20)
            te.pack(side="right", padx=5)
            include_entries2.append((e,te))
        add_include2()
        tk.Button(content_frame, text="Add to block 2", command=add_include2).pack(pady=20)

        # Exclusion block 1
        excl_inner_logic_var = tk.StringVar(value="OR")
        excl_outer_logic_var = tk.StringVar(value="OR")
        tk.Label(content_frame, text="Inner logic in each block:").pack(anchor="w", padx=15, pady=(10,0))
        for val in ("OR", "AND"):
            tk.Radiobutton(content_frame, text=val, variable=excl_inner_logic_var, value=val).pack(anchor="w", padx=20)
        tk.Label(content_frame, text="Logic between blocks:").pack(anchor="w", padx=10, pady=(10,0))
        for val in ("OR", "AND"):
            tk.Radiobutton(content_frame, text=val, variable=excl_outer_logic_var, value=val).pack(anchor="w", padx=20)

        tk.Label(content_frame, text="Exclude block 1:").pack(anchor="w", padx=10, pady=(15,0))
        exclude_frame1 = tk.Frame(content_frame)
        exclude_frame1.pack(padx=10, fill="x")
        exclude_entries1 = []

        def add_exclude1(term="", typ=""):
            row = tk.Frame(exclude_frame1)
            row.pack(fill="x", pady=2)
            term_e = tk.Entry(row, width=20)
            term_e.insert(0, term)
            term_e.pack(side="left", padx=5)
            type_e = tk.Entry(row, width=20)
            type_e.insert(0, typ)
            type_e.pack(side="right", padx=5)
            exclude_entries1.append((term_e, type_e))

        add_exclude1()
        tk.Button(content_frame, text="Add to exclude block 1", command=add_exclude1).pack(pady=5)

        # Exclusion block 2
        tk.Label(content_frame, text="Exclude block 2:").pack(anchor="w", padx=10, pady=(10,0))
        exclude_frame2 = tk.Frame(content_frame)
        exclude_frame2.pack(padx=10, fill="x")
        exclude_entries2 = []

        def add_exclude2(term="", typ=""):
            row = tk.Frame(exclude_frame2)
            row.pack(fill="x", pady=2)
            term_e = tk.Entry(row, width=20)
            term_e.insert(0, term)
            term_e.pack(side="left", padx=5)
            type_e = tk.Entry(row, width=20)
            type_e.insert(0, typ)
            type_e.pack(side="right", padx=5)
            exclude_entries2.append((term_e, type_e))

        add_exclude2()
        tk.Button(content_frame, text="Add to exclude block 2", command=add_exclude2).pack(pady=5)

        # Flags etc.
        # This is very hurtful to look at.
        flags = { '-o': {'label':'Output file (-o)', 'var':tk.StringVar()},
                  '-s': {'label':'Split field (-s)','var':tk.StringVar()},
                  '-b': {'label':'Batch size (-b)','var':tk.StringVar()},
                  '-c': {'label':'Credentials (-c)','var':tk.StringVar()} }
        bool_flags = {'-j':{'label':'JSON mode (-j)','var':tk.BooleanVar()} }
        vf_frame = tk.LabelFrame(content_frame, text="Value flags", padx=10, pady=10)
        vf_frame.pack(fill="x", padx=10, pady=5)
        for f,info in flags.items():
            row = tk.Frame(vf_frame)
            row.pack(fill="x", pady=2)
            tk.Label(row,text=info['label'],width=20,anchor='w').pack(side='left')
            tk.Entry(row,textvariable=info['var'],width=20).pack(side='right')
        bf_frame = tk.LabelFrame(content_frame, text="Toggle flags", padx=10, pady=10)
        bf_frame.pack(fill="x", padx=10, pady=5)
        for f,info in bool_flags.items():
            tk.Checkbutton(bf_frame,text=info['label'],variable=info['var']).pack(anchor='w')

        def apply_query():
            main_term = main_term_entry.get().strip()
            if not main_term:
                messagebox.showwarning("Warning","Main search term cannot be empty")
                return

            b1 = [] # Handling include entries 1
            for term_entry, type_entry in include_entries1:
                term = term_entry.get().strip()
                typ = type_entry.get().strip()
                if term:
                    if typ:
                        b1.append(f"{term}[{typ}]")
                    else:
                        b1.append(term)

            b2 = [] # Handling include entries 2
            for term_entry, type_entry in include_entries2:
                term = term_entry.get().strip()
                typ = type_entry.get().strip()
                if term:
                    if typ:
                        b2.append(f"{term}[{typ}]")
                    else:
                        b2.append(term)

            inner_op = f" {inner_logic_var.get()} "
            outer_op = f" {outer_logic_var.get()} "
            subqs = [] # Combined terms and blocks
            for terms in (b1, b2):
                if terms:
                    subqs.append("(" + inner_op.join(terms) + ")")

            query_parts = [main_term]
            if subqs:
                combined = outer_op.join(subqs)
                query_parts.append("AND (" + combined + ")")

            exclude_blocks = []
            for exclude_entries in [exclude_entries1, exclude_entries2]:
                block = []
                for term_entry, type_entry in exclude_entries:
                    term = term_entry.get().strip()
                    typ = type_entry.get().strip()
                    if term:
                        if typ:
                            block.append(f"{term}[{typ}]")
                        else:
                            block.append(term)
                if block:
                    exclude_blocks.append("(" + " OR ".join(block) + ")")

            if exclude_blocks:
                full_exclude = " OR ".join(exclude_blocks)
                query_parts.append("NOT (" + full_exclude + ")")

            # Handling flags
            for f,info in bool_flags.items():
                if info['var'].get(): query_parts.append(f)
            for f,info in flags.items():
                v = info['var'].get().strip()
                if v: query_parts += [f,v]

            final = " ".join(query_parts)
            self.query_entry.delete(0,tk.END)
            self.query_entry.insert(0,final)
            qb.destroy()

        tk.Button(content_frame,text="Apply Query",command=apply_query).pack(pady=5)





    def toggle_pause(self):
        if not self.pause_event.is_set():
            self.pause_event.set()
            self.pause_button.config(text="Resume")
            self.log("Fetching paused.")
        else:
            self.pause_event.clear()
            self.pause_button.config(text="Pause")
            self.log("Fetching resumed.")

    def stop_fetch(self):
        self.stop_event.set()
        self.log("Stopping fetch operation...")


    def log(self, message):
        self.output_box.insert(tk.END, message + "\n")
        self.output_box.see(tk.END)

    def start_fetch_thread(self):
        # Disable button to prevent multiple clicks
        self.fetch_button.config(state=tk.DISABLED)
        self.output_box.delete(1.0, tk.END)
        threading.Thread(target=self.run_fetch).start()

    def run_fetch(self):
        try:
            db = self.db_entry.get().strip().lower()
            query = self.query_entry.get().strip()

            if db in {"exit", "x"}:
                self.log("Exiting...")
                self.fetch_button.config(state=tk.NORMAL)
                return

            if not db:
                self.log("Database input cannot be empty.")
                self.fetch_button.config(state=tk.NORMAL)
                return

            if not query:
                self.log("Search query cannot be empty.")
                self.fetch_button.config(state=tk.NORMAL)
                return

            # Default settings
            output_file = "records.txt"
            total_count = 0
            split_field = None
            verbose_mode = True  # For GUI, verbose will be on.
            json_mode = False
            batch_size = 50
            credentials = []
            credential_cycle = cycle([])

            # Parse flags in query
            query, output_file_opt = clean_query_flags(query, "-o")
            if output_file_opt:
                output_file = output_file_opt

            query, split_field = clean_query_flags(query, "-s")
            query, verbose_mode_flag = clean_query_flags(query, "-v")
            verbose_mode = verbose_mode_flag or verbose_mode
            query, json_mode_flag = clean_query_flags(query, "-j")
            json_mode = bool(json_mode_flag)
            query, batch_size_str = clean_query_flags(query, "-b")
            batch_size = int(batch_size_str) if batch_size_str else 50
            query, cred_string = clean_query_flags(query, "-c")
            if cred_string:
                credentials = parse_credentials(cred_string)
                credential_cycle = cycle(credentials)

            if verbose_mode:
                self.log(f"[VERBOSE] Searching '{db}' for: '{query}'")
                if credentials:
                    self.log(f"[VERBOSE] Using {len(credentials)} credential(s) for API access.")
                self.log(f"[VERBOSE] Output file: {output_file}")
                if split_field:
                    self.log(f"[VERBOSE] Splitting output files by field: {split_field}")
                self.log(f"[VERBOSE] Output format: {'JSON' if json_mode else 'Plain text'}")

            try:
                api_key, email = next(credential_cycle)
                Entrez.email = email
                Entrez.api_key = api_key
            except StopIteration:
                Entrez.email = None
                Entrez.api_key = None
                if verbose_mode:
                    self.log("[WARN] No credentials provided. Get API key at https://www.ncbi.nlm.nih.gov/account")

            handle = Entrez.esearch(db=db, term=query, retmax=50000)
            record = Entrez.read(handle)
            ids = record.get("IdList", [])
            self.log(f"Found {len(ids)} record(s).")

            if not ids:
                self.log("No results found.")
                self.fetch_button.config(state=tk.NORMAL)
                return

            split_files = {}
            with open(output_file, "w") if not split_field else dummy_context_manager() as main_fh:
                for batch_start in range(0, len(ids), batch_size):
                    # Check if stop requested
                    if self.stop_event.is_set():
                        self.log("Fetch operation stopped by user.")
                        break

                    # Check if paused; wait here until resumed
                    # If paused, the batch will be stopped until resumed.
                    while self.pause_event.is_set():
                        if self.stop_event.is_set():
                            self.log("Fetch operation stopped by user during pause.")
                            break
                        threading.Event().wait(0.1)

                    # If stopped during pause
                    if self.stop_event.is_set():
                        break
                    try:
                        api_key, email = next(credential_cycle)
                        Entrez.email = email
                        Entrez.api_key = api_key
                    except StopIteration:
                        Entrez.email = None
                        Entrez.api_key = None

                    batch_ids = ids[batch_start:batch_start + batch_size]

                    if verbose_mode:
                        self.log(f"[VERBOSE] Fetching batch {batch_start // batch_size + 1}: {len(batch_ids)} records using credential: {email if Entrez.email else 'None'}")

                    fetch_handle = Entrez.efetch(db=db, id=",".join(batch_ids), rettype="gb", retmode="text")

                    for record in SeqIO.parse(fetch_handle, "genbank"):
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

                        if split_field:
                            key = metadata.get(split_field, "Unknown")
                            if isinstance(key, list):
                                key = "_".join(key)
                            safe_key = re.sub(r"[^\w\-\.]", "_", str(key))
                            split_filename = f"{output_file.rsplit('.', 1)[0]}_{safe_key}.{output_file.rsplit('.', 1)[1]}" 
                            if key not in split_files:
                                split_files[key] = open(split_filename, "w")
                            out_fh = split_files[key]
                        else:
                            out_fh = main_fh

                        if json_mode:
                            json.dump(metadata, out_fh, indent=4)
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

                        self.log(f"Wrote record {total_count} of {len(ids)}: {metadata['accession']}")

            for fh in split_files.values():
                fh.close()

            self.log(f"\nFinished writing {total_count} records\n")

        except KeyboardInterrupt:
            self.log("\nExiting..")
        except HTTPError as e:
            self.log(f"HTTP Error encountered: {e}\nCheck API keys and emails.")
        except Exception as e:
            self.log(f"Unexpected error: {e}")
        finally:
            self.fetch_button.config(state=tk.NORMAL)

if __name__ == "__main__":
    app = NCBIFetcherApp()
    app.mainloop()
