#!/usr/bin/env python3
"""
Get PMIDs for a keyword (via PubMed esearch), and optionally
filter to PMIDs that PubTator annotates with chosen concepts.

Examples:
  python pmids_from_pubtator.py "alkaptonuria"
  python pmids_from_pubtator.py "alkaptonuria" --concepts disease,chemical --out pmids.txt
"""

import argparse, time, sys
from typing import List, Set, Iterable
import requests
from requests.adapters import HTTPAdapter, Retry

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PT_EXPORT = "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson"

def session():
    s = requests.Session()
    retries = Retry(total=5, backoff_factor=0.5, status_forcelist=(429,500,502,503,504), allowed_methods={"GET"})
    s.mount("https://", HTTPAdapter(max_retries=retries))
    s.headers.update({"User-Agent": "PubTatorPMIDFetcher/2.0"})
    return s

def esearch_all_pmids(query: str, chunk: int = 10000) -> List[str]:
    """Collect all PMIDs for a PubMed query via esearch."""
    pmids: List[str] = []
    with session() as s:
        # first get count
        params = {"db":"pubmed","retmode":"json","term":query,"retmax":0}
        r = s.get(EUTILS, params=params, timeout=60); r.raise_for_status()
        count = int(r.json()["esearchresult"]["count"])
        for start in range(0, count, chunk):
            params = {"db":"pubmed","retmode":"json","term":query,"retmax":chunk,"retstart":start}
            r = s.get(EUTILS, params=params, timeout=60); r.raise_for_status()
            pmids.extend(r.json()["esearchresult"].get("idlist", []))
            time.sleep(0.34)
    # unique & numeric sort
    return sorted(set(pmids), key=lambda x: int(x))

def pubtator_has_concept(pmids: List[str], concepts: List[str], page: int = 1000) -> Set[str]:
    """Return subset of PMIDs that PubTator returns with ≥1 annotation for the given concepts."""
    keep: Set[str] = set()
    with session() as s:
        for i in range(0, len(pmids), page):
            batch = pmids[i:i+page]
            params = {"pmids": ",".join(batch), "concepts": ",".join(concepts)}
            r = s.get(PT_EXPORT, params=params, timeout=120); r.raise_for_status()
            data = r.json()
            for doc in data.get("documents", []):
                ann_count = 0
                for p in doc.get("passages", []):
                    ann_count += len(p.get("annotations", []))
                if ann_count > 0:
                    keep.add(doc.get("id"))
            time.sleep(0.5)
    return keep

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("keyword", help="PubMed query string (e.g., 'alkaptonuria' or full PubMed syntax).")
    ap.add_argument("--concepts", help="Comma-separated PubTator concepts to filter by (disease,gene,chemical,mutation,species,cellline).", default="")
    ap.add_argument("--out", help="Write PMIDs to file; else print to stdout.", default="")
    args = ap.parse_args()

    try:
        pmids = esearch_all_pmids(args.keyword)
        if args.concepts.strip():
            concepts = [c.strip() for c in args.concepts.split(",") if c.strip()]
            subset = pubtator_has_concept(pmids, concepts)
            pmids = sorted(subset, key=lambda x: int(x))
    except requests.HTTPError as e:
        sys.stderr.write(f"HTTP error: {e}\n"); sys.exit(1)
    except requests.RequestException as e:
        sys.stderr.write(f"Request failed: {e}\n"); sys.exit(1)

    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            for p in pmids: f.write(p+"\n")
        print(f"Wrote {len(pmids)} PMIDs to {args.out}")
    else:
        for p in pmids: print(p)
        sys.stderr.write(f"\nTotal PMIDs: {len(pmids)}\n")

if __name__ == "__main__":
    main()
