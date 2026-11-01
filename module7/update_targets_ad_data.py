"""Build the Module 7 gene table from Module 2/3 outputs and public APIs.

The Module 2 RNA and protein columns are absolute tumor-vs-normal effect
magnitudes, not signed fold changes.

Ensembl gene IDs and GWAS p-values are retrieved from the GWAS Catalog REST API
v2. The p-value is the smallest curated association for that mapped gene and
colorectal cancer (MONDO_0005575). No association means the GWAS fields stay
blank; it is not interpreted as a p-value of 1.
"""

import csv
import json
from pathlib import Path
from time import sleep
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]
MODULE2_FILE = ROOT / "module2" / "targets_ad.csv"
MODULE3_FILE = ROOT / "module3" / "targets_cc_w3.csv"
OUTPUT_FILE = ROOT / "module7" / "data" / "targets_ad_updated_data.csv"
GWAS_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2"
GWAS_TRAIT_ID = "MONDO_0005575"
GWAS_TRAIT_NAME = "colorectal cancer"
USER_AGENT = "BIOT6900-Module7-data-builder/1.0"

# Common targeted agents for metastatic colorectal cancer. EGFR is the target
# of both anti-EGFR antibodies; VEGFA encodes the VEGF-A ligand targeted by
# bevacizumab. These are supplemental target genes, not Module 2/3 candidates.
DRUG_TARGET_GENES = {
    "VEGFA": {
        "drugs": "Bevacizumab (Avastin)",
        "sources": "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=939b5d1f-9fb2-4499-80ef-0607aa6b114e",
    },
    "EGFR": {
        "drugs": "Cetuximab (Erbitux); Panitumumab (Vectibix)",
        "sources": (
            "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=8bc6397e-4bd8-4d37-a007-a327e4da34d9; "
            "https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=e0fa4bca-f245-4d92-ae29-b0c630a315c2"
        ),
    },
}

OUTPUT_COLUMNS = [
    "gene",
    "ensembl",
    "gwas_p",
    "rna_lfc",
    "prot_lfc",
    "druggability",
    "gwas_reported_trait",
    "gwas_accession",
    "gwas_variant",
    "gwas_association_url",
    "ensembl_source_url",
    "associated_drugs",
    "drug_source_urls",
]


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError(f"Required input file not found: {path}")
    with path.open(newline="", encoding="utf-8-sig") as file:
        reader = csv.DictReader(file)
        if reader.fieldnames is None:
            raise ValueError(f"Input CSV has no header: {path}")
        return reader.fieldnames, list(reader)


def gene_column(fieldnames: list[str], path: Path) -> str:
    for candidate in ("gene", "symbol", "Unnamed: 0", ""):
        if candidate in fieldnames:
            return candidate
    raise ValueError(f"No gene-symbol column found in {path}; columns: {fieldnames}")


def get_json(url: str, allow_not_found: bool = False) -> dict | None:
    request = Request(
        url,
        headers={
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        },
    )
    try:
        with urlopen(request, timeout=30) as response:
            return json.load(response)
    except HTTPError as error:
        if allow_not_found and error.code in (400, 404):
            return None
        raise RuntimeError(f"Online data request failed for {url}: {error}") from error
    except (URLError, TimeoutError) as error:
        raise RuntimeError(f"Online data request failed for {url}: {error}") from error


def gwas_gene_ensembl_ids(gene: str) -> tuple[str, str]:
    url = f"{GWAS_URL}/genes/{quote(gene)}"
    result = get_json(url, allow_not_found=True)
    if result is None:
        return "", url
    return ";".join(result.get("ensembl_gene_ids", [])), url


def gwas_associations(gene: str) -> tuple[dict, str]:
    query = urlencode(
        {
            "efo_id": GWAS_TRAIT_ID,
            "mapped_gene": gene,
            "size": 100,
        }
    )
    base_url = f"{GWAS_URL}/associations"
    associations = []
    page = 0

    while True:
        url = f"{base_url}?{query}&page={page}"
        result = get_json(url)
        if result is None:
            raise RuntimeError(f"Unexpected empty API response from {url}")
        associations.extend(
            result.get("_embedded", {}).get("associations", [])
        )
        page_info = result.get("page", {})
        if page + 1 >= page_info.get("totalPages", 0):
            break
        page += 1
        sleep(0.1)

    matching = []
    for association in associations:
        mapped_genes = association.get("mapped_genes", [])
        p_value = association.get("p_value")
        if (
            gene.casefold() in {mapped.casefold() for mapped in mapped_genes}
            and p_value is not None
        ):
            try:
                matching.append((float(p_value), association))
            except (TypeError, ValueError) as error:
                association_id = association.get("association_id", "unknown")
                raise ValueError(
                    f"Invalid p-value for GWAS association {association_id}: "
                    f"{p_value!r}"
                ) from error

    if not matching:
        return {}, f"{base_url}?{query}"
    _, best = min(matching, key=lambda match: match[0])
    return best, f"{base_url}?{query}"


def main() -> None:
    module2_fields, module2_rows = read_csv(MODULE2_FILE)
    module3_fields, module3_rows = read_csv(MODULE3_FILE)

    required_module2 = {"transcriptomic", "proteomic"}
    missing_module2 = required_module2.difference(module2_fields)
    if missing_module2:
        raise ValueError(
            f"{MODULE2_FILE} is missing required columns: "
            f"{', '.join(sorted(missing_module2))}"
        )
    if "p_sm" not in module3_fields:
        raise ValueError(f"{MODULE3_FILE} is missing required column: p_sm")

    module2_gene = gene_column(module2_fields, MODULE2_FILE)
    module3_gene = gene_column(module3_fields, MODULE3_FILE)

    tractability: dict[str, str] = {}
    for row in module3_rows:
        gene = row.get(module3_gene, "").strip()
        if gene:
            if gene in tractability:
                raise ValueError(f"Duplicate gene {gene!r} in {MODULE3_FILE}")
            tractability[gene] = row.get("p_sm", "").strip()

    output_rows = []
    seen_genes = set()
    for row in module2_rows:
        gene = row.get(module2_gene, "").strip()
        if not gene:
            continue
        if gene in seen_genes:
            raise ValueError(f"Duplicate gene {gene!r} in {MODULE2_FILE}")
        seen_genes.add(gene)
        ensembl_id, ensembl_source_url = gwas_gene_ensembl_ids(gene)
        association, gwas_source_url = gwas_associations(gene)
        if not ensembl_id:
            print(f"Warning: Ensembl returned no ID for {gene}.")
        if not association:
            print(f"Note: no GWAS Catalog {GWAS_TRAIT_NAME} association for {gene}.")

        gwas_p = association.get("p_value")
        output_rows.append(
            {
                "gene": gene,
                "ensembl": ensembl_id or row.get("ensembl", "").strip(),
                "gwas_p": str(gwas_p) if gwas_p is not None else "",
                "rna_lfc": row.get("transcriptomic", "").strip(),
                "prot_lfc": row.get("proteomic", "").strip(),
                "druggability": tractability.get(gene, ""),
                "gwas_reported_trait": "; ".join(
                    association.get("reported_trait", [])
                ),
                "gwas_accession": association.get("accession_id", ""),
                "gwas_variant": "; ".join(
                    allele.get("rs_id", "")
                    for allele in association.get("snp_allele", [])
                ),
                "gwas_association_url": association.get("_links", {})
                .get("self", {})
                .get("href", gwas_source_url if not association else ""),
                "ensembl_source_url": ensembl_source_url,
                "associated_drugs": "",
                "drug_source_urls": "",
            }
        )

    for gene, drug_info in DRUG_TARGET_GENES.items():
        if gene in seen_genes:
            raise ValueError(
                f"Supplemental drug-target gene {gene!r} is already in "
                f"{MODULE2_FILE}; review before merging these records."
            )
        ensembl_id, ensembl_source_url = gwas_gene_ensembl_ids(gene)
        association, gwas_source_url = gwas_associations(gene)
        if not ensembl_id:
            raise ValueError(f"No Ensembl gene ID returned for drug target {gene}.")
        if association:
            print(
                f"Note: a colorectal-cancer GWAS association was found for "
                f"supplemental gene {gene}; retaining its p-value."
            )
        else:
            print(
                f"Note: no GWAS Catalog {GWAS_TRAIT_NAME} association for "
                f"supplemental gene {gene}."
            )
        output_rows.append(
            {
                "gene": gene,
                "ensembl": ensembl_id,
                "gwas_p": (
                    str(association["p_value"]) if association else ""
                ),
                "rna_lfc": "",
                "prot_lfc": "",
                "druggability": "",
                "gwas_reported_trait": "; ".join(
                    association.get("reported_trait", [])
                ),
                "gwas_accession": association.get("accession_id", ""),
                "gwas_variant": "; ".join(
                    allele.get("rs_id", "")
                    for allele in association.get("snp_allele", [])
                ),
                "gwas_association_url": association.get("_links", {})
                .get("self", {})
                .get("href", gwas_source_url if not association else ""),
                "ensembl_source_url": ensembl_source_url,
                "associated_drugs": drug_info["drugs"],
                "drug_source_urls": drug_info["sources"],
            }
        )

    if not output_rows:
        raise ValueError(f"No gene rows found in {MODULE2_FILE}")

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_FILE.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(output_rows)

    missing_tractability = sum(not row["druggability"] for row in output_rows)
    missing_ensembl = sum(not row["ensembl"] for row in output_rows)
    missing_gwas = sum(not row["gwas_p"] for row in output_rows)
    print(f"Wrote {len(output_rows)} genes to {OUTPUT_FILE}")
    print(
        f"GWAS p-values are the smallest curated association per mapped gene "
        f"for {GWAS_TRAIT_NAME} ({GWAS_TRAIT_ID})."
    )
    print(
        "Note: rna_lfc and prot_lfc contain Module 2 absolute effect magnitudes "
        "(direction was discarded)."
    )
    print(
        f"Blank values: ensembl={missing_ensembl}, gwas_p={missing_gwas}, "
        f"druggability={missing_tractability}."
    )


if __name__ == "__main__":
    main()
