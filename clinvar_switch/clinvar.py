from urllib.request import urlretrieve
import gzip
import pandas as pd


def download_clinvar(date: str = "2021-11"):
    url = f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_{date}.txt.gz"
    dst = f"variant_summary_{date}.txt.gz"

    # urllib.request.ProxyHandler({"http": "myproxy:123"})
    urlretrieve(url, dst)


assembly_filter = "GRCh37"
significance_filter = ["Pathogenic", "Uncertain significance"]

df = []

with gzip.open("variant_summary_2022-11.txt.gz", "rt") as file:
    headers = next(file).strip().split("\t")

    for line in file:
        line = line.strip().split("\t")
        item = dict(zip(headers, line))

        assembly = item["Assembly"]
        chrom = item["Chromosome"]
        pos = item["PositionVCF"]
        ref = item["ReferenceAlleleVCF"]
        alt = item["AlternateAlleleVCF"]
        significance = item["ClinicalSignificance"]
        star = item["ClinSigSimple"]
        conditions = item["PhenotypeList"]
        variant_id = item["VariationID"]

        if assembly == assembly_filter and significance in significance_filter:

            df.append(
                {
                    "variant_id": variant_id,
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "significance": significance,
                    "star": star,
                    "conditions": conditions,
                }
            )

df = pd.DataFrame(df)

print(df)
