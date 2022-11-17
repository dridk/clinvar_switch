from urllib.request import urlretrieve
import gzip
import pandas as pd
import os 


def create_clinvar(date:str = "2021-11", significances = ["Pathogenic"])->pd.DataFrame:

    url = f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_{date}.txt.gz"
    dst = f"variant_summary_{date}.txt.gz"

    # urllib.request.ProxyHandler({"http": "myproxy:123"})
    
    if not os.path.exists(dst):
        print("Download file", dst)
        urlretrieve(url, dst)

    
    print("Read file ")
    assembly_filter = "GRCh37"
    vtype_filter = ["single nucleotide variant"]
    
    df = []
    
    with gzip.open(dst, "rt") as file:
        headers = next(file).strip().split("\t")
    
        for line in file:
            line = line.strip().split("\t")
            item = dict(zip(headers, line))
    
            assembly = item["Assembly"]
            vtype = item["Type"]
            chrom = item["Chromosome"]
            pos = item["PositionVCF"]
            ref = item["ReferenceAlleleVCF"]
            alt = item["AlternateAlleleVCF"]
            signif = item["ClinicalSignificance"].replace(" ","_").lower()
            star = item["ClinSigSimple"]
            conditions = item["PhenotypeList"]
            variant_id = item["VariationID"]
    
            if assembly == assembly_filter and signif in significances and vtype in vtype_filter:
    
                df.append(
                    {
                        "variant_id": variant_id,
                        "chrom": chrom,
                        "pos": int(pos),
                        "ref": ref,
                        "alt": alt,
                        "significance": signif,
                        "star": int(star),
                        "conditions": conditions,
                    }
                )
    
    df = pd.DataFrame(df)

    # Sort chrom / position 
    chrom_order = [str(i) for i in range(1,23)] + ["X","Y","MT"]
    df["chrom"]  = pd.Categorical(df["chrom"], chrom_order)
    df = df.sort_values(["chrom","pos"])

    return df


df1 = create_clinvar("2021-11", ["uncertain_significance"])
df2 = create_clinvar("2022-11", ["pathogenic","likely_pathogenic"])






df1.to_parquet("a.parquet")
df2.to_parquet("b.parquet")



