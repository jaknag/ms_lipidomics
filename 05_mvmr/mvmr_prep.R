#!/usr/bin/env python3
import pandas as pd
from pathlib import Path

def main():
    base = Path("/Users/nagrodzkij/cam/IMT/MS/lipidomics/extended")

    # ─── 1) Read in each file ───────────────────────────────────────────────────
    ms_path   = base / "severity" / "imsgc_mssev_discovery_hg38.tsv"
    bmi_path  = base / "BMI"      / "GCST90475153.h_filtered.tsv"
    vldl_path = base / "exposures" / "unfiltered" / "GCST90451306_VLDL-TG_filtered_mvmr.tsv"

    ms_dat   = pd.read_csv(ms_path,   sep="\t", dtype=str)
    bmi_dat  = pd.read_csv(bmi_path,  sep="\t", dtype=str)
    vldl_dat = pd.read_csv(vldl_path, sep="\t", dtype=str)

    # ─── 2) Pull out and rename the columns we care about ────────────────────
    # ms_dat uses "beta" and "SE"
    ms = (
        ms_dat
        .loc[:, ["SNP", "BETA", "SE"]]
        .rename(columns={"BETA": "beta_ms", "SE": "se_ms"})
        .set_index("SNP")
    )

    # bmi_dat and vldl_dat use "beta" and "standard_error"
    bmi = (
        bmi_dat
        .loc[:, ["SNP", "beta", "standard_error"]]
        .rename(columns={"beta": "beta_bmi", "standard_error": "se_bmi"})
    )

    vldl = (
        vldl_dat
        .loc[:, ["SNP", "beta", "standard_error"]]
        .rename(columns={"beta": "beta_vldl", "standard_error": "se_vldl"})
    )

    # ─── 3) Merge them together, keeping only SNPs present in all three ───────
    merged = (
        ms
        .merge(bmi,  on="SNP", how="inner")
        .merge(vldl, on="SNP", how="inner")
    )

    # ─── 4) (Optional) sort by SNP so the rows line up predictably ───────────
    merged.sort_values("SNP", inplace=True)

    # ─── 5) Write out the combined table ─────────────────────────────────────
    out_path = base / "mvmr" / "mvmr_bmi_vldl_raw.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"✅ Wrote {len(merged):,d} SNPs × {merged.shape[1]} columns to {out_path}")

if __name__ == "__main__":
    main()
