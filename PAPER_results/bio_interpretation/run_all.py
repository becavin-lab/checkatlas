#!/usr/bin/env python3
"""
CheckAtlas Biological Interpretation — Main Orchestration Script.

Reads all benchmark results from /data/analysis/data_ganguly/checkatlas_files/
and generates publication-quality figures, tables, and interpretation text.

Usage:
    python run_all.py
"""

import sys
from pathlib import Path

OUT_DIR = Path("/home/ganguly/checkatlas/PAPER_results/bio_interpretation")


def main():
    print("=" * 70)
    print("  CheckAtlas — Biological Interpretation Pipeline")
    print("=" * 70)

    # Step 1: Aggregate all results
    print("\n[1/4] Aggregating results from TSV files...")
    from aggregator import aggregate_all
    data = aggregate_all()
    print(f"  Loaded: {len(data['df_all'])} metric values across {len(data['atlas_summary'])} atlases")

    # Print summary
    print("\n  Atlas Summary:")
    summary = data["atlas_summary"]
    for _, row in summary.iterrows():
        print(f"    {row['atlas']:<18} overall={row['overall_score']:.3f}  "
              f"C={row['cluster_score']:.3f}  A={row['annotation_score']:.3f}  "
              f"B={row['batch_score']:.3f}  D={row['dimred_score']:.3f}")

    # Step 2: Generate tables
    print("\n[2/4] Generating tables...")
    from tables import generate_all_tables
    generate_all_tables(data)

    # Step 3: Generate figures
    print("\n[3/4] Generating figures...")
    from figures import generate_all_figures
    generate_all_figures(data)

    # Step 4: Verify outputs
    print("\n[4/4] Verifying outputs...")
    fig_dir = OUT_DIR / "figures"
    png_files = sorted(fig_dir.glob("*.png"))
    pdf_files = sorted(fig_dir.glob("*.pdf"))
    tex_files = sorted(OUT_DIR.glob("*.tex"))

    for name, files in [("PNG", png_files), ("PDF", pdf_files), ("LaTeX", tex_files)]:
        if files:
            print(f"  {name}: {len(files)} files")
            for f in files:
                size = f.stat().st_size
                print(f"    {f.name} ({size:,} bytes)")

    print("\n" + "=" * 70)
    print("  Pipeline complete.")
    print(f"  Figures:  {fig_dir}/")
    print(f"  Tables:   {OUT_DIR}/")
    print("=" * 70)


if __name__ == "__main__":
    main()
