## 💻 Technical structure
### Project Schema
```text
├── sources/                      # Raw NHANES 2021-2022 .XPT files
├── notebooks/                    # Step-by-step documentation of the Google framework
│   ├── 1_prepare_process.ipynb   # ETL Pipeline: Batch-loading, left-joins, missing data audit
│   ├── 2_analyze_modeling.ipynb  # Feature engineering, PR calculation
│   └── 3_share_visuals.ipynb     # Code-based plots (Seaborn/Matplotlib) exported to presentation
├── dashboards/                   # Optional: Files or templates for Tableau
│   └── obesity_treatment_gap.twbx
├── .gitignore             # Configured with wildcards to block data leakage (*.XPT, *.csv, sources/)
├── README.md              # The core asset: Your business and clinical narrative
├── requirements.txt       # Python environment dependencies (pandas, pyreadstat, statsmodels, etc.)
└── LICENSE                # Open-source license (e.g., MIT)
```