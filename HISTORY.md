Changelog
=========


(unreleased)
------------
- Minor bug correction. [drbecavin]
- Merge pull request #26 from becavin-lab/rpy2. [drbecavin]

  Rpy2
- Minor bug correction. [drbecavin]
- Fix requirements. [drbecavin]
- Fix requirements. [drbecavin]
- Fix requirements. [drbecavin]
- Run lint, clean code. [drbecavin]
- Add metrics calculation for Seurat objects (dimred not working yet)
  [drbecavin]
- Implement Seurat management via rpy2. Summary, adata, qc and umap/tsne
  is now performed with R. No conversion from seurat to scanpy is needed
  anymore. [drbecavin]
- Parse Cellranger files by searching for clustering and reductions in
  outs/analysis/ folder. [drbecavin]
- Separate search of scanpy, seurat and cellranger. Add atlas_seurat for
  managing Seurat object with rpy2. [drbecavin]
- Update README.md. [drbecavin]


0.0.16 (2022-07-26)
-------------------
- Release: version 0.0.16 🚀 [drbecavin]
- Improve metrics management. Add annot and dimred metrics. [drbecavin]
- Add Metrics management with arguments, you can specify now the metrics
  to run inside the argument list or usng config file. Fix also the use
  of config/default_config.yaml to specify program arguments.
  [drbecavin]


0.0.15 (2022-07-22)
-------------------
- Release: version 0.0.15 🚀 [drbecavin]
- Replace print by logger.info and logger.debug. [drbecavin]


0.0.14 (2022-07-22)
-------------------
- Release: version 0.0.14 🚀 [drbecavin]
- Add config file management for program arguments. [drbecavin]
- Add parameters and logger in main. [drbecavin]


0.0.13 (2022-07-19)
-------------------
- Release: version 0.0.13 🚀 [drbecavin]
- Add QC tables for Knee plots. [drbecavin]


0.0.12 (2022-07-08)
-------------------
- Release: version 0.0.12 🚀 [drbecavin]


0.0.11 (2022-07-08)
-------------------
- Release: version 0.0.11 🚀 [drbecavin]


0.0.10 (2022-07-08)
-------------------
- Release: version 0.0.10 🚀 [drbecavin]
- Prepare QC table production. [drbecavin]
- Test.txt. [drbecavin]
- Change in workflow. [drbecavin]


0.0.9 (2022-06-20)
------------------
- Release: version 0.0.9 🚀 [drbecavin]
- Fix project release. [drbecavin]


0.0.8 (2022-06-20)
------------------
- Release: version 0.0.8 🚀 [drbecavin]
- Fix project release. [drbecavin]


0.0.7 (2022-06-20)
------------------
- Release: version 0.0.7 🚀 [drbecavin]
- Release: version 0.0.6 🚀 [drbecavin]


0.0.6 (2022-06-20)
------------------
- Release: version 0.0.6 🚀 [drbecavin]
- Clean docs. [drbecavin]
- Update README.md. [drbecavin]
- Modify readme docs. [drbecavin]
- Update README.md. [drbecavin]


0.0.5 (2022-06-20)
------------------
- Release: version 0.0.5 🚀 [drbecavin]


0.0.4 (2022-06-20)
------------------
- Release: version 0.0.4 🚀 [drbecavin]
- Modify readme docs. [drbecavin]
- Update README.md. [drbecavin]
- Add examples in the doc. [drbecavin]


0.0.3 (2022-06-18)
------------------
- Release: version 0.0.3 🚀 [drbecavin]
- Merge branch 'main' of github.com:becavin-lab/checkatlas good.
  [drbecavin]
- Update README.md. [drbecavin]
- Update README.md. [drbecavin]
- Fix umap display. [drbecavin]
- Fix category management in adata.obs. [drbecavin]
- Add test on scanp version. [drbecavin]


0.0.2 (2022-06-16)
------------------
- Release: version 0.0.2 🚀 [drbecavin]
- Release: version 0.0.1 🚀 [drbecavin]
- Fix release making. [drbecavin]
- Merge branch 'main' of github.com:becavin-lab/checkatlas Good.
  [drbecavin]
- Update README.md. [drbecavin]
- Fix release making. [drbecavin]


0.0.1 (2022-06-16)
------------------
- Release: version 0.0.1 🚀 [drbecavin]
- Release: version  🚀 [drbecavin]
- Release: version  🚀 [drbecavin]
- Release: version  🚀 [drbecavin]
- Release: version  🚀 [drbecavin]
- Release: version  🚀 [drbecavin]
- Clean project before release. [drbecavin]
- Release: version  🚀 [drbecavin]
- Add corrupted h5ad management with AnnDataReadError. [drbecavin]
- Add corrupted h5ad management with AnnDataReadError. [drbecavin]
- Add resume capacity in checkatlas.py. [drbecavin]
- Fix obs_beys selection. [drbecavin]
- Add list_atlases.csv. [drbecavin]
- Fix list files. [drbecavin]
- Add HLCA index names. [drbecavin]
- Release: version  🚀 [drbecavin]
- Clean project. [drbecavin]
- Clean files. [drbecavin]
- Remove rds from type of file. [christophe Becavin]
- Prepare checkatlas for discovair. [christophe Becavin]
- Clean all codes. make fmt, make lint, and make docs. [christophe
  Becavin]
- Make fmt. [christophe Becavin]
- Add read the docs conf. [christophe Becavin]
- Add read the docs conf. [christophe Becavin]
- Add read the docs conf. [christophe Becavin]
- Add read the docs conf. [christophe Becavin]
- Add read the docs conf. [christophe Becavin]
- Add all examples. [christophe Becavin]
- Add all examples. [christophe Becavin]
- Update README.md. [drbecavin]
- Update README.md. [drbecavin]
- Add example 2 and 3. [christophe Becavin]
- Add html examples. [christophe Becavin]
- Merge branch 'main' of github.com:becavin-lab/checkatlas. [christophe
  Becavin]
- Update README.md. [drbecavin]
- Add html examples. [christophe Becavin]
- Add Example 1 and 2. [christophe Becavin]
- Add examples to checkatlas. [christophe Becavin]
- Update README.md. [drbecavin]
- Reformat code. [christophe Becavin]
- #18 Add os.sep in filename, Put multihread=False by default.
  [christophe Becavin]
- Close #18, add folders.py to save all files in checkatlas_files/
  folder. [christophe Becavin]
- Merge branch 'main' of github.com:becavin-lab/checkatlas. [christophe
  Becavin]
- Update README.md. [drbecavin]

  Add summary
- Add pycharm run configuration. [christophe Becavin]
- First big update : Labmeeting presentation. Prototyp is ready.
  [christophe Becavin]
- Add Antoine metrics, add Dask. [christophe Becavin]
- Release: version  🚀 [christophe Becavin]
- Reformat and lint corrected. [christophe Becavin]
- Modifi atlas.py return statement. [christophe Becavin]
- Merge branch 'main' of github.com:becavin-lab/checkatlas. [christophe
  Becavin]
- Update README.md. [drbecavin]
- Add files for github workflow and template. [christophe Becavin]
- Clean project files. [christophe Becavin]
- Remove github templates from project creation. [christophe Becavin]
- Add git ignore. [christophe Becavin]
- Add git ignore. [christophe Becavin]
- ✅ Ready to clone and code. [drbecavin]
- First commit with many changes. [christophe Becavin]
- Initial commit. [drbecavin]


