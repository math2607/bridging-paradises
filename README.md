# Bridging Paradises: asymmetric south-to-north dispersals in Squamate support a complex emergence of the Isthmus of Panama

[![Code: MIT](https://img.shields.io/badge/code-MIT-blue.svg)](LICENSE)
[![Data: CC BY 4.0](https://img.shields.io/badge/data-CC%20BY%204.0-lightgrey.svg)](LICENSE-DATA)
[![DOI](https://zenodo.org/badge/1301905151.svg)](https://doi.org/10.5281/zenodo.21383597)

Code and data supporting a manuscript in preparation on squamate dispersal
across the Isthmus of Panama. The full citation will be added here once the
article is published.

Archived snapshot: https://doi.org/10.5281/zenodo.21383597

---

## Overview

We reconstructed the biogeographic history of eight Neotropical squamate genera
(*Anolis*, *Bothriechis*, *Bothrops*, *Dipsas*, *Marisora*, *Micrurus*,
*Oxybelis*, *Sibon*) to test when and in which direction lineages crossed the
Isthmus of Panama.

For each genus we assembled multilocus alignments from GenBank (four
mitochondrial and four nuclear markers), delimited species with three
independent methods (GMYC, mPTP, ASAP), estimated time-calibrated trees in
BEAST, and fitted six biogeographic models in `BioGeoBEARS` across three areas.
Ancestral ranges were summarised with 1,000 Biogeographical Stochastic Maps
(BSM) per genus, from which we extracted directional dispersal counts and their
timing.

**Area B is defined biologically, not politically.** It is the union of two
biogeographic provinces of Morrone (2014) — Guatuso-Talamanca and
Puntarenas-Chiriquí — so an inferred crossing reflects movement across a
biogeographic boundary rather than a national border. Areas A and C are the
emergent land north and south of B. Full provenance in
[`data/geography/README.md`](data/geography/README.md).

## Repository structure

```
.
├── README.md
├── LICENSE                    MIT — everything under scripts/
├── LICENSE-DATA               CC BY 4.0 — everything under data/ and results/
├── CITATION.cff
├── .gitignore
│
├── scripts/
│   ├── 00_download_gbif.R     Retrieves the occurrence dataset from its DOI
│   ├── 00_functions.R         Helper functions sourced by the main script
│   └── 01_main_analysis.R     Full pipeline
│
├── data/
│   ├── README.md              PROVENANCE AND LICENSING OF EVERY FILE
│   ├── raw/				   WHERE THE RAW GBF FILE WILL BE DOWNLOADED
│   ├── geography/
│   │   ├── README.md          How A/B/C were built, and under what licence
│   │   └── A.shp, B.shp, C.shp   (+ .dbf, .shx, .prj)
│   └── biogeobears/
│       └── <Genus>/
│           ├── <Genus>_tree_final.newick
│           ├── <Genus>_MCC.TRE
|           └── <Genus>_geog_final.txt
│
├── results/
│   ├── supplementary/														For saving table for developing table s1
│   └── <Genus>/
│       ├── <Genus>_teststable_ALL_MPN.txt                                  LRT between nested models
│       ├── <Genus>_restable_ALL_MPN.txt                                    Model comparison (LnL, AIC, AICc)
│       ├── <Genus>_restable_AIC_rellike_ALL_MPN.txt                        AIC
│       ├── <Genus>_restable_AIC_rellike_ALL_formatted_MPN.txt              AIC (formatted)
│       ├── <Genus>_restable_AICc_rellike_ALL_MPN.txt                       AICc
│       ├── <Genus>_restable_AICc_rellike_ALL_formatted_MPN.txt             AICc (formatted)
│       └── BSM/
│           ├── <Genus>_ana_ALL_MASTER_table.txt                            Table with all anagenetic processes
|           ├── <Genus>_clado_ALL_MASTER_table.txt                          Table with all cladogenetic processes
|           ├── <Genus>_<model>_1000BSMs_v1.pdf                             PDF with all 1000 BSM
|           ├── <Genus>_<model>_single_stochastic_map_n1                    PDF with single (first) BSM
|           └── out/
|               ├── <Genus>_<model>_<process>_fromto_means.txt              Count means for specific processes (anagenetic, cladogenetic, a, d, foudnder) between area pairs
|               ├── <Genus>_<model>_<process>_fromto_sds.txt                Count stds for specific processes (anagenetic, cladogenetic, a, d, foudnder) between area pairs
|               ├── <Genus>_<model>_summary_counts_BSMs.txt              	Summary of all processes between all 1000 BSMs
|               ├── <Genus>_<model>_unique_<process>_counts.txt             Unique counts for specific cladogenetic processes
|               ├── <Genus>_<model>_histograms_of_event_counts.pdf          Summary graph of processes counts
|               └── <Genus>_<model>_ML_vs_BSM.pdf                           State probability as mean of BSM vs under ML model
│
├── figures/
│   ├── ALL_graph_results.pdf		                                        Age of dispersals joinning all genera
│   ├── <Genus>_graph_results.pdf		                                    Age of dispersals for each genus
│   └── <Genus>_<model>_MPN.pdf												BioGeoBEARS results for each model and 
│
└── docs/
    ├── 00_external_tools.md   Software run outside R, with versions
    └── session_info.txt       sessionInfo() from the published run
```

## Reproducing the analysis

**Requirements.** R ≥ 4.3 and the external tools listed in
[`docs/00_external_tools.md`](docs/00_external_tools.md). 

1. Clone and open `bridging-paradises` in RStudio. All paths are relative
   to the project scripts/ folder.

2. Retrieve the third-party data we cannot redistribute (see
   [`data/README.md`](data/README.md) for why):
   ```r
   source("scripts/00_download_gbif.R")   # 479 MB zipped -> ~2.3 GB
   ```

3. Run the external phylogenetic and delimitation steps described in
   [`docs/00_external_tools.md`](docs/00_external_tools.md). 

4. Run the analysis:
   ```r
   source("scripts/01_main_analysis.R")
   ```

**Determinism.** `GenSA` optimisation and stochastic mapping are seeded at the
top of `01_main_analysis.R`. Results reproduce on the same platform. Small
numerical differences across platforms and `BioGeoBEARS` versions are expected
and do not change the qualitative conclusions.

## Data sources and licensing

Every file in `data/` is documented in [`data/README.md`](data/README.md), with
its origin, its licence, and whether we redistribute it or provide a download
script instead.

**Not redistributed here:**

- **GBIF occurrence records.** The download aggregates 3,780,267 records 
  under heterogeneous licences (CC0, CC BY, CC BY-NC), so the
  aggregate has no single licence under which it could be shared. We cite the
  download DOI — GBIF.org (20 February 2026),
  https://doi.org/10.15468/dl.n295k4 — and provide
  `scripts/00_download_gbif.R`.
- **GADM administrative units**, downloaded at run time by `geodata::gadm()`
  into a temporary directory, under GADM's own terms.

**Redistributed here:** the area polygons A, B and C, under CC BY 4.0. They
derive from Natural Earth (public domain) and from Löwenberg-Neto (2014,
CC BY 3.0). If you reuse them, cite Löwenberg-Neto (2014) and Morrone (2014).

Code in `scripts/` is MIT. Data and results generated by us are CC BY 4.0.

## Citation

Please cite both the article and the archived repository:

NOTE: the citation will be provided once the paper is published

```bibtex
@article{[KEY],
  author  = {[AUTHORS]},
  title   = {[TITLE]},
  journal = {[JOURNAL]},
  year    = {[YEAR]},
  doi     = {[ARTICLE DOI]}
}

@software{[KEY]_code,
  author    = {[AUTHORS]},
  title     = {Code and data for: [TITLE]},
  year      = {[YEAR]},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.XXXXXXX}
}
```

## Contact

Matheus Pontes Nogueira, https://orcid.org/0000-0002-2146-1346, mpnogueira26@gmail.com — Programa de Pós-Graduação em Evolução e
Diversidade, Universidade Federal do ABC, Brazil.

Issues and questions are welcome through the repository issue tracker.

## Acknowledgements

This research was funded by Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES 88887.645493/2021-00), Fundação de Amparo à Pesquisa do Estado de São Paulo (FAPESP 2020/12658-4, 2022/05543-1, and 2023/06676-8), and Conselho Nacional de Desenvolvimento Científico e Tecnológico (CNPq 307956/2022-9). This study was financed in part by the Coordenação de Aperfeiçoamento de Pessoal de Nível Superior - Brasil (CAPES) - Finance Code 001. We are deeply grateful to Laura R.V. Alencar for her invaluable guidance, inspiring discussions, and significant contributions throughout the development of this study. We thank Felipe Grazziotin for helping with R scripts for plotting the results. We also thank the Laboratório de Evolução e Diversidade 1 and the Graduate Program in Evolution and Diversity of Universidade Federal do ABC for the support. No permits were required in the realization of this work.

