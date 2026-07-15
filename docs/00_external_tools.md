# External tools

Parts of this pipeline run outside R. `scripts/01_main_analysis.R` *reads* the
outputs of these tools but does not *produce* them. This document records what
was run, so the chain from GenBank accession to biogeographic inference is
complete.

Fill in every bracketed field before publication. A reader who cannot reproduce
a step will assume the step is arbitrary.

---

## 1. Sequence alignment

| | |
|---|---|
| Tool | ClustalW |
| Version | 2.1 |
| Invoked from | R, inside `01_main_analysis.R` |
| Parameters | defaults |

Alignments were trimmed to 50% site occupancy.

Markers: four mitochondrial (`cytb`, `cox1`, `nd2`, `nd4`) and four nuclear
(`c-mos`, `jun`, `nt3`, `rag1`).

## 2. Phylogenetic inference

| | |
|---|---|
| Tool | BEAST |
| Version | [v2.6] |
| Citation | Bouckaert et al. (2019) |
| Substitution model | Selected per gene (see Table S1 in main document) |
| Clock model | lognormal relaxed |
| Tree prior | Coalescent Constant Population |
| Calibrations | see Table S1 in main document |
| Chain length | see Table S1 in main document |
| Convergence | ESS > 100 for all parameters, checked in Tracer v1.7.2 |
| Burn-in | 10% |

XML files: `docs/beast_xml/<genus>.xml`.
Output consumed by R: `<genus>_MCC.tre`.

## 3. Species delimitation

Three methods were applied. A tip was collapsed into an OTU when **at least two
of the three methods agreed** (implemented as connected components of an
agreement graph, in `01_main_analysis.R`).

### GMYC

| | |
|---|---|
| Tool | R package `splits` |
| Version | 1.0.20 |
| Model | single |
| Run inside | `01_main_analysis.R` |

### mPTP

| | |
|---|---|
| Tool | `mptp` |
| Version | 0.2.5 |
| Command | `mptp --mcmc 50000000 --multi --mcmc_sample 1000000 --mcmc_burnin 1000000 --tree_file <genus>_MCC.newick --output_file <genus>_result_mptp.txt --outgroup <outgroup>` |
| Output consumed by R | `<genus>_result_mptp.txt` |

### ASAP

| | |
|---|---|
| Tool | `SpartExplorer` |
| Source | https://spartexplorer.mnhn.fr/delimitation |
| Input | `<Genus>_merged_trimmed_sequences_0.50_NEW.fas` |
| Output consumed by R | `<Genus>_result_asap.txt` |

## 4. Geographic areas

Areas A, B and C were delimited in QGIS v3.40.14. See
[`../data/geography/README.md`](../data/geography/README.md).

---

