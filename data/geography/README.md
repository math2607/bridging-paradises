# Biogeographic areas A, B and C

The three areas used in every `BioGeoBEARS` analysis. They are a derivative
work; this file records how they were made and under what terms they are
redistributed.

## Definition

**Area B** is the union of two biogeographic provinces of Morrone's (2014)
regionalisation of the Neotropical region:

| Province | `NUMB` | Dominion | Subregion |
|---|---|---|---|
| Guatuso-Talamanca | 18 | Pacific | Brazilian |
| Puntarenas-Chiriquí | 19 | Pacific | Brazilian |

**Area A** is the emergent land north of B.
**Area C** is the emergent land south of B.

A, B and C are **not** administrative units. B is a biogeographic unit defined
by endemicity, and it does not coincide with the political Isthmus of Panama.
This is deliberate: the boundary is biological rather than national.

## How the polygons were built

In QGIS:

1. Load `ne_50m_land` (Natural Earth, land polygons including major islands).
2. Load `Lowenberg_Neto_2014.shp` from Löwenberg-Neto, P. (2014) and select provinces 18 and 19 (`NUMB`
   field). Their union is **B**.
3. Clip `ne_50m_land` by B. Land north of B becomes **A**; land south of B
   becomes **C**.
4. Export A, B and C as separate shapefiles, EPSG:4326.

## Sources and licences

| Layer | Source | Licence |
|---|---|---|
| Land polygons | Natural Earth, `ne_50m_land` v4.0.0, 1:50m physical vectors | Public domain |
| Biogeographic provinces | Löwenberg-Neto, P. (2014) Neotropical region: a shapefile of Morrone's (2014) biogeographical regionalisation. *Zootaxa* 3802(2): 300. https://doi.org/10.11646/zootaxa.3802.2.12 | CC BY 3.0 |
| Underlying regionalisation | Morrone, J.J. (2014) Biogeographical regionalisation of the Neotropical region. *Zootaxa* 3782(1): 1–110. https://doi.org/10.11646/zootaxa.3782.1.1 | — |

The Löwenberg-Neto shapefile's own FGDC metadata states access constraints:
*none*; use constraints: *authorship credit*. The Zootaxa correspondence is
published under CC BY 3.0.

`A.shp`, `B.shp` and `C.shp` are therefore redistributed under **CC BY 4.0**, a
compatible licence, with attribution. CC BY 3.0 is not copyleft, so a derivative
may carry a different attribution licence.

**If you reuse these polygons, cite Löwenberg-Neto (2014) and Morrone (2014),
not only this repository.**
