# Site facts verified against DLNR (2013), and the fence reconciliation

Read this before writing or checking any sentence of the study-system
description, before citing DLNR (2013), and before touching the ungulate
exclusion dates. Page numbers are pages of the PDF at
`references/Puu-Makaala-LRMP-Managment-Plan.pdf` (also at
`provided/NEON Docs/PlotCharacterization/Docs/`), read with `pdftotext -layout`
on 2026-10-19.

## What the management plan actually says

| Fact | Value | Page |
|---|---|---|
| Reserve area | 18,706 ac = **7,570 ha** (4,899 ha original + 2,671 ha of the former Kulani property added Nov 2010 by EO 4338) | 4, 6, 21 |
| Established | 1981, Governor's Executive Order 3102 | 4, 5 |
| Elevation | 2,800-6,229 ft = 853-1,899 m | 8 |
| Rainfall | 79-157 in = 2,000-4,000 mm | 8 |
| Oldest flows | older than 10,000 yr | 8 |
| Kulani Correctional Facility | ~280 ac = 113 ha, interior to the boundary, used by DOD-YCA under revocable permit | 6, 21 |
| Prison camp history | Kulani portion used as a prison camp since 1946; certain areas used for logging and ranching | 28 |
| Cattle | removed in 2005 from the **adjoining former pasture** (342 ac) proposed for addition to the reserve, not from tracts already inside it | 21 |
| Fence units | nine ungulate-proof units, ~9,600 ac = 3,885 ha, built 1991-2003 in the Kulani and Pu'u Maka'ala sections | 35 |
| Units named | Lava Flow, Na Lua Mahoe, Kulani Cone, Pu'u Kipu, Wright Road, 'Aku, Mauna Loa Boys School, South Boundary, North Boundary | 35 |
| Proposed units | Upper 'Ola'a (1,000 ac) and Army Road (1,000 ac) listed as **proposed**, not built | 44, 55 |

**The manuscript carried 7,580 ha through v0_4; the correct figure is 7,570 ha.**
Corrected in v0_5.

## The fence date reconciliation

The apparent conflict between the 1991-2003 exclusion window and the enclosure
completion years in the state layer is **not a conflict**. DLNR (2013) gives no
per-unit completion dates at all. The layer
(`data/PuuMakaala_Fence_Units_20240517.zip`) holds 13 units:

| Unit | ha | Completed | UngFree |
|---|---|---|---|
| 500 | 372.7 | 1992 | Yes |
| Mauna Loa Boys School | 362.2 | 1995 | Yes |
| Puu Kipu | 1064.5 | 1997 | Yes |
| Aku | 83.9 | 1998 | Yes |
| Kulani Cone | 596.8 | 1998 | Yes |
| North Boundary | 294.8 | 2002 | Yes |
| Wright Road | 970.2 | 2003 | Yes |
| South Boundary | 422.2 | 2004 | Yes |
| Stainback | 421.1 | 2014 | Yes |
| Upper Olaa | 756.4 | 2016 | **No** |
| Army Road | 412.0 | 2019 | **No** |
| Kata-Stein | 600.9 | *none* | Yes |
| Unencumbered Kulani Pasture | 116.1 | *none* | Yes |

Seven of the nine units the plan names appear here with completion years
1995-2004. Stainback, Upper 'Ola'a and Army Road postdate the plan, and the
plan lists the latter two as proposed. Six enclosures contain NEON plots
(Army Road, Kata-Stein, Mauna Loa Boys School, North Boundary, Puu Kipu,
Stainback); three plots fall outside any enclosure.

## Citation corrections

- The cultural study cited in drafts through v0_4 as `Cultural-Study-Puu-Makala`
  is **Maly, K. and O. Maly (2004)**, *He Mo'olelo 'Aina: A Cultural Study of
  the Pu'u Maka'ala Natural Area Reserve, Districts of Hilo and Puna, Island of
  Hawai'i*, Kumu Pono Associates LLC, 31 May 2004. Copy in `references/`.
- Lava flow ages are **Sherrod et al. (2021)**, USGS SIM 3143, already in
  `MS2/references.docx`.
