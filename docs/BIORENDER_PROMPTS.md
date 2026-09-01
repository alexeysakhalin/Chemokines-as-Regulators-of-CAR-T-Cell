# BioRender specifications for revised and planned figures

## Scope and use

The prompts below cover the five revised conceptual figures, the revised graphical abstract, and the planned data-dependent Figure 6. They are written for direct use in BioRender. All labels are in English.

Figures 1-5 and the graphical abstract are evidence-constrained conceptual schematics. Figure 6 is a hybrid experimental figure: BioRender may be used for its non-quantitative study-design strip and final assembly, but every graph, point, confidence interval, spatial coordinate, sample count, and statistical annotation must come from the reproducible Python workflow and verified source data.

The demonstration outputs in `results/` are synthetic software tests. They must not be inserted into the manuscript, used to infer a biological direction, or used as visual templates for invented values.

**Numbering note:** the submitted PDF declares five numbered figures and provides legends only for Figures 1-5, but the graphical abstract is attached under the filename `Figure 6.JPEG`. Before adding the planned analytical Figure 6, rename that attachment to `Graphical_Abstract` or the journal's required graphical-abstract filename. Otherwise the new figure and the graphical abstract will conflict.

## Shared visual standard

- Use a white background, flat vector graphics, no 3D rendering, no photorealism, no decorative textures, and no unnecessary shadows.
- Use Arial or Helvetica throughout. Use sentence case. Keep labels at least 8 pt at final journal size and panel letters 12-14 pt bold.
- Use short direct labels. Define abbreviations in the figure legend, not inside crowded panels.
- Use the same icon and color for a biological entity in every figure.
- CAR-T cell: navy `#1F4E79`; endogenous CD8 T cell: teal `#1B8A8F`; NK cell: cyan `#27A9B8`; cDC1 or other dendritic cell: green `#4C9F70`; tumor cell: salmon `#E8877D`; CAF or CXCL12-rich stroma: purple `#8E6BBE`; suppressive myeloid cell: orange `#D98C3F`; blood vessel: muted red `#B85C5C`; lymphatic vessel: pale green `#8DBFA8`; extracellular matrix: gray-blue `#9AA8B5`.
- Use solid arrowheads for a supported direction of movement or signaling, dotted arrows for a proposed or context-dependent route, blunt-ended lines for inhibition, and dashed outlines for conceptual groupings. Never use arrow thickness as a quantitative encoding.
- Use direct labels and distinct shapes in addition to color. The figure must remain interpretable in grayscale and for readers with color-vision deficiency.
- Do not add percentages, fold changes, P values, survival curves, heatmaps, UMAP clusters, spatial densities, or relative arrow widths unless they are supplied from verified data.
- Do not convert an association into causality, an RNA signal into a protein gradient, a total protein measurement into a proteoform-specific result, or a preclinical result into a clinical claim.

---

## Revised Figure 1 - Chemokine-guided entry, positioning, retention, and egress

**Figure type:** conceptual mechanism; complete replacement of the current Figure 1.

**Caption intent:** contrast a chemokine-permissive and an immune-excluded tumor while distinguishing entry, intratumoral positioning, retention, and egress; show that cellular source, sequence, concentration, and tissue location determine the effect of a chemokine axis.

**Evidence anchors:** Böttcher et al., *Cell* (2018), [doi:10.1016/j.cell.2018.01.004](https://doi.org/10.1016/j.cell.2018.01.004); Spranger et al., *Cancer Cell* (2017), [doi:10.1016/j.ccell.2017.04.003](https://doi.org/10.1016/j.ccell.2017.04.003); Steele et al., *Nature Immunology* (2023), [doi:10.1038/s41590-023-01443-y](https://doi.org/10.1038/s41590-023-01443-y); Donnadieu et al., *Journal of Leukocyte Biology* (2020), [doi:10.1002/JLB.1MR0520-746R](https://doi.org/10.1002/JLB.1MR0520-746R); Foeng et al., *Cell Reports Medicine* (2022), [doi:10.1016/j.xcrm.2022.100543](https://doi.org/10.1016/j.xcrm.2022.100543).

### Prompt

```text
Create a publication-quality landscape scientific schematic for a peer-reviewed immunology review. Use a clean flat-vector BioRender style, white background, no photorealism, no 3D effects, and an aspect ratio of approximately 1.65:1. Use English labels only and make all text legible at two-column journal width.

TITLE: Chemokine-guided CAR-T-cell entry, positioning, retention, and egress

LAYOUT: Draw two matched tumor cross-sections separated by a thin gray divider. Label the left panel A, “Inflamed / chemokine-permissive tumor”, and the right panel B, “Immune-excluded tumor”. In both panels place a blood vessel at the top, tumor-cell nests in the lower center, stromal tissue around the nests, and a peripheral lymphatic vessel at the outer edge. Use the same biological icons in both panels.

PANEL A - INFLAMED / CHEMOKINE-PERMISSIVE TUMOR:
1. Show activated tumor endothelium with ICAM1 and VCAM1 and local presentation of CXCL9, CXCL10, and CXCL11.
2. Show CXCR3-positive CAR-T cells leaving the vessel and following a directional CXCL9/CXCL10/CXCL11 field toward antigen-positive tumor nests. Label “CXCR3-dependent entry and positioning”.
3. Add a compact sequential relay inside the tumor: an NK cell releases CCL5 and XCL1; arrows recruit a cDC1; the cDC1 releases CXCL9 and CXCL10; arrows recruit and position CXCR3-positive cytotoxic T cells. Label “NK cell -> cDC1 -> cytotoxic T-cell relay”. Do not draw a direct NK-cell-to-T-cell bypass arrow.
4. At the tumor nest show CAR-antigen recognition followed by local IFN-gamma and TNF release. Add a small feedback arrow from local inflammatory signaling to CXCL9/CXCL10 production.
5. Show one tissue-resident-memory-like T cell with low S1PR1. Label “S1PR1 downregulation: general tissue-retention program”. Keep this separate from CXCL12-dependent lymphatic egress.

PANEL B - IMMUNE-EXCLUDED TUMOR:
1. Show abnormal, poorly activated endothelium with reduced ICAM1/VCAM1 and limited T-cell extravasation.
2. Surround tumor nests with dense collagen-rich extracellular matrix and purple CAFs. Show CXCR4-positive CAR-T cells accumulating in stroma or moving along collagen fibers without reaching tumor cells.
3. Show high peritumoral CAF-derived CXCL12. Use curved or reversed migration arrows labeled “stromal sequestration / fugetaxis at high local CXCL12”. Do not depict CXCL12 as universally repellent.
4. At the lymphatic vessel show CXCL12 positioning CXCR4-positive endogenous CD8 T cells near the tumor edge and a dotted arrow into the lymphatic vessel labeled “context-dependent lymphatic egress”. Add a small inset: “weak or absent antigen encounter: CXCR4-high, egress permissive” versus “high-affinity antigen encounter: CXCR4 down, ACKR3 up, egress reduced”. Mark the inset “endogenous CD8 T-cell evidence; not CAR-T-specific”.
5. Add regulatory T cells and suppressive myeloid cells in stroma, but keep them secondary to the spatial trafficking message.

BOTTOM STRIP: Add a four-step axis: “Entry -> intratumoral positioning -> retention -> egress”. Below it write: “Outcome depends on ligand source, concentration, spatial presentation, receptor state, and timing.”

COLOR AND LABEL RULES: CAR-T cells navy; endogenous CD8 T cells teal; NK cells cyan; cDC1 green; tumor cells salmon; CAFs purple; suppressive myeloid cells orange; blood vessels muted red; lymphatics pale green; CXCL9/10/11 field blue; CXCL12 field violet. Use direct labels and a small legend.

EXCLUSIONS: Do not show all axes as simultaneously active in every tumor. Do not imply that CXCR4 is uniformly protumor or antitumor. Do not label the Steele mechanism as CAR-T-specific. Do not depict S1PR1 as a chemokine receptor. Do not add cell counts, percentages, fold changes, or quantitative arrow widths.
```

---

## Revised Figure 2 - CAR-T activation, tumor-cell death competence, and immune amplification

**Figure type:** conceptual mechanism; complete replacement of the current Figure 2.

**Caption intent:** show that chemokine induction, tumor-cell death, DAMP release, and recruitment are parallel but conditionally coupled processes; distinguish tumor-cell death pathways from CAR-T-cell loss and state explicitly that the integrated circuit is proposed rather than demonstrated end to end in one CAR-T solid-tumor model.

**Evidence anchors:** Galluzzi et al., *Nature Reviews Immunology* (2017), [doi:10.1038/nri.2016.107](https://doi.org/10.1038/nri.2016.107); Workenhe et al., *OncoImmunology* (2021), [doi:10.1080/2162402X.2021.1893466](https://doi.org/10.1080/2162402X.2021.1893466); Karki et al., *Cell* (2021), [doi:10.1016/j.cell.2020.11.025](https://doi.org/10.1016/j.cell.2020.11.025); Liu et al., *Science Immunology* (2020), [doi:10.1126/sciimmunol.aax7969](https://doi.org/10.1126/sciimmunol.aax7969).

### Prompt

```text
Create a publication-quality horizontal mechanistic figure for a peer-reviewed immunology review. Use a clean flat-vector BioRender style, white background, minimal text, consistent arrow grammar, English labels only, and an aspect ratio of approximately 1.75:1. Use panel letters A-D.

TITLE: CAR-T-cell activation, tumor-cell death competence, and local immune amplification

PANEL A - CAR IMMUNE SYNAPSE:
Show a navy CAR-T cell engaging an antigen-positive salmon tumor cell. At the synapse show perforin and granzyme B delivery. Outside the synapse show IFN-gamma and TNF acting on the target and neighboring tumor or stromal cells. Label the routes “contact-dependent cytotoxicity” and “paracrine cytokine exposure”. Add a small separate FAS-FASL branch directed toward the CAR-T cell and label “context-dependent effector-cell loss”; do not merge it with tumor-cell death.

PANEL B - TUMOR-CELL INFLAMMATORY SIGNALING:
Inside the tumor cell draw two parallel receptor pathways. First: IFN-gamma -> IFN-gamma receptor -> JAK -> STAT1 -> IRF1, leading to “antigen processing and presentation” and “CXCL9 / CXCL10 / CXCL11”. Show these ligands recruiting CXCR3-positive CAR-T cells, endogenous CD8 T cells, and NK cells. Second: TNF -> TNFR1, branching to “NF-kappaB inflammatory / survival program” and “FADD-caspase-8 death-complex formation”. Under the survival branch show c-FLIP, BCL-XL, and cIAPs as representative factors.

PANEL C - CONDITIONAL DEATH-PATHWAY COMPETENCE:
Create three separate outcome boxes, not one universal combined pathway.
1. “Extrinsic apoptosis”: FADD and active caspase-8.
2. “Necroptosis”: RIPK1 -> RIPK3 -> MLKL, with the note “requires functional RIPK3 and MLKL”.
3. “Pyroptosis”: granzyme B -> caspase-3 -> GSDME cleavage, with the note “reported in selected CAR-T target systems; requires expressed, cleavable GSDME”.
Place a dotted bracket around the three boxes labeled “integrated PANoptotic signaling can occur in defined contexts”. Do not show PANoptosis as the obligatory outcome.

PANEL D - TWO TISSUE-LEVEL CONSEQUENCES:
Upper route: “Productive tumor-cell death”. Show release or exposure of tumor antigens, ATP, HMGB1, and nucleic acids; a dendritic cell acquires tumor material and cross-presents antigen; endogenous CD8 T cells are recruited or activated. Label “DAMP-dependent activation + chemokine-directed recruitment”.
Lower route: “Death-resistant inflammatory tumor”. Show a surviving tumor cell continuing inflammatory and chemokine output, repeated CAR-T contacts, suppressive myeloid-cell accumulation, and progressive CAR-T dysfunction. Label “chronic antigen and inflammatory signaling without effective clearance”.

BOTTOM NOTE: “The individual steps are supported, but the full CAR-T-specific circuit remains a proposed integrated model.”

COLOR AND LABEL RULES: CAR-T cells navy; endogenous CD8 T cells teal; NK cells cyan; dendritic cells green; tumor cells salmon; DAMPs gold; chemokine arrows blue; death modules red/magenta; NF-kappaB survival module amber; suppressive myeloid cells orange.

EXCLUSIONS: Do not label apoptosis, necroptosis, pyroptosis, or DAMP release as uniformly immunogenic. Do not state that chemokine induction proves tumor killing. Do not show necroptosis without RIPK3 and MLKL. Do not add pathway activation values, response rates, or quantitative arrow widths.
```

---

## Revised Figure 3 - Suppressive chemokine circuits, proteoforms, and spatial context

**Figure type:** conceptual mechanism; expanded replacement of the current Figure 3.

**Caption intent:** explain that suppressive recruitment and T-cell positioning depend on ligand source and tissue context, and that chemically or proteolytically distinct chemokine forms can have different functions.

**Evidence anchors:** Molon et al., *Journal of Experimental Medicine* (2011), [doi:10.1084/jem.20101956](https://doi.org/10.1084/jem.20101956); Abel et al., *Journal of Immunology* (2004), [doi:10.4049/jimmunol.172.10.6362](https://doi.org/10.4049/jimmunol.172.10.6362); Lesch et al., *Nature Biomedical Engineering* (2021), [doi:10.1038/s41551-021-00737-6](https://doi.org/10.1038/s41551-021-00737-6); Foeng et al., *Cell Reports Medicine* (2022), [doi:10.1016/j.xcrm.2022.100543](https://doi.org/10.1016/j.xcrm.2022.100543).

### Prompt

```text
Create a publication-quality landscape scientific schematic for an immunology review. Use a clean flat-vector BioRender style, a white background, English labels only, and a 2 x 2 modular layout with an aspect ratio of approximately 1.70:1. Use panel letters A-D and one shared conclusion strip.

TITLE: Suppressive chemokine circuits, proteoforms, and spatial context in solid tumors

PANEL A - SUPPRESSIVE-CELL RECRUITMENT:
Show three parallel ligand-receptor axes entering a tumor niche: CCL2 -> CCR2-positive monocytes followed by TAM or monocytic MDSC accumulation; CXCL8 -> CXCR1/CXCR2-positive neutrophils and polymorphonuclear MDSCs; CCL22 -> CCR4-positive regulatory T cells. Show the shared outcomes “suppressive-cell accumulation”, “reduced CAR-T access and function”, and “chronic inflammatory resistance”. Add the qualifier “context dependent”.

PANEL B - CCL2 PROTEOFORM EFFECT:
Place native CCL2 and nitrated CCL2 side by side. Native CCL2 forms a conventional CCR2-directed field. Nitrated CCL2 has a small chemical-modification marker and a blocked or weakened migration arrow for an antigen-specific endogenous T cell, while a monocyte/myeloid migration arrow remains. Label “reduced antigen-specific T-cell chemotaxis” and “myeloid recruitment relatively preserved”. Add “mechanistic T-cell evidence; not CAR-T-specific”. Do not portray nitration as complete loss of function for every leukocyte.

PANEL C - MEMBRANE AND SOLUBLE CXCL16:
Show a vascular, stromal, or myeloid source cell expressing transmembrane CXCL16. A CXCR6-positive T cell binds the membrane form; label “adhesion / local retention”. Show ADAM10 cleaving the extracellular domain and generating soluble CXCL16. Show a soluble CXCL16 field attracting a CXCR6-positive activated T cell; label “soluble chemoattractant”. Add a spatial warning: “vascular or myeloid localization may retain cells away from tumor nests”. Do not imply that one cited study experimentally isolated every membrane-versus-soluble effect.

PANEL D - CONTEXT-DEPENDENT CXCL12-CXCR4:
Show purple CAFs surrounding a tumor nest and producing CXCL12. Depict two different spatial configurations: an organized directional field that can recruit a CXCR4-positive effector T cell, and a high peritumoral concentration associated with stromal sequestration or fugetactic movement away from the tumor nest. Add a peripheral lymphatic vessel with a dotted egress arrow for a CXCR4-positive endogenous CD8 T cell. Label “recruitment, retention, fugetaxis, or egress depend on spatial organization”.

BOTTOM STRIP: “Chemokine abundance alone is insufficient: measure cellular source, chemical or proteolytic processing, matrix binding, concentration, and spatial localization.” Add small icons for spatial protein imaging, soluble-protein measurement, and functional migration testing.

COLOR AND LABEL RULES: effector T cells navy/teal; regulatory T cells violet; monocytes and TAMs brown-orange; MDSCs orange; neutrophils pale orange; CAFs purple; tumor cells salmon; CCL2 green; CXCL8 orange; CCL22 violet; CXCL16 cyan; CXCL12 purple.

EXCLUSIONS: Do not call total chemokine mRNA a functional gradient. Do not show membrane and soluble CXCL16 as interchangeable or directly comparable on one quantitative scale. Do not depict CXCL12-CXCR4 as uniformly suppressive. Do not add invented concentrations, migration percentages, or proteoform ratios.
```

---

## Revised Figure 4 - Product heterogeneity, persistence, exhaustion, and senescence-like dysfunction

**Figure type:** conceptual mechanism; complete replacement of the current Figure 4.

**Caption intent:** link infusion-product heterogeneity and repeated tumor interaction to distinct functional trajectories while separating exhaustion-associated and senescence-like programs.

**Evidence anchors:** Wang et al., *OncoImmunology* (2021), [doi:10.1080/2162402X.2020.1866287](https://doi.org/10.1080/2162402X.2020.1866287); Sarén et al., *Clinical Cancer Research* (2023), [doi:10.1158/1078-0432.CCR-23-0178](https://doi.org/10.1158/1078-0432.CCR-23-0178); Bai et al., *Science Advances* (2022), [doi:10.1126/sciadv.abj2820](https://doi.org/10.1126/sciadv.abj2820); Bai et al., *Nature* (2024), [doi:10.1038/s41586-024-07762-w](https://doi.org/10.1038/s41586-024-07762-w); Herzberg et al., *Journal for ImmunoTherapy of Cancer* (2025), [doi:10.1136/jitc-2024-010709](https://doi.org/10.1136/jitc-2024-010709); Liu et al., *Nature Communications* (2018), [doi:10.1038/s41467-017-02689-5](https://doi.org/10.1038/s41467-017-02689-5).

### Prompt

```text
Create a publication-quality landscape mechanistic figure for a peer-reviewed CAR-T and tumor-immunology review. Use a clean flat-vector BioRender style, white background, English labels only, and an aspect ratio of approximately 1.70:1. Use a top-to-bottom flow with panel letters A-C. Do not add quantitative data.

TITLE: CAR-T product state, persistence, and distinct routes to dysfunction

PANEL A - HETEROGENEOUS INFUSION PRODUCT:
Show a mixed CAR-positive T-cell population with conceptual states labeled “early-memory-like”, “effector / polyfunctional”, “activated”, and “exhaustion-associated”. Add a visually separate “senescence-like / low proliferative reserve” state. Use equal-sized groups without percentages. Add two upstream influences: “prior treatment and starting-cell quality” and “manufacturing duration / repeated ex vivo stimulation”.

PANEL B - TWO TUMOR-INTERACTION TRAJECTORIES:
Left route, “productive clearance and persistence”: show an early-memory-like or polyfunctional CAR-T cell entering a tumor, killing an antigen-positive tumor cell, reducing antigen burden, and leaving a small pool of persistent memory-like cells. Add restrained local chemokine recruitment and DAMP-assisted dendritic-cell activation. Label “tumor clearance with preserved proliferative reserve”.
Right route, “chronic failure”: show a death-resistant tumor cell that remains responsive to TNF and IFN-gamma but is not efficiently eliminated. Include repeated CAR-T contacts, persistent antigen, chronic inflammatory chemokines, hypoxia, suppressive cytokines, TAMs, and MDSCs. Show progressive loss of cytotoxicity and cytokine production. Label “repeated stimulation without effective clearance”.

PANEL C - DISTINCT DYSFUNCTIONAL PROGRAMS:
Create two separate boxes with no merging arrow.
1. “Exhaustion-associated program”: chronic antigen signaling; TOX, NR4A family, and BATF; increased inhibitory-receptor programs; reduced cytokine production and killing.
2. “Senescence-like program”: proliferative history, metabolic competition, telomere or DNA-damage stress, and reduced proliferative reserve.
Between the boxes write “senescence is biologically distinct from exhaustion”. Show inflammatory and suppressive chemokine niches as contextual factors that can sustain repeated interactions, not as sole causes of either state.

SINGLE-CELL INSET: Add a small conceptual cluster map explicitly labeled “schematic, not quantitative”. Connect it to four readouts: “memory”, “polyfunctionality”, “cytotoxicity”, and “inhibitory / dysfunction programs”. Add “No single chemokine, including CXCL13, defines CAR-T-cell exhaustion by itself.”

COLOR AND LABEL RULES: early-memory-like cells navy; effector/polyfunctional cells bright blue; activated cells cyan; exhaustion-associated cells orange-red; senescence-like cells gray-purple; tumor cells salmon; suppressive niche orange/purple. Use direct labels rather than a color-only legend.

EXCLUSIONS: Do not use PD-1, LAG-3, TIM-3, CXCL13, or one transcription factor as a stand-alone definition of exhaustion. Do not merge exhaustion and senescence. Do not draw a UMAP with axes, cluster sizes, percentages, or trajectories that could be mistaken for observed data. Do not imply that chemokines alone cause dysfunction.
```

---

## Revised Figure 5 - Chemokine engineering, local delivery, and safety-controlled CAR-T design

**Figure type:** translational strategy schematic; complete replacement of the current Figure 5.

**Caption intent:** shift emphasis from a catalogue of chemokine receptors to chemokine engineering, source and timing, negative translational lessons, matched comparators, and explicit safety controls.

**Evidence anchors:** Lugassy et al., *PNAS* (2025), [doi:10.1073/pnas.2501791122](https://doi.org/10.1073/pnas.2501791122); Gerlza et al., *Protein Engineering, Design and Selection* (2019), [doi:10.1093/protein/gzz043](https://doi.org/10.1093/protein/gzz043); Adachi et al., *Nature Biotechnology* (2018), [doi:10.1038/nbt.4086](https://doi.org/10.1038/nbt.4086); Luo et al., *Clinical Cancer Research* (2020), [doi:10.1158/1078-0432.CCR-20-0777](https://doi.org/10.1158/1078-0432.CCR-20-0777); Eckert et al., *Molecular Therapy Oncolytics* (2020), [doi:10.1016/j.omto.2019.12.003](https://doi.org/10.1016/j.omto.2019.12.003); Moon et al., *OncoImmunology* (2018), [doi:10.1080/2162402X.2017.1395997](https://doi.org/10.1080/2162402X.2017.1395997).

### Prompt

```text
Create a publication-quality landscape strategy figure for a peer-reviewed immunology review. Use a clean flat-vector BioRender style, white background, English labels only, four vertical modules across the upper two-thirds, and one safety-and-monitoring layer across the bottom. Use an aspect ratio of approximately 1.80:1. Mark evidence level with small neutral badges: “preclinical”, “early clinical”, or “conceptual”.

TITLE: Chemokine engineering and controllable trafficking strategies for CAR-T cells

MODULE 1 - MATCHED RECEPTOR-LIGAND TRAFFICKING:
Show representative alternative products, not one universal cell. Include a CCR2b-engineered CAR-T cell following tumor-derived CCL2 and a separate CXCR1- or CXCR2-engineered CAR-T cell following a CXCL8-rich signal. Add a caution label: “benefit depends on ligand source and competing suppressive cells”. Keep this module compact.

MODULE 2 - ENGINEER THE CHEMOKINE:
Subpanel A, “protease resistance”: show DPP-4 removing the first two N-terminal residues of native CXCL9 or CXCL10, producing a receptor-binding but non-agonistic form. Next show N-terminal-glutamine CXCL9-Fc and CXCL10-Fc protected from DPP-4 cleavage, maintaining CXCR3 agonism. Label “DPP-4-resistant agonist”.
Subpanel B, “matrix binding”: show CXCL10 N20K with enhanced glycosaminoglycan binding on endothelium or extracellular matrix, supporting stable local presentation and T-cell migration in functional assays. Label “enhanced GAG binding”.
Mark both strategies “preclinical”.

MODULE 3 - LOCAL OR ACTIVATION-DEPENDENT PAYLOADS:
Show separate CAR-T designs expressing IL-7 plus CCL19, IL-7 plus CCL21, IL-7 plus CCL3, or CCR5 plus IL-12. Use separate cell icons and label “alternative constructs”. For IL-7/CCL19 show NFAT-dependent expression after CAR engagement rather than constitutive systemic secretion. State “CCL19 and CCL21 both signal through CCR7; comparative efficacy is model- and conditioning-dependent”. Add an oncolytic-virus icon delivering a local chemokine as an alternative intratumoral source.

MODULE 4 - SOURCE, TIMING, AND NEGATIVE LESSONS:
Show two evidence contrasts.
1. Oncolytic VSV-CXCL9 creates a tumor-to-blood chemokine field, followed by a neutral or blocked trafficking arrow labeled “did not measurably increase T-cell trafficking in the tested models”.
2. Local vaccinia-virus CXCL11 increases recruitment and efficacy, whereas constitutive CXCL11 secretion by modified T cells fails to recruit a subsequent T-cell dose and reduces T-cell function in the same study. Label “local viral delivery effective; constitutive T-cell delivery detrimental in this model”.
Present these as design cautions, not universal rules.

BOTTOM SAFETY AND MONITORING LAYER:
Show four independent control options: “activation-dependent promoter”, “reversible ON/OFF switch”, “inducible caspase-9 suicide switch”, and “hypoxia-responsive expression”. Add monitoring icons for organ-specific inflammation, cytokine release syndrome, neurotoxicity, biodistribution, and serial blood/tumor sampling. Write “Control systems can reduce exposure but do not abolish on-target off-tumor toxicity.” Add a small paired-design icon: “same donor; receptor-matched product versus otherwise identical CAR-T control”.

FINAL OUTCOME BAR: “Goal: improve tumor access while preserving function, spatial specificity, and controllability.” Add “No chemokine-engineered CAR-T product is currently approved; controlled comparative evidence is limited.”

COLOR AND LABEL RULES: CAR-T cells navy; engineered receptors cyan; engineered chemokines green; cytokine/chemokine payloads violet; oncolytic viruses orange; tumor cells salmon; safety systems navy/gray; supported beneficial arrows green; negative or caution findings amber/red.

EXCLUSIONS: Do not place CCR2b, CXCR1, CXCR2, CXCR3, and CCR5 on one universal product. Do not claim universal superiority of CCL21 over CCL19. Do not depict orthogonal IL-2/IL-2R-beta as a chemokine receptor system. Do not present conference abstracts as definitive efficacy. Do not imply clinical approval, increased survival, or improved response rates without verified clinical data.
```

---

## Planned Figure 6 - Integrated single-cell, spatial, protein, and functional validation

**Figure type:** hybrid experimental figure; data dependent. Do not add it to the manuscript as a results figure until source data, biological replication, quality control, and statistical outputs are available.

**Caption intent:** show whether a chemokine-modified product differs from an otherwise identical CAR-T control across patient- or donor-level state, spatial localization, CXCL16 form-specific protein measurements, trafficking function, and matched effect estimates without collapsing the evidence into an unvalidated composite score.

**Python outputs to use:** `figure6a_single_cell_states`, `figure6b_spatial_map`, `figure6c_cxcl16_protein_validation`, `figure6d_functional_assays`, and `figure6e_matched_product_effects`. Only outputs from an evidence-eligible manifest may enter a manuscript figure. Synthetic outputs must retain their warning and remain outside the manuscript.

**Evidence anchors for design:** Deng et al., *Nature Medicine* (2020), [doi:10.1038/s41591-020-1061-7](https://doi.org/10.1038/s41591-020-1061-7); Sheih et al., *Nature Communications* (2020), [doi:10.1038/s41467-019-13880-1](https://doi.org/10.1038/s41467-019-13880-1); Sarén et al., *Clinical Cancer Research* (2023), [doi:10.1158/1078-0432.CCR-23-0178](https://doi.org/10.1158/1078-0432.CCR-23-0178); Steffin et al., *Nature* (2025), [doi:10.1038/s41586-024-08261-8](https://doi.org/10.1038/s41586-024-08261-8).

### BioRender prompt for the non-quantitative study-design strip and assembly frame

```text
Create a publication-quality assembly frame for a hybrid experimental figure in a peer-reviewed CAR-T and tumor-immunology manuscript. Use a clean flat-vector BioRender style, white background, English labels only, and an aspect ratio of approximately 1.85:1. BioRender must create only the non-quantitative study-design strip, panel headings, and neutral borders. Leave all quantitative panel interiors blank for insertion of verified Python-generated vector plots. Do not generate data points, axes, P values, confidence intervals, sample counts, heatmaps, UMAPs, spatial coordinates, or effect sizes.

TITLE: Integrated validation of chemokine-guided CAR-T-cell state and trafficking

TOP STUDY-DESIGN STRIP:
Show one patient or independent donor as the biological unit. Split the same starting material into two otherwise matched products: “identical CAR-T control” and “chemokine-module CAR-T”. Show three coordinated sampling streams: “final product”, “serial blood”, and “pretreatment / on-treatment tumor”. Connect them to four orthogonal assays: “scRNA-seq + CITE-seq”, “spatial transcriptomics / multiplex imaging”, “soluble and membrane protein”, and “chemotaxis / retention / egress”. Add “technical wells, cells, fields, and regions are nested measurements”.

LOWER ASSEMBLY GRID: Create five clean labeled frames without plotted content.
A. “Patient-level CAR-T and endogenous T-cell states”. Reserve space for observed patient/donor points or paired lines across PRODUCT and protocol-defined timepoints. State labels: memory-like, effector, proliferating, interferon-responsive, exhaustion-associated, and senescence-associated only when prespecified and validated upstream.
B. “Spatial sources, CAR-T localization, and tissue landmarks”. Reserve space for a real tissue map plus a smaller patient-level distance or neighborhood summary. Landmarks: perfused vessel, tumor nest, stroma, lymphatic vessel, ligand source, and CAR-T cell.
C. “CXCL16 protein forms on assay-specific scales”. Reserve two separate subframes: membrane CXCL16 by non-permeabilized flow cytometry or quantitative imaging, and soluble CXCL16 by validated immunoassay or targeted mass spectrometry. Do not place both forms on a shared numerical axis unless a common calibrated unit is experimentally justified.
D. “Chemotaxis, retention, and egress”. Reserve separate subframes for these assays. Indicate expected controls with small labels: no ligand, forward gradient, uniform ligand, reverse gradient, receptor blockade, viability, and recovery.
E. “Matched product effect estimates”. Reserve a forest-plot frame with a vertical zero line, effect direction “favors identical control” to the left and “favors chemokine-module product” to the right, patient- or donor-level effect estimates, 95% confidence intervals, and endpoint-specific units.

BOTTOM INFERENCE STRIP: “Patient or independent donor is the unit of inference; effect size and uncertainty are primary; multiplicity is controlled within endpoint families.” Add “Orthogonal concordance is interpreted mechanistically, not as one composite score.”

COLOR AND LABEL RULES: identical control gray-blue; chemokine-module product navy; endogenous T cells teal; tumor cells salmon; ligand sources green; vessels muted red; stroma purple; uncertainty intervals dark gray. Use the same group colors in every inserted quantitative panel. Keep panel letters A-E large and aligned.

EXCLUSIONS: Do not draw placeholder curves or example dots. Do not fabricate sample counts, biological directions, response rates, significance stars, confidence intervals, patient identifiers, or spatial distributions. Do not use cells, spots, technical wells, image fields, or tissue regions as independent patients. Do not combine membrane CXCL16 MFI with soluble CXCL16 concentration on one inferential scale. Do not present synthetic software-test outputs as biological evidence.
```

### Quantitative panel rules for Python-generated content

- **Panel A:** show patient/donor-level observations and within-patient trajectories. Do not infer a cell state from one marker. CAR-positive and endogenous T cells must be distinguished by a validated CAR gate, transgene feature, clonotype evidence, or an explicitly stated limitation.
- **Panel B:** display distances in micrometres and summarize regions within patient or donor before cohort inference. A representative map must be paired with a specimen-level or patient-level summary and cannot stand alone as cohort evidence.
- **Panel C:** show membrane and soluble CXCL16 in separate panels and retain assay, unit, limit-of-quantification, and normalization information. Any soluble-to-surface ratio is exploratory.
- **Panel D:** separate chemotaxis, retention, and egress because they answer different biological questions. Low recovery cannot be called retention without viability, nonspecific loss, and barrier-integrity controls.
- **Panel E:** use the matched contrast “chemokine-module product minus identical control”, state the unit for every endpoint, show 95% confidence intervals, and avoid significance stars without exact adjusted P values in the accompanying table.

---

## Revised graphical abstract - Spatial chemokine logic and CAR-T outcome

**Figure type:** concise conceptual summary; replacement of the current graphical abstract after Figures 1-5 are revised.

**Caption intent:** integrate spatial chemokine organization, tumor-cell death competence, rational engineering, and biomarker-guided validation in one visual without implying that every pathway operates simultaneously.

**Evidence anchors:** Ozga et al., *Immunity* (2021), [doi:10.1016/j.immuni.2021.01.012](https://doi.org/10.1016/j.immuni.2021.01.012); Foeng et al., *Cell Reports Medicine* (2022), [doi:10.1016/j.xcrm.2022.100543](https://doi.org/10.1016/j.xcrm.2022.100543); Workenhe et al., *OncoImmunology* (2021), [doi:10.1080/2162402X.2021.1893466](https://doi.org/10.1080/2162402X.2021.1893466).

### Prompt

```text
Create a concise publication-quality graphical abstract for a peer-reviewed immunology review. Use a clean flat-vector BioRender style, white background, English labels only, and a wide aspect ratio of approximately 2.20:1. Use three connected stages and keep all text legible at single-column preview size. Use only short labels and one-line conclusions.

MAIN TITLE: Chemokine-guided CAR-T-cell responses in solid tumors

LEFT THIRD - SPATIAL CHEMOKINE LANDSCAPE:
Show a blood vessel, tumor nest, CAF-rich stroma, and peripheral lymphatic vessel. Place CAR-T cells at three possible locations: entering through activated endothelium, retained or excluded in stroma, and leaving through a lymphatic vessel. Label “entry”, “positioning / retention”, and “egress”. Show a blue CXCL9/CXCL10/CXCL11-CXCR3 route toward an inflamed tumor and a violet CXCL12-CXCR4 route with context-dependent recruitment, sequestration, fugetaxis, or egress. Add a small CXCL16 icon split into “membrane adhesion” and “soluble chemotaxis”. Add the sequential relay “NK cell - CCL5 / XCL1 -> cDC1 - CXCL9 / CXCL10 -> CXCR3-positive T cells”. Do not draw a direct NK-cell-to-T-cell bypass. Place “source x location x time x proteoform = functional cue” below the tissue.

CENTER THIRD - TUMOR-CELL DEATH COMPETENCE:
Show a CAR-T cell engaging an antigen-positive tumor cell and releasing perforin/granzyme B, IFN-gamma, and TNF. Split the tumor outcome into “death executed” and “death resisted”. Under death executed show conditional apoptosis, necroptosis, or pyroptosis, followed by DAMP and antigen release, dendritic-cell activation, and broader effector recruitment. Under death resisted show continued inflammatory chemokine production, repeated CAR-T contact, suppressive myeloid cells, and CAR-T dysfunction. Add the conclusion “Inflammation is productive only when spatial recruitment is coupled to effective tumor-cell killing.”

RIGHT THIRD - RATIONAL ENGINEERING AND VALIDATION:
Show three compact design icons: “matched receptor”, “engineered chemokine / local payload”, and “microenvironment combination”. Add a safety-control icon labeled “inducible expression / OFF or suicide switch”. Below show paired sample icons for final product, serial blood, and tumor biopsy leading to single-cell, spatial, protein, and functional migration assays. End with “biomarker-selected, controllable trafficking”. Add a small note: “Most solid-tumor strategies remain preclinical or early clinical; no chemokine-engineered CAR-T product is approved.”

COLOR AND LABEL RULES: CAR-T cells navy; endogenous T cells teal; tumor cells salmon; dendritic cells green; CAFs purple; suppressive myeloid cells orange; productive route green; chronic-failure route amber/red; engineered chemokines cyan-green; safety controls navy-gray.

EXCLUSIONS: Do not depict every axis as active simultaneously. Do not infer a functional gradient from bulk RNA or total protein. Do not show CXCL12-CXCR4 as uniformly suppressive. Do not call every lytic death immunogenic. Do not use CXCL13 as a stand-alone exhaustion marker. Do not claim clinical superiority, response rates, or quantitative improvement for any module.
```

## Final assembly checklist

- Keep the conceptual figures and quantitative Figure 6 visually distinct: conceptual arrows describe hypotheses or established mechanisms; quantitative panels show observed data and uncertainty.
- Cross-check every molecule, receptor, cell type, and direction of action against the final manuscript text and legend.
- Export conceptual artwork as editable SVG or PDF and archive the source file. Export final raster copies at 300 dpi or higher only after journal dimensions are fixed.
- Preserve panel letters and text as vector elements. Do not rasterize small labels.
- Record the figure source-file checksum and final export checksum in the project provenance table.
- For Figure 6, verify that every sample is traceable to the manifest and that `evidence_eligible=true`; otherwise keep the panel outside the manuscript.
- Remove synthetic warning labels only when the figure has been regenerated from verified non-synthetic inputs. Never cover, crop, or manually delete the warning from a synthetic output.
