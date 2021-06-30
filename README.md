### Manuscript information (drafting phase)
---

###### about..

Gene-expression regulation is key to understanding species capacity to adapt and acclimate to climate-change. Post-larval acclimation to hypercapnic seawater improves the growth and oxidative status of juvenile Pacific geoduck Panopea generosa ([Gurr et al. 2021; pulished in Journal of Experimental Biology](https://pubmed.ncbi.nlm.nih.gov/34027545/)),  opening interest in the transcriptome profiles underpinning of modified phenotypes.

**Geoduck samples**
- Hatchery propegated Pacific geoduck clams *Panopea generosa*
- **Location**: Jamestown Point Whitney Shellfish Hatchery (Jamestown Sk'llalam Tribe), Brinnon, Washingtown (USA)

- **Age**: pediveliger to juvenile stage (30 days to 3-4 months post-fertlization) with different stress histories under repeated reciprocal hypercapnia.


**Within this repo...**

  (1) Stepwise pipeline for the raw read filtering, mapping, and counts (using.. fastp, multiQC, Stringtie2, HISAT2, prep.py). Inludes a detailed R markdown

  (2) Co-expression network analysis (WGCNA) and pairwise differential expression (DESeq2) of filtered read matrices in R

  (3) Gene function (goseq and goSlim) and pathway (KEGGprofiler) enrichment analysis in R

## Title (*draft...*)
---

### Experience-mediated divergene of transcriptome profiles under intermittent hypercapnia


## Authors
---
Samuel J. Gurr<sup>1</sup>, Shelly A. Trigg<sup>3</sup>, Brent Vadopalas<sup>2</sup>, Steven B. Roberts2<sup>3</sup>, and Hollie M. Putnam<sup>1</sup>

<sup>1</sup>University of Rhode Island, College of the Environment and Life Sciences, 120 Flagg Rd, Kingston, RI 02881, USA.

 <sup>2</sup>University of Washington, School of Aquatic and Fishery Sciences, 1122 NE Boat St, Seattle, WA 98105, USA.

<sup>3</sup>University of Washington, Washington Sea Grant, 3716 Brooklyn Ave NE, Seattle, WA 98105, USA.

## Abstract (*draft...*)
---

Gene-expression regulation is key to understanding species’ capacity to adapt and/or acclimate to climate-change. Post-larval acclimation to hypercapnic seawater improves growth and oxidative status of juvenile Pacific geoduck Panopea generosa, opening interest in the transcriptome profiles underpinning modified phenotypes. Following three-months of conditioning immediately post-settlement under ambient and moderately-elevated pCO2, repeated hypercapnia and an ambient depuration period elicited variation of transcriptome profiles between stress-acclimated and naïve juvenile geoducks. Representative of the main finding, stress-acclimated geoducks were enriched for quality control of mitochondria, immune defence during hypercapnia, and opportunistically increased transcripts involved in energy metabolism and biosynthesis during ambient recovery. Furthermore, continuous gene ontology enrichment included histone methyltransferases and transcription factors, illustrating that moderate-stress history may fine-tune transcription. In contrast, the naïve profile showed greater transcriptional demand and continuous enrichment for fatty-acid degradation and glutathione components suggesting unsustainable energetic requirements if changes in carbonate chemistry exacerbated or persisted. Altogether, transcriptomic findings complement the emergent physiological phenotypes in Gurr et al. (2021), highlighting beneficial gene-expression regulation and cellular maintenance by the acclimatized phenotype as opposed to putative depletion of endogenous fuels to supply broad transcription in absence of prior experience. Post-larval acclimatory periods, predecessor to episodic changes, can enhance robustness to environmental stress in juvenile P. generosa.


#### more..
---

##### Raw TagSeq data on NCBI:
[Accession #: PRJNA740307; "Transcriptome profiles of Panopea generosa under hypercapnic seawater"](https://www.ncbi.nlm.nih.gov/bioproject/740307)

##### Experimental design and phenotype manuscript:
[JEB, Gurr et al. 2021](https://pubmed.ncbi.nlm.nih.gov/34027545/)

##### RNA Extraction protocol:
[online notebook
post](https://samgurr.github.io/SamJGurr_Lab_Notebook/Updated-protocol-DNA-RNA-Extraction-of-geoduck-samples-(Zymo-kit)/)
