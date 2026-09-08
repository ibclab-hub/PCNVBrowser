# PCNVBrowser
**PCNVBrowser** is a web-based platform for exploring, analyzing, and visualizing copy number variations (CNVs) across populations.  
It integrates CNV datasets from multiple tools and populations, enabling dynamic comparison, variant annotation, and enrichment analysis.

#### [http://biweb.konkuk.ac.kr/PCNVBrowser/](http://biweb.konkuk.ac.kr/PCNVBrowser/).
<p align="center">
<img width=800 alt="image" src="https://github.com/user-attachments/assets/61dcddfd-dcaa-4070-997d-5bfd5814fd70" />
</p>

---

## This repository contains:

- `CNV_calling_pipelines/` — CNV calling workflows for the four callers (cn.MOPS, CNVkit, Control-FREEC, readDepth)
- `CNV_Integration_Pipeline/` — integration of per-tool calls into population-level call-set support frequencies
- `Carrier_Frequency/` — individual-based carrier frequency calculation (DEL/DUP/ALL, tool-support thresholds)
- `Validation_Concordance/` — cross-caller and cross-resolution concordance, external truth-set benchmarking, and 1000 Genomes population-level concordance
- `Applications/` — analysis scripts for Applications 1-3 (chr8 carrier-frequency comparison, PD CNV filtering, Korean-specific CNV comparison)
- `Application_data/` — CNV datasets used as inputs for Applications 2 and 3
- `utils/` — general-purpose utility for summarizing user-defined region overlap across continental population groups

---

## License & Data Policy

PCNVBrowser is freely available for both academic and commercial use under the CC BY 4.0 License.
No user data is stored on the server; all uploaded content is discarded upon session termination.

---

## Contact

For questions or contributions, contact:

- Nayoung Park: `p3159@konkuk.ac.kr`  
- Dayeon Kim: `dayeonkim@hufs.ac.kr `
