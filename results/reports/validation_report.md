# Validation Report

## Execution Context

- Git commit: `140272f`
- Git status: `clean except ignored local artifacts`
- Python: `3.13.7`
- Platform: `macOS-26.5.1-arm64-arm-64bit-Mach-O`

## Primary Model

- Samples: 9 total.
- Infection counts: {'VacV': 6, 'mock': 3}.
- Model-ready genes: 21835.
- Count/metadata alignment: True.
- Full DGE rows: 21835.
- Required DGE columns present: True.
- Translation-factor matches: 465.
- Significant translation-factor impacts at padj < 0.05: 184.

## DHX29 Primary Signal

- Gene ID: `ENSG00000067248.11`.
- log2 fold change, VacV versus mock: 0.974837.
- Adjusted p-value: 1.56335e-05.

## Metadata Conflict Handling

The primary model uses supplementary count filenames as the fraction source of truth. The following GEO title/file conflicts were recorded and retained in the manifest:

- `GSM8544828`: `title_fraction=polysome;filename_fraction=total`
- `GSM8544830`: `title_fraction=polysome;filename_fraction=total`

## Strict-Title Sensitivity Model

- Excluded samples: ['GSM8544828', 'GSM8544830'].
- Samples: 7 total.
- Infection counts: {'VacV': 5, 'mock': 2}.
- Model-ready genes: 20777.
- Significant strict-title translation-factor impacts at padj < 0.05: 129.
- DHX29 primary log2FC/padj: 0.974837 / 1.56335e-05.
- DHX29 strict-title log2FC/padj: 0.912853 / 0.00935895.
- Direction agreement among translation factors significant in either primary or strict model: 1.000.

## Cross-Dataset Novelty Pilot

- Meta-signature rows: 607.
- Pilot registry rows: 2.
- `GSE287860`: 23151 genes modeled; 554 target genes; 158 significant targets.
- `GSE288000_NTC`: 21722 genes modeled; 354 target genes; 131 significant targets.
- DHX29 dataset count: 3.
- DHX29 significant dataset count: 2.
- DHX29 mean log2 fold change across pilot contexts: 0.543678.
- DHX29 direction consistency: 1.000.

## Mechanistic Evidence Layer

- Host model genes: 21617.
- Viral model genes excluded from host-only curation: 218.
- Host protein-coding model genes by GENCODE v44: 14148.
- Curated host protein-coding translation targets: 197.
- Curated significant targets with abs(log2FC) >= 0.5: 140.
- DHX29 overall mechanistic priority rank: 52.
- DHX29 DHX/DDX helicase priority rank: 4.
- Minimum curated module enrichment padj: 1.35958e-85.
- Top mechanistic candidate: `EIF3L` (score 9.188; log2FC -2.230; padj 3.82e-11).
- DHX29 viral-load correlation in infected samples: Pearson r=0.779, p=0.0677; exploratory because n=6 infected samples.
- Strongest curated module enrichment: `translation_targets_significant_with_effect` / `60S_ribosomal_subunits` (overlap 43; padj 1.36e-85).

## Figure Checks

- `results/figures/translation_heatmap.pdf`: 31397 bytes; SHA256 `d89c7fbc51ac5b66da4dbbdaa5e3a1ededea474fdbc3708dd376d43898b894a6`.
- `results/figures/translation_heatmap.png`: 442218 bytes; SHA256 `3ccb6291d5dc6242f91739d9bafbbd73e2f4ffa92c04183a653f9d31f9421616`.
- `results/figures/translation_heatmap.svg`: 169164 bytes; SHA256 `e1c8988cc1c9611fac432abb756959fec5ef20b3df4e159b07edbe070139490a`.
- `results/figures/volcano_plot.pdf`: 380911 bytes; SHA256 `45e475e5b43dd4077146edd0555227de8dbeb55e20a0ba51ae1ca804b7d54ef8`.
- `results/figures/volcano_plot.png`: 382407 bytes; SHA256 `2a9931d5d05fa86891490b4c85b5fa49ef8e2ea3e5c391a8ee4cf9d3fd3b80a6`.
- `results/figures/volcano_plot.svg`: 3566085 bytes; SHA256 `88b35da6ab81ceeaad69c9d2cf8350f20ba835d48f14753a16291a470deaa538`.
- `results/meta/poxvirus_translation_factor_meta_signature.pdf`: 26387 bytes; SHA256 `d79dcb8b2d689ba54b3bce7eab6229c53a9151554479ef6b05d41db14213cd90`.
- `results/meta/poxvirus_translation_factor_meta_signature.png`: 167004 bytes; SHA256 `6d490a91f4de2f801bf95b2bade66c0bee4b04b8f8700199de3396ad9c4a43c6`.
- `results/meta/poxvirus_translation_factor_meta_signature.svg`: 88662 bytes; SHA256 `4e9d6a58b3116f5586aa2dbcd3f47c7fd99cccb74e0a47609bf448ec78c0fd86`.
- `results/mechanistic/candidate_evidence_matrix.pdf`: 26314 bytes; SHA256 `a8be70e04268fb25287e6c2b30fd483342bb46686cf895ce4ea351a517efb8e1`.
- `results/mechanistic/candidate_evidence_matrix.png`: 183679 bytes; SHA256 `48c96edd52c21090ef7da2f814913c05c76678016b57c48efe750b30ba0500a7`.
- `results/mechanistic/candidate_evidence_matrix.svg`: 95674 bytes; SHA256 `dbd247aa24f45f7162d0369c17e49aff26e32cffc7fa2e915cd0e03ed77307b7`.
- `results/mechanistic/dhx29_viral_load_correlation.pdf`: 22029 bytes; SHA256 `fa4a105c01180e8110b4a69fba50f65cc443c16da9af173734518ce2f92f7936`.
- `results/mechanistic/dhx29_viral_load_correlation.png`: 163804 bytes; SHA256 `e8b3d94429a695cd4c39c7a5fd3f13f3e1a5178f75b012633f4d35059f9ffff0`.
- `results/mechanistic/dhx29_viral_load_correlation.svg`: 53820 bytes; SHA256 `0784dc4af55efd6c648a96020e197d5f09a564f52dbc6804d7cd659ce00bbef8`.
- `results/mechanistic/dhx_ddx_helicase_priority.pdf`: 17554 bytes; SHA256 `1025d29a33ba9ac4f469400b8225aa2fc60662f6252a56d4ef8935d41c9f4bd7`.
- `results/mechanistic/dhx_ddx_helicase_priority.png`: 275816 bytes; SHA256 `53195feb2b6abf56d656ab0ec1d1b9cddf3de1c1e27e145224aec2f0191e6b68`.
- `results/mechanistic/dhx_ddx_helicase_priority.svg`: 95779 bytes; SHA256 `f1512418bb4234d456479f356fd733ccd4bfbc57c830070f326772164884e6a3`.
- `results/mechanistic/host_only_curated_volcano.pdf`: 343544 bytes; SHA256 `97561f06bd572428316dbcf27630dc491a8b31a427c20b86f4163f932413165d`.
- `results/mechanistic/host_only_curated_volcano.png`: 376334 bytes; SHA256 `365d8ec4af50c1d3fc5b917ed62f43b5f6e21b71bd1e58e83388976aa61b8c51`.
- `results/mechanistic/host_only_curated_volcano.svg`: 2435488 bytes; SHA256 `b9cb2339e8921c7bd27c474256a5b8e1c530421156b2fcae6487b01eaafdf93f`.
- `results/mechanistic/module_enrichment_barplot.pdf`: 19198 bytes; SHA256 `40d95bb45964e29e04cb8a1e389c5a6794ae485c461ab2d4878fd1edde48229c`.
- `results/mechanistic/module_enrichment_barplot.png`: 236579 bytes; SHA256 `20cda4d00b2a2fb4fbec4b6c2f9c9a2c1b6e0f843d60dca8530408c92df499bd`.
- `results/mechanistic/module_enrichment_barplot.svg`: 103253 bytes; SHA256 `db99dccfbe83a09b2288e1a46f4737b15b9726d09b1cbbfdd8a430363f3c9db4`.
- `results/mechanistic/string_priority_network.pdf`: 29738 bytes; SHA256 `928d81a57015ebc9b49ed0658464698c81082d4a0890f7a21758425fb9e55ab9`.
- `results/mechanistic/string_priority_network.png`: 1524929 bytes; SHA256 `198781269ee00cdec1c3cee2ef98e9d5a50bd0af13c6d61b2b979bab5ae08c34`.
- `results/mechanistic/string_priority_network.svg`: 107117 bytes; SHA256 `8862d118391f3718498b6c608fe1ba0d3e468cdaeee570e8f53e996c343d1c8b`.

## Hash Receipt

A machine-readable hash table was written to `results/reports/file_hashes_sha256.csv`.

## Interpretation Boundary

These outputs support a transcript-level secondary analysis of a public dataset. They do not establish direct protein-level mechanism without additional translational, proteomic, or perturbation validation.
