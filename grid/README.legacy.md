
# Methods notes: automatic ancestry

The preferred PRS-CSx replication uses a classifier trained with labelled reference individuals. 
The files supplied with DiscoDivas include projected-PCA weights and four population medians, but not the complete labelled individual-level 1KG PC training table. 
Therefore v9 implements an automatic, reproducible UKB-specific approximation rather than silently treating nearest-center assignment as formal ancestry.

The automatic classifier uses two independent signals:

1. released 1KG AFR/EAS/EUR/SAS centers in the DiscoDivas PC coordinate system;
2. UKB self-reported ethnicity, used only to identify high-purity training anchors when it agrees with the nearest center.

Balanced LDA is trained on projected PC1-PC20. Equal class priors prevent the very large EUR group from dominating. 
Posterior probability below 0.90 produces OTH, retaining ambiguous and admixed participants for DiscoDivas continuum analyses. 
High-confidence, genetic/self-report-concordant participants are used for fine-tuning cohorts in evaluation.

This is methodologically stronger than nearest-center grouping and requires no manually curated ancestry file. 
It is not claimed to be numerically identical to the random-forest classifier trained directly on labelled 1KG individuals in the original PRS-CSx paper. 
A user-supplied validated ancestry file still overrides the automatic procedure when exact external classification is available.
