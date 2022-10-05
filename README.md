# 7 Letter Scoring Scheme

The full paper can be found [here](https://github.com/SimonJunyouLi/7-letter-scoring-scheme/blob/main/7%20letter%20scoring%20scheme.pdf). 

# **Motives**

Although the COVID-19 vaccines are "95% effective", they often have to be stored under extreme conditions. Moreover, certain vaccinated patients have shown signals of blood clots. Therefore, it is imperative to improve the stability and safety of the vaccine. 

However, there is no publicly available research on improving vaccine designs. Nevertheless, research in related fields can be applied to optimize vaccines. 

Thus, this project aims to investigate whether [the Nussinov and Jacobson model (1980)](https://www.pnas.org/content/pnas/77/11/6309.full.pdf), implemented with the 7-letter scoring scheme suggested by Prof. Forbes J. Burkowski, is a viable method to develop a mRNA vaccine.

# **Data**
The [Wuahn reference genome](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?VirusLineage_ss=Severe%20acute%20respiratory%20syndrome%20coronavirus%202%20(SARS-CoV-2),%20taxid:2697049&SeqType_s=Nucleotide) was chosen as the reference sequence to be modified in this study, while the Pfizer vaccine 
was used as a comparison for the resulting sequence. 

# **Methodology**
This project evalutes how successful the mRNA vaccine desing is the by two standard, the [free energy](https://en.wikipedia.org/wiki/Thermodynamic_free_energy) of the sequence [(MFE)](http://eternawiki.org/wiki/index.php5/Minimum_Free_Energy_Structure) and the [codon adaptation index](https://en.wikipedia.org/wiki/Codon_Adaptation_Index) (CAI). For an optimized mRNA sequence, one must have a minimized MFE for the stability of the sequence and a maximized CAI for the maximum effective codons, i.e. the effectiveness of the sequence. 

The result of the method was tested with the [RNAfold WebServer](http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi) and [EMBOSS](https://www.bioinformatics.nl/cgi-bin/emboss/cai) and compared with the Pfizer mRNA vaccine.

# **Result**
![image](https://user-images.githubusercontent.com/47229668/155698096-8b2db12c-9042-4e43-982b-c1216804493c.png)

The resulting sequence achieved a lower MFE than both the reference sequence and the Pfizer sequence. However, the method was less successful in maximizing CAI as the resulting sequence's CAI score was lower than Pfizer's but higher than that of the reference sequence. The resulting sequence's CAI score of 0.829 is more than adequate, as [“CAI of >0.8 is regarded as good”](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4724736/).
