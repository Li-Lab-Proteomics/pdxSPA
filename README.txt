pdxSPA v1.0 : Shared-peptides Allocation Methods for Patient derived xenograft (PDX)

Introduction:
pdxSPA includes three algorithm for quantifying the MS data from the PDX models, especially focusing on the quantification of the shared peptides. X-SPA was to assign the shared peptides based on the ratio of the sum of human-unique (HU) and mouse-unique (HU)MU peptide intensity at individual protein level, only if both HU and MU peptides were identified.The ratio will be changed when different protein is considered with X-SPA. The second one was S-SPA, using X-SPA when HU or MU peptides were identified; otherwise, the shared peptides were assigned based on the over the overall ratio of the sum of HU and MU peptide intensity. The third one is SPA, by which the shared peptides were assigned based on the overall ratio of the sum of HU and MU peptide intensity.

This readme file covers the following topics:
---------------------------------------------
1. Installation
git clone https://github.com/cathy0237/pdxSPA.git
cd pdxSPA/
pip install .


2.Running pxdSPA
Usage: pdxSPA [OPTIONS]

Options:
  -i, --input_file PATH  Pass in the evidence file after quantifile PDX samples by MaxQuant.
  -o, --outdir PATH      Pass in the output directory.  [default: .]
  -m, --method INTEGER   Choose methods between 0 (SPA), 1 (X-SPA) and 2 (S-SPA).  [default: 0]
  --version              Show the version and exit.
  -h, --help             Show this message and exit.

Example:
Runs pdxSPA on the evidence file from MaxQuant. SPA was used for shared-peptide allocation.
pdxSPA -i evidence.txt -m 0


3.pdxSPA output
The output gives the following result files. 
Summary file - gives the summary on peptides allocation and ratio. 
Expression file - gives the HUMAN protein expression in PDX by using pdxSPA.




