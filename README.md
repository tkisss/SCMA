# SCMA: Single-cell metabolomics analysis 

## Installation  	
#### install from PyPI

    python -m pip install --index-url https://test.pypi.org/simple/ --no-deps scma==0.1.1
    
#### install from GitHub

	git clone https://github.com/tkisss/SCMA.git
	cd SCMA
	python setup.py install
    
## Quick Start

### Command line

    scma.py -f INPUT_FILE -o OUTPUT_PATH
    
    INPUT_FILE: file of single-cell metabolomics, each two columns are a group which contains the ratio of nuclear to mass and its signal value, cell names should be groupname-xx, for example: A549-1,A549-2, gefitinib-A549-6,gefitinib-A549-8..., name before the last '-' will be considered as the group name, so the groups are A549 and gefitinib-A549.
    OUTPUT_PATH: output path

#### other parameters 
* ppm_threshold_peak: peak error threshold for peak selection in the same cell. Default:10
* ppm_threshold_cell: peak error threshold when combining different cells. Default: 20
* decrease: peak selection. Default:True
* peak: whether the input file has selected peak. Default:False
* filter: peaks appearing in less than a% cells will be filtered. Default:0.5
* knn: use knn for missing value imputation. Default:True
* n_neighbors: knn algorithm parameter. Default:5
* method: dimensional reduction methods, including 'PLS','PCA','UMAP','tSNE'. Default:'PLS'
* p_value: differential analysis parameter. Default:0.05
* log2fold: differential analysis parameter. Default:0.5

#### Output
Output will be saved in the output folder including:
* **processed.csv**: preprocessed data after peak selection and coordinate alignment
* **filtered.csv**: results after peak filtering, peaks appearing in less than a% cells will be filtered (parameter: --filter)
* **knn.csv**: results after knn imputation
* **PLS.txt**: dimensional reduction result
* **de.csv**: differential peaks among different groups (parameter: --p_value and --log2fold)
* **violinplot.png**: violin plot correspond to processed.csv, show how many cells each peak involved in
* **embedding.pdf**: dimensional reduction embedding
* **heatmap.png**: heatmap of differential peaks among different groups

#### Help
Look for more usage of scMA

	scma.py --help 