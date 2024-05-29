# ScRNA-project
## 1. Script structure
**Load students ttest pvalue and pearson correlation value file to MySQL:**
load_file.py
MySQL.py
**Merge gene pairs' pvalue and correlation value as long table and output as .arrow files**:
jointable.py
**Filter gene pairs according to pavlue, correlation (and jaccard index which calculate with stringdb real protein interaction):**
Filter.py
**Data visualization:**
DrawTool.py
drawNetwork.py
**Common Tool**
utils folder

## 2. figure folder Structure
```
corr_pval
├── tissue_level_40-60Age_corr{corr_value}_p{p_value}  # tissue level data from <40 and >60 years-old donars
│   ├── {study}|{group}|{tissue}.graphml  # correlation network, graphml format
│   ├── {study}|{group}|{tissue}.pdf # correlation network, pdf format
│   ├── agingMarkers-related_tissue_heatmap_c.png  # *_c.png->heatmap contains full organism, full markers and top candidates
│   ├── agingMarkers-related_tissue_heatmap_f.png  # *_f.png->heatmap contains full organism, full markers
│   ├── agingMarkers-related_tissue_heatmap_m.png  # *_m.png->heatmap contains full markers
│   ├── agingMarkers-related_tissue_heatmap_mc.png  # *_mc.png->heatmap contains full markers and top candidates
│   ├── marker_degree.csv  # connection degrees screened by correlation shrefold and p-value shrefold
│   └── temp_result.csv  # reshaped data used to construct network
├── tissue_level_medianAge_corr{corr_value}_p{p_value}
├── cell_level_40-60Age_corr{corr_value}_p{p_value}  # single-cell level data from <40 and >60 years-old donars
│   ├── {study}|{group}|{tissue}|{cellType}.graphml  # correlation network, graphml format
│   ├── {study}|{group}|{tissue}|{cellType}.pdf # correlation network, pdf format
│   ├── agingMarkers-related_tissue_heatmap_{c/f/m/mc}.png # Summation cellType counts by tissues
│   ├── agingMarkers-related_cell_heatmap_{c/f/m/mc}.png # Take the average of tissue counts by cellType 
│   ├── agingMarkers-related_tissueCell_heatmap_{c/f/m/mc}.png  # Show as tissue_cellType, same tissue_cellType from different studies would conduct summation
│   ├── marker_degree.csv 
│   └── temp_result.csv
└── cell_level_medianAge_corr{corr_value}_p{p_value}
```
