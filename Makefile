.PHONY: all env data signatures lincs gsfm annotation ranking clean

# Run the full pipeline end-to-end
all: signatures lincs gsfm annotation ranking

# Set up conda environment
env:
	conda env create -f environment.yml
	conda run -n exercise-mimetics pip install -e .

# Extract MoTrPAC data via R (requires MotrpacRatTraining6moData)
data:
	Rscript scripts/extract_motrpac.R --outdir data/raw --timepoints 8w

# Notebook 01 — MoTrPAC signatures
signatures:
	jupyter nbconvert --to notebook --execute \
		--ExecutePreprocessor.timeout=600 \
		--output notebooks/01_motrpac_signatures_executed.ipynb \
		notebooks/01_motrpac_signatures.ipynb

# Notebook 02 — LINCS query
lincs: signatures
	jupyter nbconvert --to notebook --execute \
		--ExecutePreprocessor.timeout=1200 \
		--output notebooks/02_lincs_query_executed.ipynb \
		notebooks/02_lincs_query.ipynb

# Notebook 03 — GSFM embedding
gsfm: lincs
	jupyter nbconvert --to notebook --execute \
		--ExecutePreprocessor.timeout=3600 \
		--output notebooks/03_gsfm_embedding_executed.ipynb \
		notebooks/03_gsfm_embedding.ipynb

# Notebook 04 — Annotation
annotation: lincs
	jupyter nbconvert --to notebook --execute \
		--ExecutePreprocessor.timeout=1200 \
		--output notebooks/04_annotation_executed.ipynb \
		notebooks/04_annotation.ipynb

# Notebook 05 — Final ranking
ranking: lincs gsfm annotation
	jupyter nbconvert --to notebook --execute \
		--ExecutePreprocessor.timeout=600 \
		--output notebooks/05_final_ranking_executed.ipynb \
		notebooks/05_final_ranking.ipynb

# Remove generated outputs (not raw data)
clean:
	rm -f notebooks/*_executed.ipynb
	rm -f data/processed/*.parquet
	rm -f data/processed/*.txt
	rm -f results/figures/*.png
	rm -f results/*.csv

# Remove everything including downloaded data
clean-all: clean
	rm -rf data/raw/*
	rm -rf data/external/*
