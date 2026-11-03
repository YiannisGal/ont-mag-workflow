#!/usr/bin/env bash
# One-off database downloads. Run inside an environment that has each tool,
# e.g.  conda run -n <env> bash ...  or copy the lines you need.
# Then put the resulting paths into config/config.yaml -> databases.
# Sizes are approximate; GTDB-Tk, Kraken2 and Bakta are the big ones.
set -euo pipefail
DB=${1:?usage: download_databases.sh /path/to/databases}
mkdir -p "$DB"; cd "$DB"

# CheckM2 (~3 GB)
checkm2 database --download --path "$DB/checkm2"

# GUNC proGenomes 2.1 (~13 GB)
mkdir -p gunc && gunc download_db gunc

# GTDB-Tk R232 (~110 GB). Check https://ecogenomics.github.io/GTDBTk/ for the current URL
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz
# tar xzf gtdbtk_r232_data.tar.gz

# KOfam (~4 GB)
mkdir -p kofam && cd kofam
wget -nc https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
wget -nc https://www.genome.jp/ftp/db/kofam/ko_list.gz
tar xzf profiles.tar.gz && gunzip -f ko_list.gz
cd "$DB"

# dbCAN (~2 GB)
dbcan_build --cpus 8 --db-dir "$DB/dbcan" --clean

# geNomad (~1.5 GB)
genomad download-database "$DB"

# SingleM default metapackage
singlem data --output-directory "$DB/singlem"

# Optional
# bakta_db download --output "$DB/bakta" --type full
# download-antismash-databases --database-dir "$DB/antismash"
# Kraken2/Bracken: use a prebuilt index, e.g. https://benlangmead.github.io/aws-indexes/k2
# Kaiju proGenomes: https://bioinformatics-centre.github.io/kaiju/downloads.html
