outdir=tmp
k2db=gg2_db
# rm -rf $k2db
mkdir -p $k2db $k2db/taxonomy $k2db/library

cp $outdir/*nodes.dmp $k2db/taxonomy/nodes.dmp
cp $outdir/*names.dmp $k2db/taxonomy/names.dmp
cp $outdir/*seqs.fna $k2db/library/seqs.fna

# fix delimiters! 
sed -i -e 's/$/\t|/g' $k2db/taxonomy/names.dmp
sed -i -e 's/$/\t|/g' $k2db/taxonomy/nodes.dmp

echo 'Kraken build'
kraken2-build --add-to-library $k2db/library/seqs.fna --db $k2db --no-masking

# full thing takes ~90 mins
kraken2-build --build --db $k2db --threads 224


# look
kraken2-inspect --db $k2db | head -n20

