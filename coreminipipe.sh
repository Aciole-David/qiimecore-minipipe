#!/bin/bash
#Several Core Microbiome Tests
echo "#Several Core Microbiome Tests"
echo 

run=0
source activate qiime1 #Ativar qiime1
echo "qiime1 activated"

for otutab in `ls *.biom`; do #otu-table

now=$(date) #printar data e hora
group=SC7-CC7 #nome do teste (depende de quais amostras sao usadas)
echo "factor of test: "$group
map=$group/MAP-$group.txt #arquivo de mapa
echo "map file: "$map
for s in SC CC
do
 echo "factors: " "$s"
done

echo "otu table: "$otutab
run=1

if [ $run -eq 1 ]; then
echo "run= $run"
for s in SC CC #loop para executar filtragem de amostras e calcular core-microbiome para todos os casos
do
	echo "$now #Filtering samples $s..." 
	filter_samples_from_otu_table.py -i ./$otutab -o ./$group/filtered-$s.biom -m ./$map -s $group:$s
	echo filter_samples_from_otu_table.py -i ./$otutab -o ./$group/filtered-$s.biom -m ./$map -s $group:$s
	echo "$now #Samples $s filtered..."
	echo 

	echo "$now #Computing $s core microbiome..."
	echo compute_core_microbiome.py -i ./$group/filtered-$s.biom -o ./$group/$s --mapping_fp ./$map --valid_states $group:$s
	compute_core_microbiome.py -i ./$group/filtered-$s.biom -o ./$group/$s --mapping_fp ./$map --valid_states $group:$s
  
	echo "$now #$s core microbiome computed..."
	echo
done
fi
run=2
echo "run= $run"

if [ $run -eq 2 ]; then #para ter certeza que vai rodar somente apos os comandos acima

 #pegar nome do primeiro caso, apenas para usar os nomes dos arquivos biom para pegar as fracoes (% de prevalencia em amostras)
 caso1=SC

 mkdir ./$group/$otutab-fractions
 f0=`ls -1v ./$group/$caso1/*.biom`
 for filenames in $f0; do
  f=`echo "$filenames"`
  frac=`echo "${f#*table_}"`
  fracsteps=`echo "${frac%.*}"`
  echo $fracsteps
  mkdir ./$group/$otutab-fractions/fracstep-$fracsteps #criar pasta para guardar arquivo lista usado na etapa de merge otu-tables
  a=`find ./$group/*/*$fracsteps.biom | tr '\n' ','` #listar, separando por virgula, as otu tables geradas para fazer merge delas
  echo $a >./$group/$otutab-fractions/fracstep-$fracsteps/$fracsteps.txt #guardar lista num arquivo .txt
  b=`rev ./$group/$otutab-fractions/fracstep-$fracsteps/$fracsteps.txt | cut -c2- |rev` #Remover ultimo caractere - virgula - e guardar na variavel "b"
  echo $b >./$group/$otutab-fractions/fracstep-$fracsteps/tables-to-merge-$fracsteps.txt
  rm ./$group/$otutab-fractions/fracstep-$fracsteps/$fracsteps.txt

 #echo "this is a: $a" #mostrar a lista com virgula no final
 #echo "this is b: $b" #MOSTRAR A LISTA SEM VIRGULA NO FINAL

  echo "$now #Merging core otu-tables $fracsteps..." #Unir tabelas core nessa fracao
  merge_otu_tables.py -i $b -o ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom
  echo merge_otu_tables.py -i $b -o ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom>./$group/$otutab-fractions/fracstep-$fracsteps/merge-command-$fracsteps.txt
  echo "$now #Core otu-tables $fracsteps merged..."
  echo

  echo "$now #Creating summary for $fracsteps"
  biom summarize-table -i ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom -o ./$group/$otutab-fractions/fracstep-$fracsteps/summary-$otutab-$group-$fracsteps-CORE.txt
  echo biom summarize-table -i ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom -o ./$group/$otutab-fractions/fracstep-$fracsteps/summary-$otutab-$group-$fracsteps-CORE.txt >./$group/$otutab-fractions/fracstep-$fracsteps/summary-command-$fracsteps.txt
  echo "$now #Summary $fracsteps created"
  echo 

  echo "$now #Converting otu-tables $fracsteps to tsv"
  biom convert -i ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom -o ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage"
  echo biom convert -i ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom -o ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.txt --to-tsv --header-key taxonomy --output-metadata-id "ConsensusLineage">./$group/$otutab-fractions/fracstep-$fracsteps/convert-command-$fracsteps.txt
  echo "$now #Merged table $fracsteps converted"
  echo 
  
  ################################################################################################################
  #Adding a second convert step to change taxonomy "k__Fungi" to "k__Bacteria" in order to use Microbiome Analyst!
  echo "$now #Second convertion of otu-tables $fracsteps to tsv"
  biom convert -i ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom -o ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.txt --to-tsv --header-key taxonomy
  echo biom convert -i ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom -o ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.txt --to-tsv --header-key taxonomy>./$group/$otutab-fractions/fracstep-$fracsteps/second-convert-command-$fracsteps.txt
  echo "$now #Second convertion of merged tables $fracsteps done"
  echo 
  
  #Now edit the txt with sed
  sed -i.bak 's/k__Fungi/k__Bacteria/g' ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.txt
  echo sed -i.bak 's/k__Fungi/k__Bacteria/g' ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.txt>./$group/$otutab-fractions/fracstep-$fracsteps/sed-command-$fracsteps.txt
  
  #Convert adjusted txt table back to biom format
  biom convert -i ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.txt -o ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy
  echo biom convert -i ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.txt -o ./$group/$otutab-fractions/fracstep-$fracsteps/$otutab-$group-$fracsteps-CORE.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy>./$group/$otutab-fractions/fracstep-$fracsteps/last-convertion-command-$fracsteps.txt
  ################################################################################################################

  echo "$now #Core computing $fracsteps done"
  
	
  echo
 done

fi
run=3
echo "run= $run"



if [ $run -eq 3 ]; then #para ter certeza que vai rodar somente apos os comandos acima
 echo
 echo "$now #Core computing for all prevalence fractions done"
 echo
fi

for s in SC CC; do
	rm ./$group/*.biom
	rm -R ./$group/$s/;
done

done
