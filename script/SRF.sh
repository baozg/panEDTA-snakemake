#!/bin/bash
fa=$1
kmer=151
threads=$2
tel="TTTAGGG"

srf="/ebio/abt6/zbao/software/srf"
kmc=
tidk=
seqkit=
barrnap=
minimap2=


sample=$(basename -s .fa $fa)

outpre=$sample\_k${kmer}

# kmc
mkdir -p $sample.fa.mod.EDTA.raw/SRF
cd $sample.fa.mod.EDTA.raw/SRF
ln -s ../../$fa ./
mkdir -p tmp/${sample}_SRF

$kmc/kmc -fm -k$kmer -t$threads -ci10 -cs100000 $fa ${outpre}_count.kmc tmp/${sample}_SRF
$kmc/kmc_dump ${outpre}_count.kmc ${outpre}_count.txt

# srf
${srf}/srf -p ${sample} ${outpre}_count.txt > ${outpre}_SRF.fa
${srf}/srfutils.js enlong ${outpre}_SRF.fa > ${outpre}_SRF.enlong.fa

# analyze

## telomere
tellen=${#tel}
echo $tellen
$tidk/tidk search -d ./ -s $tel -o ${sample} ${outpre}_SRF.enlong.fa
sed "1d" ${sample}_telomeric_repeat_windows.tsv|awk -v tellen=$tellen '(($3+$4)*tellen)/$2>0.8' |cut -f1 > tel.SRF.list

## rDNA
$seqkit/seqkit grep -v -r -f tel.SRF.list ${outpre}_SRF.enlong.fa > ${outpre}_SRF.enlong_exclude_tel.fa 
$barrnap/barrnap -t $threads --kingdom euk ${outpre}_SRF.enlong_exclude_tel.fa > ${outpre}_SRF.enlong_exclude_tel.rdna.gff3
grep -v "^#" ${outpre}_SRF.enlong_exclude_tel.rdna.gff3|cut -f1|sort|uniq > rDNA.SRF.list

## TE
$seqkit/seqkit grep -v -r -f rDNA.SRF.list ${outpre}_SRF.enlong_exclude_tel.fa > ${outpre}_SRF.enlong_exclude_tel_rdna.fa
$tesorter/TEsorter -dp2 -p ${threads} ${outpre}_SRF.enlong_exclude_tel_rdna.fa
grep -v "^#" ${outpre}_SRF.enlong_exclude_tel_rdna.fa.rexdb.cls.tsv|cut -f1|sort|uniq > TE.SRF.list

# SRF filter
$seqkit/seqkit grep -v -r -f TE.SRF.list ${outpre}_SRF.enlong.fa > ${outpre}_SRF.enlong.TEfree.fa

$minimap2/minimap2 -t $threads -c -N1000000 -f1000 -r100,100 ${outpre}_SRF.enlong.TEfree.fa $fa > ${outpre}_SRF_enlong.TEfree.paf
${srf}/srfutils.js paf2bed -l 50 ${outpre}_SRF_enlong.TEfree.paf > ${outpre}_SRF_enlong.TEfree.bed
${srf}/srfutils.js filter -l 50 -c 0.9 -d 0.9 ${outpre}_SRF_enlong.TEfree.paf > ${outpre}_SRF_enlong.TEfree.filter.paf
sort -k1,1V -k3,3n ${outpre}_SRF_enlong.TEfree.filter.paf|cut -f6|sort|uniq -c|awk '$1 > 50'|awk '{print $2}' > final.include.SRF.list

# masking
grep -f final.include.SRF.list ${outpre}_SRF_enlong.TEfree.bed > ${outpre}_SRF.mask.bed

$bedtools/bedtools maskfasta -fi $fa -bed ${outpre}_SRF.mask.bed -fo $sample.masked.SRF.fa

# EDTA_raw (unmask for LTR and Helitron, mask for TIR)

