#!/bin/bash
###############################################################################
# PERFORMANCE OF FALCON HANDLING INCREASING SUBSTITUTIONAL MUTATIONS RATES
# REQUIREMENTS: linux OS with cmake, git and gnuplot
# (sudo apt-get install cmake git gnuplot)
###############################################################################
# PARAMETERS ==================================================================
INSTALL=1;
SIMULATE=1;
SHUFFLE=1;
FALCON=1;
MUMMER=1;
MUMMER20=1;
FILTER=1;
PLOT=1;
#==============================================================================
MLIMIT=41;
FQLINE=200;
FQNREADS=20000;
#==============================================================================
DISTRIBUTION="0.3,0.2,0.2,0.3,0.000";
EXTRAMUT=" ";
#EXTRAMUT=" -ir 0.01 -dr 0.01 ";
#FPARAM=" -m 14:500:1:5/50 -c 15 -g 0.9 ";
FPARAM=" -m 20:500:1:5/50 -m 14:100:1:0/0 -m 12:1:0:0/0 -m 4:1:0:0/0 \
#-c 15 -g 0.95 ";
###############################################################################
if [[ "$INSTALL" -eq "1" ]]; then
# CLEAN & INSTALL =============================================================
rm -fr goose-* FALCON XS goose/ falcon/ xs/ SAMPLE* DB-mfa;
# GET GOOSE FRAMEWORK =========================================================
git clone https://github.com/pratas/goose.git
cd goose/src/
make
cd ../../
cp goose/src/goose-* .
cp goose/scripts/ShufFASTQReads.sh .
cp goose/scripts/GlobalMUMmer.sh .
# GET FALCON ==================================================================
git clone https://github.com/pratas/falcon.git
cd falcon/src/
cmake .
make
cp FALCON ../../
cd ../../
# GET XS ======================================================================
git clone https://github.com/pratas/xs.git
cd xs
make
cp XS ../
cd ..
# GET MUMmer 3.23 =============================================================
rm -rf MUMmer.tar.gz MUMmer3.23/
wget -O MUMmer.tar.gz \
http://downloads.sourceforge.net/project/mummer/mummer/3.23/MUMmer3.23.tar.gz
tar -xvzf MUMmer.tar.gz
cd MUMmer3.23
make check
make install
cd ..
cp MUMmer3.23/nucmer .
cp MUMmer3.23/show-coords .
cp MUMmer3.23/delta-filter .
fi
###############################################################################
# SIMULATE ====================================================================
if [[ "$SIMULATE" -eq "1" ]]; then
./XS -v -ls $FQLINE -n $FQNREADS -f $DISTRIBUTION -s 0 SAMPLE.fq
# MUTATE ======================================================================
./goose-fastq2fasta < SAMPLE.fq > SAMPLE.fa
./goose-fasta2seq   < SAMPLE.fa > SAMPLE
./goose-seq2fasta -n "Substitution0" < SAMPLE > SAMPLE0.fa
cat SAMPLE0.fa > DB.mfa;
for((x=1 ; x<$MLIMIT ; ++x));
  do
  MRATE=`echo "scale=3;$x/100" | bc -l`;
  echo "Substitutions rate: $MRATE";
  ./goose-mutatedna -mr $MRATE $EXTRAMUT < SAMPLE > SAMPLE$x;
  ./goose-seq2fasta -n "Substitution$x" < SAMPLE$x > SAMPLE$x.fa 
  cat SAMPLE$x.fa >> DB.mfa;
  done
fi
###############################################################################
# SHUFFLE READS ===============================================================
if [[ "$SHUFFLE" -eq "1" ]]; then
. ShufFASTQReads.sh SAMPLE.fq > SAMPLE-SHUF.fq
mv SAMPLE-SHUF.fq SAMPLE.fq
./goose-fastq2fasta < SAMPLE.fq > SAMPLE.fa
fi
###############################################################################
# RUN FALCON ==================================================================
if [[ "$FALCON" -eq "1" ]]; then
./FALCON -v -F $FPARAM -n 4 -t $MLIMIT -x TOP-SUBS SAMPLE.fq DB.mfa
#./FALCON -v -F $FPARAM -n 4 -t $MLIMIT -x TOP-SUBS SAMPLE0.fa DB.mfa
fi
###############################################################################
# RUN MUMMER ==================================================================
if [[ "$MUMMER" -eq "1" ]]; then
rm -f TOP-MUMMER;
for((x=0 ; x<$MLIMIT ; ++x));
  do
  printf "%u\t" "$x" >> TOP-MUMMER;
  ./nucmer -p mummer-tmp SAMPLE.fa SAMPLE$x.fa
  #./nucmer -p mummer-tmp SAMPLE0.fa SAMPLE$x.fa
  ./delta-filter -1 mummer-tmp.delta > mummer-tmp.delta2
  ./show-coords -clr mummer-tmp.delta2 > mummer-tmp.delta3
  echo "Running Global similarity for MUMmer ...";
  . GlobalMUMmer.sh mummer-tmp.delta3 >> TOP-MUMMER;
  done
fi
###############################################################################
# RUN MUMMER 20 ===============================================================
if [[ "$MUMMER20" -eq "1" ]]; then
rm -f TOP-MUMMER20;
for((x=0 ; x<$MLIMIT ; ++x));
  do
  printf "%u\t" "$x" >> TOP-MUMMER20;
  ./nucmer -c 20 -p mummer-tmp SAMPLE.fa SAMPLE$x.fa
#  ./nucmer -c 20 -p mummer-tmp SAMPLE0.fa SAMPLE$x.fa
  ./delta-filter -1 mummer-tmp.delta > mummer-tmp.delta2
  ./show-coords -clr mummer-tmp.delta2 > mummer-tmp.delta3
  echo "Running Global similarity for MUMmer -c 20 ...";
  . GlobalMUMmer.sh mummer-tmp.delta3 >> TOP-MUMMER20;
  done
fi
###############################################################################
# FILTER ======================================================================
if [[ "$FILTER" -eq "1" ]]; then
cat TOP-SUBS | awk '{ print $4"\t"$3;}' | sed 's/\Substitution//g' | sort -n \
| awk '{ print $1"\t"$2;}' > TOP-SUBS-FILT;
fi
###############################################################################
# PLOT ========================================================================
if [[ "$PLOT" -eq "1" ]]; then
gnuplot << EOF
set terminal pdfcairo enhanced color
set output "mut.pdf"
set auto
set key right top
set yrange [0:100] 
set grid
#unset key
set ylabel "Similarity"
set xlabel "Mutation rate"
plot "TOP-SUBS-FILT" u 1:2 w lines title "FALCON", \
 "TOP-MUMMER" u 1:2 w lines title "MUMmer", \
 "TOP-MUMMER20" u 1:2 w lines title "MUMmer -c 20"
EOF
fi
#==============================================================================

