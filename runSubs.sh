#!/bin/bash
###############################################################################
# PERFORMANCE OF FALCON HANDLING INCREASING SUBSTITUTIONAL MUTATIONS RATES
# REQUIREMENTS: linux OS with cmake, git and gnuplot
# (sudo apt-get install cmake git gnuplot)
###############################################################################
# PARAMETERS ==================================================================
MLIMIT=31;
FQLINE=200;
FQNREADS=10000;
INSTALL=1;
SIMULATE=1;
FALCON=1;
FILTER=1;
PLOT=1;
DISTRIBUTION="0.3,0.2,0.2,0.3,0.001";
EXTRAMUT="" # -ir 0.01 -dr 0.01 ";
FPARAM=" -m 20:500:1:3/50 -m 14:100:1:0/0 -m 12:1:0:0/0 -m 4:1:0:0/0 \
-c 10 -g 0.95 ";
###############################################################################
if [[ "$INSTALL" -eq "1" ]]; then
# CLEAN & INSTALL =============================================================
rm -fr goose-* FALCON XS COUNT goose/ falcon/ xs/ count/ SAMPLE* TOP* DB-mfa;
# GET GOOSE FRAMEWORK =========================================================
git clone https://github.com/pratas/goose.git
cd goose/src/
make
cp goose-* ../../
cd ../../
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
# GET COUNT ===================================================================
git clone https://github.com/pratas/count.git
cd count
cmake .
make
cp COUNT ../
cd ..
fi
###############################################################################
# SIMULATE ====================================================================
if [[ "$SIMULATE" -eq "1" ]]; then
./XS -v -ls $FQLINE -n $FQNREADS -f $DISTRIBUTION -s 0 SAMPLE.fq
###############################################################################
# MUTATE ======================================================================
./goose-fastq2fasta < SAMPLE.fq > SAMPLE.fa
./goose-fasta2seq   < SAMPLE.fa > SAMPLE
rm -f DB.mfa SPACE;
echo "\n\n" > SPACE;
for((x=0 ; x< $MLIMIT ; ++x));
  do
  MRATE=`echo "scale=2;$x/100" | bc -l`;
  echo "Substitutions rate: $MRATE";
  ./goose-mutatedna -mr $MRATE $EXTRAMUT < SAMPLE > SAMPLE$x;
  ./goose-seq2fasta -n "Substitution$x" < SAMPLE$x > SAMPLE$x.fa 
  cat SAMPLE$x.fa SPACE >> DB.mfa;
  done
fi
###############################################################################
# RUN FALCON ==================================================================
if [[ "$FALCON" -eq "1" ]]; then
./FALCON -v -F $FPARAM -n 4 -t $MLIMIT -x TOP-SUBS SAMPLE.fq DB.mfa
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
set key right bottom
set yrange [0:100] 
#set xrange [0:50]  
set grid
unset key
set ylabel "Similarity"
set xlabel "Mutation rate"
plot "TOP-SUBS-FILT" u 1:2 w lines
EOF
fi
#==============================================================================
