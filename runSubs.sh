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
FILTER=1;
PLOT=1;
#==============================================================================
MLIMIT=31;
FQLINE=200;
FQNREADS=10000;
#==============================================================================
DISTRIBUTION="0.3,0.2,0.2,0.3,0.001";
EXTRAMUT=" -ir 0.01 -dr 0.01 ";
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
cd ../../
cp goose/src/goose-* .
cp goose/scripts/ShufFASTQReads.sh .
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
fi
###############################################################################
# SIMULATE ====================================================================
if [[ "$SIMULATE" -eq "1" ]]; then
./XS -v -ls $FQLINE -n $FQNREADS -f $DISTRIBUTION -s 0 SAMPLE.fq
# SHUF
###############################################################################
# MUTATE ======================================================================
./goose-fastq2fasta < SAMPLE.fq > SAMPLE.fa
./goose-fasta2seq   < SAMPLE.fa > SAMPLE
echo "\n\n" > SPACE;
./goose-seq2fasta -n "Substitution0" < SAMPLE > SAMPLE0.fa
cat SAMPLE0.fa SPACE > DB.mfa;
for((x=1 ; x< $MLIMIT ; ++x));
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
if [[ "$SHUFFLE" -eq "1" ]]; then
. ShufFASTQReads.sh SAMPLE.fq > SAMPLE-SHUF.fq
mv SAMPLE-SHUF.fq SAMPLE.fq
fi
###############################################################################
# RUN FALCON ==================================================================
if [[ "$FALCON" -eq "1" ]]; then
./FALCON -v -F $FPARAM -n 4 -t $MLIMIT -x TOP-SUBS SAMPLE.fq DB.mfa
fi
###############################################################################
# RUN MUMMER
if [[ "$MUMMER" -eq "1" ]]; then
# ./nucmer -maxmatch -c 30 -p test SAMPLE.fa SAMPLE1.fa
# ./show-coords -clr test.delta | awk '{print $10;'}  | tail -n 1
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
set grid
unset key
set ylabel "Similarity"
set xlabel "Mutation rate"
plot "TOP-SUBS-FILT" u 1:2 w lines
EOF
fi
#==============================================================================

