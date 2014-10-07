#!/bin/bash
PWD=`pwd`
DATA=$PWD/data
RES=$PWD/results

echo $DATA
echo $RES

export CLASSPATH=$PWD:$PWD/lib/commons-cli-1.2.jar

for DATATYPE in "easy" "hard" "bio" "nlp" "finance" "speech" "vision"

do
java cs475.Classify -mode train -algorithm pegasos -model_file $RES/${DATATYPE}_pegasos.model -data $DATA/$DATATYPE.train

java cs475.Classify -mode test -model_file $RES/${DATATYPE}_pegasos.model -data $DATA/${DATATYPE}.dev -predictions_file $RES/${DATATYPE}_pegasos_predictions.txt

echo "Pegasos: $DATATYPE.dev"
python $DATA/compute_accuracy.py $DATA/${DATATYPE}.dev $RES/${DATATYPE}_pegasos_predictions.txt 
echo ""

done
