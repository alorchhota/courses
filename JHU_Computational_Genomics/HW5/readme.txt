* Answer to Q1-Q5 are available in hw5-answer.txt file.
* The programs are written in python 2.7.
* Command to run codes for Q4:

python franken.py data/frankengene1.fasta.txt data/trainingData1.txt results/pred.txt data/testData1.txt

The general form of the command:

python franken.py GENOME_FILE_NAME TRAIN_DATA_FILE_NAME PREDICTION_FILE_NAME TEST_DATA_FILE_NAME

Here, franken.py takes 4 arguments: 
1) the file name of genome data.
2) the file name of training data.
3) the file name of the ouput file where the predicions will be saved.
4) the file name of test data.

Prediction of franken genome has been save in results/prediction.txt file.

* Command to run codes for Q5.

python 5.orf.py < INPUT_FILE > OUTPUT_FILE

Note: This program needs the codon-map file which has been saved in inputs/rna-codon.txt file. If you use different directory, then please pass the codon-map file name as 1st argument.
