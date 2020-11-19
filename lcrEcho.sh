pip install -r requirements.txt --upgrade

export PYTHONPATH=$PYTHONPATH:.

python assembly/app/lcr2/lcr2_pipeline.py \
	https://ice.synbiochem.co.uk \
	davedave@manchester.ac.uk \
	****** \
	LCR \
	True \
	SBC008019 \
	SBC008028

javac lcr5.java
java lcr5


SUB=$(echo $(date +%Y) | cut -c3-4)
SUB2=$(date +%m)
SUB3=$(date +%d)
DIRECT=$SUB$SUB2$SUB3"LCR"


mv "out/"$DIRECT"/1/worklist.csv" "out/"$DIRECT"/1/1_PCR_worklist.csv"
mv "out/"$DIRECT"/2/worklist.csv" "out/"$DIRECT"/2/2_DIGEST_worklist.csv"
mv "out/"$DIRECT"/3/worklist.csv" "out/"$DIRECT"/3/3_QC_worklist.csv"
mv "out/"$DIRECT"/4/worklist.csv" "out/"$DIRECT"/4/4_BRIDGING_OLIGOS_worklist.csv"
mv "out/"$DIRECT"/5/worklist.csv" "out/"$DIRECT"/5/5_DIGEST_POOLING_worklist.csv"
mv "out/"$DIRECT"/6/worklist.csv" "out/"$DIRECT"/6/6_LCR_worklist.csv"




mv "out/"$DIRECT"/1" "out/"$DIRECT"/1_PCR"
mv "out/"$DIRECT"/2" "out/"$DIRECT"/2_DIGEST"
mv "out/"$DIRECT"/3" "out/"$DIRECT"/3_QC"
mv "out/"$DIRECT"/4" "out/"$DIRECT"/4_BRIDGING_OLIGOS"
mv "out/"$DIRECT"/5" "out/"$DIRECT"/5_DIGEST_POOLING"
mv "out/"$DIRECT"/6" "out/"$DIRECT"/6_LCR"


