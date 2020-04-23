all:
	echo Making PINArray_TECal.
	(cd TECal && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_TECAL.
	
	echo Making PINArray_TECalSum.
	(cd TECalSum && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_TECALSUM.
	
	echo Making PINArray_ECalTFitRaw.
	(cd ECalTFitRaw && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_ECalTFitRaw.
	
	echo Making PINArray_ECalTFitRaw_Sum.
	(cd ECalTFitRaw_Sum && make -is clean && make && make -is clean && cd ..)
	echo DONE: PINARRAY_ECalFitRaw_Sum.
clean:
	(cd TECal       && make -is clean && cd ..)
	(cd TECalSum    && make -is clean && cd ..)
	(cd ECalTFitRaw && make -is clean && cd ..)
	(cd ECalTFitRaw_Sum && make -is clean && cd ..)
	rm -f *~
