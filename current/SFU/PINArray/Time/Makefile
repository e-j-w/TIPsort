all:
	echo Making PINArray_TRaw.
	(cd TRaw && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_TRAW.
	
	echo Making PINArray_TRawSum.
	(cd TRawSum && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_TRAWSUM.
	
	echo Making PINArray_TCal.
	(cd TCal && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_TCAL.
	
	echo Making PINArray_TCalSum.
	(cd TCalSum && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_TCALSUM.
	
	echo Making PINArray_TCalLimSum.
	(cd TCalLimSum && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_TCALLIMSUM.
	
#	TFitRaw is causing problems compiling. I don't need it for now. So i'm not compiling it.
#	echo Making PINArray_TFitRaw.
#	(cd TFitRaw && make -is clean && make  && make -is clean && cd ..)
#	echo DONE: PINARRAY_TFITRAW.
	
	echo Making PINArray_RFUnwrapping.
	(cd RFUnwrapping && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_RFUNWRAPPING.

clean:
	(cd TRaw         && make -is clean && cd ..)
	(cd TCal         && make -is clean && cd ..)
	(cd TRawSum      && make -is clean && cd ..)
	(cd TCalSum      && make -is clean && cd ..)
	(cd TCalLimSum   && make -is clean && cd ..)
	(cd TFitRaw      && make -is clean && cd ..)
	(cd RFUnwrapping && make -is clean && cd ..)
	rm -f *~
