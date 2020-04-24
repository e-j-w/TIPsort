all:
	echo Making PINArray_ERaw.
	(cd ERaw && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_ERAW.
	
	echo Making PINArray_ERawSum.
	(cd ERawSum && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_ERAWSUM.
	
	echo Making PINArray_ECal.
	(cd ECal && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_ECAL.
	
	echo Making PINArray_ECalSum.
	(cd ECalSum && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_ECALSUM.
	
	# echo Making PINArray_ECalSumRing.
	# (cd ECalSumRing && make -is clean && make  && make -is clean && cd ..)
	# echo DONE: PINARRAY_ECALSUMRING.
	
	echo Making PINArray_ECalLimSum.
	(cd ECalLimSum && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_ECALLIMSUM.
	
	echo Making PINArray_ECalLim.
	(cd ECalLim && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_ECALLIM.
	
	echo Making PINArray_ECalSum_TLim.
	(cd ECalSum_TLim && make -is clean && make  && make -is clean && cd ..)
	echo DONE: PINARRAY_ECALSUM_TLIM.

clean:
	(cd ERaw          && make -is clean && cd ..)
	(cd ERawSum       && make -is clean && cd ..)
	(cd ECal          && make -is clean && cd ..)
	(cd ECalSum       && make -is clean && cd ..)
	(cd ECalSumRing   && make -is clean && cd ..)
	(cd ECalLim       && make -is clean && cd ..)
	(cd ECalLimSum    && make -is clean && cd ..)
	(cd ECalSum_TLim  && make -is clean && cd ..)
	rm -f *~
