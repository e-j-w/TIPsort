all:
	echo Making TCal.
	(cd TCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: TCAL.
	
	echo Making TvsRFPhase.
	(cd TvsRFPhase && make -is clean && make && make -is clean && cd ..)
	echo DONE: TvsRFPHASE.
	
	echo Making TTCal.
	(cd TTCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: TTCAL.
	
	echo Making TFitCFD.
	(cd TFitCFD && make -is clean && make && make -is clean && cd ..)
	echo DONE: TFITCFD.
	
	echo Making TFitCFDChannel.
	(cd TFitCFDChannel && make -is clean && make && make -is clean && cd ..)
	echo DONE: TFITCFDChannel.
	
	echo Making TRaw.
	(cd TRaw && make -is clean && make && make -is clean && cd ..)
	echo DONE: TRAW.
	
	echo Making TRaw_v_TCal.
	(cd TRaw_v_TCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: TRAW_V_TCAL.

clean:
	(cd TCal       && make -is clean && cd ..)
	(cd TvsRFPhase && make -is clean && cd ..)
	(cd TTCal      && make -is clean && cd ..)
	(cd TFitCFD    && make -is clean && cd ..)
	(cd TFitCFDChannel    && make -is clean && cd ..)
	(cd TRaw       && make -is clean && cd ..)
	(cd TRaw_v_TCal       && make -is clean && cd ..)

	rm -f *~
