all:
	echo Making display_RawEventAboveTimestampDelta.
	(cd RawEventAboveTimestampDelta && make -is clean && make && make -is clean && cd ..)
	echo DONE: DISPLAY_RawEventAboveTimestampDelta.
	
	echo Making display_RawEvent.
	(cd RawEvent && make -is clean && make && make -is clean && cd ..)
	echo DONE: DISPLAY_RawEvent.
	
	echo Making display_RawEventForTigressCore.
	(cd RawEventForTigressCore && make -is clean && make && make -is clean && cd ..)
	echo DONE: DISPLAY_RawEventForTigressCore.
	
	echo Making display_Tigress.
	(cd Tigress && make -is clean && make && make -is clean && cd ..)
	echo DONE: DISPLAY_TIGRESS.
	
	echo Making display_CsIArray.
	(cd CsIArray && make -is clean && make && make -is clean && cd ..)
	echo DONE: DISPLAY_CsIArray.
	
	echo Making display_RF.
	(cd RF && make -is clean && make && make -is clean && cd ..)
	echo DONE: DISPLAY_RF.

clean:
	(cd RawEvent && make -is clean && cd ..)
	(cd RawEventAboveTimestampDelta && make -is clean && cd ..)
	(cd RawEventForTigressCore && make -is clean && cd ..)
	(cd Tigress  && make -is clean && cd ..)
	(cd CsIArray && make -is clean && cd ..)
	(cd RF       && make -is clean && cd ..)
	rm -f *~
