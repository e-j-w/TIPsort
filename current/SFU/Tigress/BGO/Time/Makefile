all:
	echo Making BGO_TRaw.
	(cd TRaw && make -is clean && make && make -is clean && cd ..)
	echo DONE: BGO_ERAW.
	
	echo Making BGO_TCal.
	(cd TCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: BGO_ECAL.
	
	echo Making BGO_TRawSum.
	(cd TRawSum && make -is clean && make && make -is clean && cd ..)
	echo DONE: BGO_ERAWSUM.
	
	echo Making BGO_TvsRFPhase.
	(cd TvsRFPhase && make -is clean && make && make -is clean && cd ..)
	echo DONE: BGO_TvsRFPhase.

clean:
	(cd TRaw    && make -is clean && cd ..)
	(cd TCal    && make -is clean && cd ..)
	(cd TRawSum && make -is clean && cd ..)
	(cd TvsRFPhase && make -is clean && cd ..)

	rm -f *~
