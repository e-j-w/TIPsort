all:
	echo making Tigress_TvsRFPhase.
	(cd TvsRFPhase && make -is clean && make && make -is clean && cd ..)
	echo DONE: Tigress_TVSRFPHASE.
	
	echo making Tigress_TCal.
	(cd TCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: Tigress_TCal.
	
	echo making Tigress_TCalSum.
	(cd TCalSum && make -is clean && make && make -is clean && cd ..)
	echo DONE: Tigress_TCalSum.
	
	echo making Tigress_TTCal.
	(cd TTCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: Tigress_TTCal.
	
	echo making Tigress_TCalDiffABSuppSum.
	(cd TCalDiffABSuppSum && make -is clean && make && make -is clean && cd ..)
	echo DONE: Tigress_TCalDiffABSuppSum.
	
	echo making Tigress_TCalDiffABSupp.
	(cd TCalDiffABSupp && make -is clean && make && make -is clean && cd ..)
	echo DONE: Tigress_TCalDiffABSupp.
	
	echo making Tigress_TRaw.
	(cd TRaw && make -is clean && make && make -is clean && cd ..)
	echo DONE: Tigress_TRaw.
	
	echo making Tigress_TRaw_v_TCal.
	(cd TRaw_v_TCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: Tigress_TRaw_v_TCal.

clean:
	(cd TvsRFPhase && make -is clean && cd ..)
	(cd TCal && make -is clean && cd ..)
	(cd TTCal && make -is clean && cd ..)
	(cd TCalSum && make -is clean && cd ..)
	(cd TTCal && make -is clean && cd ..)
	(cd TCalDiffABSuppSum && make -is clean && cd ..)
	(cd TCalDiffABSupp && make -is clean && cd ..)
	(cd TRaw && make -is clean && cd ..)
	(cd TRaw_v_TCal && make -is clean && cd ..)


	rm -f *~
