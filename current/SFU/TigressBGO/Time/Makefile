all:
	echo making TigressBGO_TTCalCFD.
	(cd TTCalCFD && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressBGO_TTCalCFD.
	
	echo making TigressBGO_TTCal.
	(cd TTCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressBGO_TTCal.
	
	echo making TigressBGO_TTHP.
	(cd TTHP && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressBGO_TTHP.
	
	echo making TigressBGO_TTHPCal.
	(cd TTHPCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressBGO_TTHPCAL.
	
	echo making TigressBGO_TTHPCal_Sum.
	(cd TTHPCal_Sum && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressBGO_TTHPCal_Sum.
	
	#echo making TigressBGO_TTHP_Map.
	#(cd TTHP_Map && make -is clean && make && make -is clean && cd ..)
	#echo DONE: TigressBGO_TTHP_Map.
	
	echo making TigressBGO_TTHPCal_Map.
	(cd TTHPCal_Map && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressBGO_TTHPCal_Map.
	
	echo making TigressBGO_TDiff.
	(cd TDiff && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressBGO_TDiff.
	
	echo making TigressBGO_TDiffSum.
	(cd TDiffSum && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressBGO_TDiffSum.

clean:
	(cd TTCalCFD    && make -is clean && cd ..)
	(cd TTCal       && make -is clean && cd ..)
	(cd TTHP        && make -is clean && cd ..)
	(cd TTHPCal     && make -is clean && cd ..)
	(cd TTHPCal_Sum && make -is clean && cd ..)
	#(cd TTHP_Map    && make -is clean && cd ..)
	(cd TTHPCal_Map && make -is clean && cd ..)
	(cd TDiff       && make -is clean && cd ..)
	(cd TDiffSum    && make -is clean && cd ..)


	rm -f *~
