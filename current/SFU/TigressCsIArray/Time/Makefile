all:
	echo making TigressCsIArray_TTCal.
	(cd TTCal && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressCsIArray_TTCal.
	
	echo making TigressCsIArray_TTCalFirstHit.
	(cd TTCalFirstHit && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressCsIArray_TTCalFirstHit.
	
	echo making TigressCsIArray_TDiff.
	(cd TDiff && make -is clean && make && make -is clean && cd ..)
	echo DONE: TigressCsIArray_TDiff.

clean:
	(cd TTCal && make -is clean && cd ..)
	rm -f *~
	
	(cd TTCalFirstHit && make -is clean && cd ..)
	rm -f *~
	
	(cd TDiff && make -is clean && cd ..)
	rm -f *~
