all:
	echo Making FitEVsECal.
	(cd FitEVsECal && make -is clean && make && make -is clean && cd ..)
	echo DONE: FitEVsECal.
	
	echo Making PIDVsECal.
	(cd PIDVsECal && make -is clean && make && make -is clean && cd ..)
	echo DONE: PIDVsECal.

clean:
	(cd FitEVsECal           && make -is clean && cd ..)
	(cd PIDVsECal           && make -is clean && cd ..)

	rm -f *~
