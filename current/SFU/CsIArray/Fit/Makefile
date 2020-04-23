all:
	echo Making PID.
	(cd PID && make -is clean && make && make -is clean && cd ..)
	echo DONE: PID.
	
	echo Making GoodnessOfFit.
	(cd GoodnessOfFit && make -is clean && make && make -is clean && cd ..)
	echo DONE: GoodnessOfFit.
	
	echo Making Energy.
	(cd Energy && make -is clean && make && make -is clean && cd ..)
	echo DONE: Energy.

clean:
	(cd PID           && make -is clean && cd ..)
	(cd GoodnessOfFit && make -is clean && cd ..)
	(cd Energy	  && make -is clean && cd ..)

	rm -f *~
