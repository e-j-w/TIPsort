all:
	echo Making Time.
	(cd Time && make -is clean && make && make -is clean && cd ..)
	echo DONE: TIGRESSCSIARRAY_TIME.
	
	echo Making Energy.
	(cd Energy && make -is clean && make && make -is clean && cd ..)
	echo DONE: TIGRESSCSIARRAY_Energy.

clean:
	(cd Time          && make -is clean && cd ..)
	(cd Energy          && make -is clean && cd ..)
	rm -f *~
