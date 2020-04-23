all:
	echo Making Energy.
	(cd Energy && make -is clean && make && make -is clean && cd ..)
	echo DONE: BGO_ENERGY.
	
	echo Making Time.
	(cd Time && make -is clean && make && make -is clean && cd ..)
	echo DONE: BGO_TIME.

clean:
	(cd Energy && make -is clean && cd ..)
	(cd Time   && make -is clean && cd ..)
	rm -f *~
