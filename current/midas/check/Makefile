all:
	echo Making check_MidasTrigger.
	(cd trigger && make -is clean && make  && make -is clean && cd ..)
	echo DONE: check_MidasTrigger
	
	echo Making check_MidasTrigFragment.
	(cd trigFragment && make -is clean && make  && make -is clean && cd ..)
	echo DONE: check_MidasTrigFragment

clean:
	(cd trigger      && make -is clean && cd ..)
	(cd trigFragment      && make -is clean && cd ..)
	rm -f *~
