all:
	#Makes the tprtest executable
	g++ -I./include -DTPRDEBUG main.C -o tprtest
test:
	#Check to see if we read things in more or less correctly.
	for tpr in `ls testtprs/*tpr` ; do \
		./tprtest $$tpr ; \
	done
