all:
	#Makes the tprtest executable
	g++ -I./include -DTPRDEBUG -DTPRTEST tprplugin.C -o tprtest
pythontest:
	python3 comparison.py
test:
	#Check to see if we read things in more or less correctly.
	for tpr in `ls testtprs/*tpr` ; do \
		./tprtest $$tpr ; \
	done
