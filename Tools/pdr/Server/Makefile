

all: install

install: 
	mkdir -p /usr/local/pdr/bin
	cp *.py /usr/local/pdr/bin
	chmod u+x /usr/local/pdr/bin/PDR_SERVER.py

test: install
	cd /usr/local/pdr;python bin/PDR_SERVER.py
	
