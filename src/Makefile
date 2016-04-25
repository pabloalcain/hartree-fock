.PHONY: default library executable all clean
FC = gfortran
F2PY = f2py

default: all

all: library executable

library: johnson.so

executable: johnson.e

johnson.so: frontend.for johnson.for
	f2py -c $^ -m johnson
	cp johnson.so ..

johnson.e: johnson.for
	gfortran $^ -o $@

clean:
	rm -rfv *.so *.e ../johnson.so
