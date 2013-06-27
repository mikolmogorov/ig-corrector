PY="python"

all: cpptools calign

cpptools:
	cd src; make

calign:
	cd third-party/calign; $(PY) setup.py build --build-platlib=.

clean:
	cd src; make clean
	cd third-party/calign; rm align.so; rm -rf build
