all: sddsdemo

SDDSINC=-I$(firstword $(wildcard /usr/include/SDDS))
SDDSLIB=-L$(firstword $(wildcard /usr/lib/x86_64-linux-gnu/SDDS /usr/lib/i386-linux-gnu/SDDS /usr/lib64/SDDS /usr/lib/SDDS))

sddsdemo: sddsdemo.c
	gcc $(SDDSINC) -o sddsdemo sddsdemo.c $(SDDSLIB) -lmdbcommon -lSDDS1 -lmdblib -llzma -lz

clean:
	rm -f sddsdemo demo.sdds


















