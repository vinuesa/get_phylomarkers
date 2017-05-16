# compiler 
CXX        = g++
CXXFLAGS   = -O3 

# targets
default : ppma_2_bmp

ppma_2_bmp: 
	$(CXX) $(CXXFLAGS)  -o ppma_2_bmp ppma_io.C bmp_io.C ppma_2_bmp.C -lm

clean : FORCE
	rm -f *.o
	rm -f Phi

FORCE:



