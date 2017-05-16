
# FastTree [ http://microbesonline.org/fasttree/#Install ]

cd FastTree_v2.1.10
gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
# or if that fails
gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
cd ..

# clustal-omega

cd clustal-omega-1.2.4
./configure
make
cd ..

# pexec-1.0rc8

cd pexec-1.0rc8
./configure
make cd ..

# PhiPack

cd PhiPack/src
make
cd ../..

# PAUP

# Please visit https://people.sc.fsu.edu/~dswofford/paup_test/ , download the appropriate version and put in bin/
