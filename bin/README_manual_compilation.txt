# last update: 2021-09-14
# FastTree [ http://microbesonline.org/fasttree/#Install ]

cd FastTree_v2.1.11

wget -c wget -c http://www.microbesonline.org/fasttree/FastTree-2.1.11.c

# compile with double precision to resolve short branches!
gcc -static -s -DUSE_DOUBLE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree-2.1.11.c -lm
#gcc -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree.c -lm
# or if that fails
gcc -DNO_SSE -O3 -finline-functions -funroll-loops -Wall -o FastTree FastTree-2.1.11.c -lm
cd ..

# clustal-omega

cd clustal-omega-1.2.4
./configure
make
cd ..

# parallel

cd parallel-20171022
./configure
make
cd ..

# PhiPack

cd PhiPack/src
make
cd ../..

# consense

cd phylip-3.695/src/ 
make -f Makefile.unx consense pars seqboot
cd ../..

# PAUP

# Please visit https://people.sc.fsu.edu/~dswofford/paup_test/ , download the appropriate version and put in bin/

# bc

please install from your package manager (apt, yum ...)

