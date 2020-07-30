### RAxML

#https://github.com/stamatak/standard-RAxML.git
#go to wd and type in the terminal:
#cd standard-RAxML/
#ls Makefile.* #let's see the makefiles
#sudo make -f Makefile.AVX.mac #this will generate a RAxML executable called: raxmlHPC-AVX
#Now, to also compile the PThreads version:
#rm *.o
#and then type:
#sudo make -f Makefile.PTHREADS.mac #which will generate an executable called: raxmlHPC-PTHREADS
###################

# When to use which Version?
# The use of the sequential version is intended for small to medium datasets and for initial experiments 
# to determine appropriate search parameters.
# The PThreads version will work well for very long alignments, but performance is extremely hardware-dependent!

#make the tree with RAxML
system("~/Desktop/standard-RAxML/raxmlHPC-PTHREADS -T8 -p 12345 -s ~/Desktop/SHCS/Input/OLD/trimmed_A.fas -m GTRCAT -n TreeA_RAxMLPThreads.tre")
