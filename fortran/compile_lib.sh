f95 -c -fPIC *.f90
ar rcv libmy_lib.a *.o
ranlib libmy_lib.a 
rm *.o
