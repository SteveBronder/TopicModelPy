f95 -c -fPIC *.f90
ar rcv libmy_lib2.a *.o
ranlib libmy_lib2.a 
rm *.o
