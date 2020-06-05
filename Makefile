CC = g++
CFLAGS = -g -Wall -std=c++11

.PHONY: debug

integral2d_test: saye_utils.o saye_algorithm.o integral2d_test.o
	g++ -g -std=c++11 integral2d_test.o saye_utils.o saye_algorithm.o -o integral2d_test.x

integral3d_test: saye_utils.o saye_algorithm.o integral3d_test.o
	g++ -g -std=c++11 integral3d_test.o saye_utils.o saye_algorithm.o -o integral3d_test.x

QMOF_test: QMOF_test.o saye_utils.o saye_algorithm.o grad_descent.o
	g++ -g -std=c++11 QMOF_test.o saye_utils.o saye_algorithm.o grad_descent.o -o QMOF_test.x

integral2d_test.o: integral2d_test.cpp
	$(CC) $(CFLAGS) -c integral2d_test.cpp

integral3d_test.o: integral3d_test.cpp
	$(CC) $(CFLAGS) -c integral3d_test.cpp

QMOF_test.o: QMOF_test.cpp
	$(CC) $(CFLAGS) -c QMOF_test.cpp

grad_descent.o: grad_descent.cpp grad_descent.h
	$(CC) $(CFLAGS) -c grad_descent.cpp

saye_utils.o: saye_utils.cpp saye_utils.h
	$(CC) $(CFLAGS) -c saye_utils.cpp

saye_algorithm.o: saye_algorithm.cpp saye_algorithm.h
	$(CC) $(CFLAGS) -c saye_algorithm.cpp
