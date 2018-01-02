all:
	gcc Poisson.cpp -std=c++11 -lstdc++ -o Poisson.out
run:
	./Poisson.out
draw:
	python drawPoisson.py Poisson.txt
draw-white:
	python drawPoisson.py Poisson.txt WHITE
clean:
	rm -rf Poisson.out Poisson.txt Poisson.bmp
