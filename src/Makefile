libaminda.so : update_labeling.o image_ops.o
	g++ -shared -fPIC -fopenmp -o $@ $^

%.o : %.cpp
	g++ -c -O3 -fPIC -fopenmp -o $@ $^

.PHONY : clean
clean :
	-rm libaminda.so update_labeling.o image_ops.o
