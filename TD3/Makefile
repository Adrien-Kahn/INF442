kmeans: kmeans.cpp
	g++ -std=c++11 kmeans.cpp -o kmeans `pkg-config gtkmm-3.0 --cflags --libs`

itest: itest.cpp
	g++ -std=c++11 itest.cpp -o itest

clean:
	rm kmeans

.PHONY: clean
