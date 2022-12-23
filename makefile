NTESTS = -DNTESTS=1000000 # Number of tests

All:framework.cpp Methods.h
	g++ framework.cpp Methods.h -DMETHOD="Baydoun" $(NTESTS) -o Baydoun -std=c++20
	g++ framework.cpp Methods.h -DMETHOD="Vieta" $(NTESTS) -o Vieta -std=c++20
Baydoun:framework.cpp Methods.h
	g++ framework.cpp Methods.h -DMETHOD="Baydoun" $(NTESTS) -o Baydoun -std=c++20
Vieta:framework.cpp Methods.h
	g++ framework.cpp Methods.h -DMETHOD="Vieta" $(NTESTS) -o Vieta -std=c++20
Clear:
	$(RM) Baydoun
	$(RM) Vieta
