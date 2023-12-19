LIBS = -pthread -fopenmp
CXX = g++
CXXFLAGS = -g -O0

test_logger: tests/test_logger.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) aux_functions.cpp config_reader.cpp tests/test_logger.cpp -o tests/test_logger.exe
	./tests/test_logger.exe

test_kmers: tests/test_kmers.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) aux_functions.cpp minimizer.cpp tests/test_kmers.cpp -o tests/test_kmers.exe
	./tests/test_kmers.exe 

test_file_parsing: tests/test_file_parsing.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) tests/test_file_parsing.cpp global_copan.cpp aux_functions.cpp pe_ext_reader.cpp config_reader.cpp ext_contig.cpp -o tests/test_file_parsing.exe
	./tests/test_file_parsing.exe tests/test_config.ini

test_alignments: tests/test_alignments.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) tests/test_alignments.cpp global_copan.cpp aux_functions.cpp pe_ext_reader.cpp config_reader.cpp ext_contig.cpp minimizer.cpp -o tests/test_alignments.exe
	./tests/test_alignments.exe tests/test_config.ini

test_minimizer: tests/test_minimizer.cpp
	$(CXX) $(LIBS) $(CXXFLAGS) tests/test_minimizer.cpp global_copan.cpp aux_functions.cpp pe_ext_reader.cpp config_reader.cpp ext_contig.cpp minimizer.cpp -o tests/test_minimizer.exe
	./tests/test_minimizer.exe tests/test_config.ini

clean: 
	rm -f tests/*.exe