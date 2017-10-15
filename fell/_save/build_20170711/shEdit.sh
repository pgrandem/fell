### shEdit.sh
### Pierre Grandemange
### 11/07/2017

### rep Objects/Namespaces .h files
gedit --new-window $repObjects/include/*.h	$repNamespaces/include/*.h &

### rep Objects/Namespaces .cc files
gedit --new-window $repObjects/src/*.cc	$repNamespaces/src/*.cc &

### makefile main.cc and local functions
gedit --new-window	shEdit.sh Makefile ../source/main.cc \
										../source/src/*.cc  ../source/include/*.h &

