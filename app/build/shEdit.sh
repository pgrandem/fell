### shEdit.sh
### ----------------------------------------------------------------------------
### set programming environnement
### Pierre Grandemange
### 11/07/2017

### rep Objects/Namespaces .h files, then sleep for 300ms
gedit --new-window $repObjects/include/*.h	$repNamespaces/include/*.h &
sleep 0.3

### rep Objects/Namespaces .cc files,then sleep for 300ms
gedit --new-window $repObjects/src/*.cc	$repNamespaces/src/*.cc &
sleep 0.3

### makefile main.cc and local functions
gedit --new-window	../source/readme.txt shEdit.sh Makefile ../source/main.cc \
										../source/src/*.cc  ../source/include/*.h &

