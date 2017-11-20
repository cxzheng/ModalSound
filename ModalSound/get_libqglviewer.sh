
wget http://www.libqglviewer.com/src/libQGLViewer-2.7.1.tar.gz
tar -xzf libQGLViewer-2.7.1.tar.gz
cd libQGLViewer-2.7.1/QGLViewer
qmake
make -j
sudo make install
