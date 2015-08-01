Installation
------------

There are two methods to install Japsa in your computer:


1. Install from the pre-compiled package
########################################

Pre-compiled package of Japsa is made available under each release. Installation
from this will not require extra build tools such as javac, git, and make.

Just run the install.sh script (install.bat for Windows) in the release directory::
 
   tar zxvf JapsaRelease.tar.gz
   cd JapsaRelease
   ./install.sh

The installation will ask for specific details to install the package. If you
agree with its suggestion, just type Enter. The questions are:

* *Directory to install japsa:* Enter a directory to install japsa

* *Default memory allocated to jvm:* Enter a default amount of memory allocated
  to the Java Virtual Machine. This value should be smaller than the size of
  your computer. This value, however, can be changed for each specific invocation
  of a program.

* *Enforce your jvm to run on server mode:* Type *y* if your java support running
  in server mode.

* *Path to HDF library:* Enter path to HDF library. Generally, you need to have
  HDFViewer (https://www.hdfgroup.org/products/java/release/download.html)
  installed, and enter the path to file *libjhdf5.so* (on Linux/Unix/Mac) or
  to *jhdf5.dll* (Windows). This is only required if you intend to use npReader(
  jsa.np.f5reader).


* *Path to JRI library:* Enter path to JRI library. This is required only for
  running real-time species typing (jsa.np.speciesTyping).

2. Obtain source code and compile
#################################

This method is recommended as japsa will be compiled and with the same Java
version used to run. However, Java Development Kit and Make, and preferably
git are required. This method has not been tested with Windows.

First download the latest source code::

   git clone https://github.com/mdcao/japsa
   cd japsa

or download a release::
   wget    https://github.com/mdcao/japsa/releases/download/v15.8a/JapsaRelease.tar.gz
   tar zxvf JapsaRelease
   cd JapsaRelease
   
and compile and install::      
   make install [INSTALL_DIR=~/.usr/local] [MXMEM=7000m] [SERVER=true] [JLP=/usr/lib/jni:/usr/lib/R/site-library/rJava/jri]

This will install japsa according the directives:

* *INSTALL_DIR*: specifies the directory to install japsa
* *MXMEM*: specifies the default memory allocated to the java virtual machine
* *SERVER*: specifies whether to launch the java virtual machine in server mode
* *JLP*: specifies paths to *libjhdf5*  (needed for npReader) and *libjri.so*
  (needed for real-time species typing).

If any of the above directives are not specified, the installation will ask
during the installation.

To uninstall Japsa, run the following in the japsa directory::

   make uninstall INSTALL_DIR=~/.usr/local
   
 where INSTALL_DIR points the directory Japsa was installed.

