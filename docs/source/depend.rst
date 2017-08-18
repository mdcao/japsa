============
Dependencies
============

Japsa has the following dependencies, which are included in the package.

* Jfreechart >= 1.0.19 (http://www.jfree.org/jfreechart/) and Jcommon (http://www.jfree.org/jcommon/)
* colloquial.jar, now java-arithcode (https://github.com/bob-carpenter/java-arithcode)
* common-math.jar >=3.3.0 (http://commons.apache.org/proper/commons-math/)
* htsjdk >=1.126
* guava >=16.0
* jhdf5 >= 18

Naturally, it also requires a Java Runtime Environment >=1.8 (java) installed
to compile and run the package.

Some tools in the package have additional dependencies (not included in the package), as follows:

* *jsa.np.npreader* (npReader): requires HDF-Java library (https://www.hdfgroup.org/products/java/release/download.html). 

Japsa is provided with a ready to run package at every stable release.
These pre-built releases are compiled with javac 1.8 to ensure compatibility.
If you wish to use the latest version of Japsa, or use japsa compiled with your
version of Java runtime, you will need to build from source, which requires:

* Java Development Kit (javac) >= 1.8
* make
* git/wget (optional)
* mvn (optional if you want to build with maven)

Note that a specific tool may require extra dependencies (such as *bwa* etc).
Check the documentation for indivisial tools for more detailed information.
