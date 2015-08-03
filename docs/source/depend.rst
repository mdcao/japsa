------------
Dependencies
------------

Japsa has the following dependencies, which are included in the package.

* Jfreechart >= 1.0.19 (http://www.jfree.org/jfreechart/) and Jcommon (http://www.jfree.org/jcommon/)
* colloquial.jar, now java-arithcode (https://github.com/bob-carpenter/java-arithcode)
* common-math.jar >=3.3.0 (http://commons.apache.org/proper/commons-math/)
* htsjdk >=1.126
* sam >=1.84 This will be replaced by htsjdk above
* guava >=16.0
* jhdf5 >= 18
* JRI -- Dependency on JRI will be removed in the future

Naturally, it also requires a Java Runtime Environment >=1.6 (java) installed
to run the package.

Some tools in the package have additional dependencies (not included in the package), as follows:

* *jsa.np.f5reader* (npReader): requires HDF-Java library (https://www.hdfgroup.org/products/java/release/download.html)
* *jsa.np.speciesTyping* (realtime species typing): R, RJava and jri (these dependency will be removed in a future release) and kalign2

Japsa is provided with a ready to run package at every stable release.
These pre-built releases are compiled with javac 1.6 to ensure compatibility.
If you wish to use the latest version of Japsa, or use japsa compiled with your
version of Java runtime, you will need to build from source, which requires:

* Java Development Kit (javac) >= 1.6
* make
* git/wget (optional)


Note that a specific tool may require extra dependencies (such as *bwa* etc).
Check the documentation for indivisial tools for more detailed information.
