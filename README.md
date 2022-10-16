# Normalized Spectral Clustering
Implementation of a version of the normalized spectral clustering algorithm

# Building and Running
For the project to run you need to build the CPython extension.  
In Terminal run the following command:  
`python3 setup.py build_ext --inplace`

# Compilation
In Terminal run the following command:  
`bash comp.sh`

# Testing Algorithm
In Terminal run the following command:  
`bash tester.sh testfiles <interface> <regular> <leaks> <efficiency>`

**Syntax**:  
- interface [c|py|both]: Specifies the interface you want to test.
- regular [yes|no]:	Specifies whether to run the regular tests or not.
- leaks	[yes|no]: Specifies if you want memory leak tests.
- efficiency [yes|no]: Specifies whether you want an efficiency test on the interfaces you chose.

**Running Example**:  
`bash tester.sh testfiles both yes yes no`

