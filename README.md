# A general approach of extending Python with Fortran

## Required: 
* Python 3.7 (Recommended)
* ctypes 
* numpy
* matplotlib
* A Fortran compiler supporting Fortran2003 standard

## Instructions:

* Make the shell script "compile" exectuable

```
chmod 755 compile
```

* Open the shell script "compile"; Uncomment/comment compiler options

* The default is Intel Fortran

* Run the "compile"

```
./compile

```
* Run Python script 

```
python calc_pair_correlation.py 

```
