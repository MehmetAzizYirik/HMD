# HMD

Copyright 2018 Mehmet Aziz Yirik

## Introduction

HMD (Hyperfast Molecule Deriver) is a molecular structure generator. The project started as the re-implementation of Faulon's study [1].The project is still in progress.

## Method

The method relies on the grouping of the atoms. To avoid duplicates, the equivalence classes of the atoms are calculated. For that, Canon class of CDK is used. For details about the equivalence class calculation (2). By taking the molecular information of the input atom container, first, the initial atomc ontainer is built. The molecular information is a String including atom and implicit hydrogen information. For instance, "C3C3C3" means the initial atom container includes three carbon atoms to all of which three implicit hydrogens are attached. 

For the exetnsion of an atom container from a chosen atom, at each intermediate extension, first, the equivalence classes are calculated. Then, by taking one element from each classes, the bonds are added between the atoms. This extension step (atomext) is recursively performed in atomsat function. The function saturates a chosen atom. For the saturation of all the atoms of the initial atom container, genall function is performed. The short description of main functions are given below: 
```
genall()  : It is the main structure generator function saturating all the atoms of the molecule and considering all the possible extensions.
atomsat() : This function extends the molecule until the chosen index, atom, is saturated.
atomext() : This function extends the atom in the atomcontainer by adding a new bond ( or increasing bond order) between the atom and the others.
```

## Download Source Code

It is assumed that users have git on their system and have initialised their local directory. For more information [set-up-git](https://help.github.com/articles/set-up-git/ )

To download HMD source code:

```
$ git clone https://github.com/MehmetAzizYirik/HMD.git
```
## Compiling

To compile HMD, Apache Maven and Java 1.8 (or later) are required.
```
HMD/$ mvn package
```
This command will create jar file named specifically as "jar-with-dependencies" under target folder.

## Usage

HMD.jar can be run from command line with the specified arguments. An example command is given below.

```
java -jar HMD.jar -i C3C3C2C2C1C1 -v -d C:\Users\UserName\Desktop\hmd
```

The definitions of the arguments are given below:

```
usage: java -jar maygen.jar -i <arg> [-v] -d <arg>

Generates structures for a given molecular information. The input is the
string of atom symbols with their number of implicit hydrogen.For example
'C3C3C3' means three carbon atoms each of which has three implicit
hydrogens.Besides this molecular information, the directory is needed to
be specified for the outputfile.

 -i,--molecularinfo <arg>   String of atoms with their implicit hydrogen
                            information (required)
 -v,--verbose               Print messages about the duration time of the
                            generator
 -d,--filedir <arg>         Creates and store the output sdf file in the
                            directory (required)

Please report issues at https://github.com/MehmetAzizYirik/HMD
```

## Running the Tests

For the Generator class, a test class called Test-Generator is built. This test class includes the tests of the main functions. The outputs of the the functions are tested based on the size ( or the length) of the expected output files. 

## License
This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/MehmetAzizYirik/HMD/blob/master/LICENSE) file for details

## Authors

 - Mehmet Aziz Yirik - [MehmetAzizYirik](https://github.com/MehmetAzizYirik)
 
## Acknowledgements
![YourKit](https://camo.githubusercontent.com/97fa03cac759a772255b93c64ab1c9f76a103681/68747470733a2f2f7777772e796f75726b69742e636f6d2f696d616765732f796b6c6f676f2e706e67)

The developer uses YourKit to profile and optimise code.

YourKit supports open source projects with its full-featured Java Profiler. YourKit, LLC is the creator of YourKit Java Profiler and YourKit .NET Profiler, innovative and intelligent tools for profiling Java and .NET applications.

![cdk](https://github.com/MehmetAzizYirik/HMD/blob/master/cdk.png)

This project relies on the Chemistry Development Project (CDK), hosted under [CDK GitHub](http://cdk.github.io/). Please refer to these pages for updated information and the latest version of the CDK. CDK's API documentation is available though our [Github site](http://cdk.github.io/cdk/).

## References

1- Faulon, J.L., 1992. On using graph-equivalent classes for the structure elucidation of large molecules. Journal of chemical information and computer sciences, 32(4), pp.338-348.

2- Weininger, D., Weininger, A. and Weininger, J.L., 1989. SMILES. 2. Algorithm for generation of unique SMILES notation. Journal of Chemical Information and Computer Sciences, 29(2), pp.97-101.


