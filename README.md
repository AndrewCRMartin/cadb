CADB
====

makecadb/searchcadb
-------------------

### V1.0 08.10.98

### (c) 1998 UCL, Dr. Andrew C.R. Martin


`makecadb`/`searchcadb` reimplement the methodology from my D.Phil. for
creating a CA distance matrix database and searching it for loop
regions defined by distances from the N- and C-termini of the loop.

INSTALLATION
------------

To build the programs, you must have my Bioplib library installed and
a Unix DBM library (the GNU GDBM is fine).

Modify the LFLAGS and IFLAGS variables in the Makefile to point to the
Bioplib library and include files.

Build the programs by typing:
```
   make
```

CREATING A DATABASE
-------------------

A database is created from the PDB directory, by using the command:
```
   makecadb pdbdir dbfile
```   
where `pdbdir` is the PDB directory and `dbfile` is the database file 
to write. The program defaults to creating a database file with 20
distances calculated in each direction. You can change this using the
`-d` flag:
```
   makecadb -d 30 pdbdir dbfile
```
would create a database with 30 distances.

Type:
```
   makecadb -h
```
for a usage summary.

The database file will be approximately 25% bigger than the raw PDB
files from which it was calculated.

If Bioplib has been compiled with GUNZIP support enabled, then the PDB
files in the PDB directory may be gzipped.


SEARCHING THE DATABASE
----------------------

The database is searched using the searchcadb command which is driven
by a control file (or by commands typed at the prompt).

Type:
```
   searchcadb -h
```
for a usage summary and type:
```
   seachcadb
```
then:
```
   help
   quit
```
to get a summary of commands which may be used in the control file.

To use the program, your control file must specify the database and
the loop length for which you are searching. e.g.
```
   database pdb.081098.20
   length 17
```
This is then followed by the distance constraints you wish to apply.
Normally you have a set of `dp` constraints --- positive constraints
from the N-terminus and `dp` constraints --- negative constraints from
the C-terminus.

The `dp` or `dm` commands are followed by the offset to the residue for
which you are calculating a distance (e.g. 2 will be 2 residues along)
then the minimum and maximum allowed distances for this constraint.

Each of these is normally in 2 sets, those which are local to that
terminus and those which span to the other side of the loop. The
offset for the local constraints will not change with loop length. The
offset for those which span the loop will vary with the loop length.

For example you might use the following set of constraints to pull out
17 residue loops:
```
   ! +ve constraints within Nter
   dp 2 5.98 7.74
   dp 3 6.52 11.89
   ! +ve constraints to Cter
   dp 16 10.8 14.3
   dp 15 8.45 12.27
   ! -ve constraints within Cter
   dm 2 6.1 7.56
   dm 3 7.04 12.22
   ! -ve constraints to Nter
   dm 16 10.8 14.3
   dm 15 10.07 14.99
```
while you would apply the following set for 10 residue loops:
```
   ! +ve constraints within Nter
   dp 2 5.98 7.74
   dp 3 6.52 11.89
   ! +ve constraints to Cter
   dp 9 10.8 14.3
   dp 8 8.45 12.27
   ! -ve constraints within Cter
   dm 2 6.1 7.56
   dm 3 7.04 12.22
   ! -ve constraints to Nter
   dm 9 10.8 14.3
   dm 8 10.07 14.99
```
Note that the actual distances haven't changed in this case, just the
offsets. Note also that only the offsets which span the loop have
changed. 

At the end of the control file, you put the command:
```
   end
```
which will actually run the search.

See the paper: Martin et al. PNAS 86(1989),9269-9272 for details of
this method.


To run the program having written a control file, use the following:
```
   searchcadb controlfile resultsfile
```
where `controlfile` is the controlfile you have created and
`resultsfile` is your output reslts file.

