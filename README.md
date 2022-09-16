# lmpanlys
Tools for LAMMPS. 

## Installation

Clone the repository

```
  1. git clone https://github.com/vt87/lmpanlys.git
```

  Add the following in ~/.bashrc
```  
  2. export PYTHONPATH="{PATH}/lmpanlys/src:$PYTHONPATH"
```
  PATH is where git directory is cloned. 

  Source the bashrc using the following command
```
  3. source ~/.bashrc
```


*Example*
```
cd ~/test

git clone https://github.com/vt87/lmpanlys.git

export PYTHONPATH="~/test/lmpanlys/src:$PYTHONPATH"

source ~/.bashrc
```
## Usage

### Obtaining the pair coefficients in LAMMPS format for Martini by reading GROMACS Martini database file

Create martini.inp file containing lammps type and corresponding Martini label separated by at least one space. 

*Example*

1 C5

2 TN1a

3 P1

*NOTE* : Different lines for different types.

Run the following script on Python

```
from paircoeffs_martini import paircoeffs_martini
m = paircoeffs_martini(fname = "martini.inp",
                       database = $GROMACS_DATABASE_FILE$)
```
*GROMACS_DATABASE_FILE* is the location of Martini 3.0.0 database file in your computer.

For reference, file is provided in utils/databases/martini_v3.0.0.itp

Output of the script will generate paircoeffs.txt in LAMMPS format. Example below

```
# MARTINI force field
# The bead id and its martini label are as follows
# 1 C5
# 2 TN1a
# 3 P1
pair_style lj/gromacs 9.0 12.0
pair_coeff 1 1 0.810229 4.700000 9.0 12.0
pair_coeff 1 2 0.552103 3.950000 9.0 12.0
pair_coeff 1 3 0.733748 4.700000 9.0 12.0
pair_coeff 2 2 0.360899 3.400000 9.0 12.0
pair_coeff 2 3 0.552103 3.950000 9.0 12.0
pair_coeff 3 3 0.927342 4.700000 9.0 12.0
```

*NOTE*: lj/gromacs is part of the EXTRA-PAIR package

## Additional options

Script with additonal options

```
from paircoeffs_martini import paircoeffs_martini
m = paircoeffs_martini(fname = "martini.inp",
                       database = "../utils/database_martini/martini_v3.0.0.itp",
                       cutoff = 1.2,
                       inner_cutoff = 0.9,
                       fftype = "lj/gromacs",
                       ofile = "paircoeffs.txt")
```
```fname``` : Name of input file

```database``` : Location of Gromacs Martini database file

```cutoff``` : LJ cutoff *(in nanometers)*

```inner_cutoff``` : LJ inner cutoff *(in nanometers)* . Used in ```lj/gromacs```

```fftype``` : Type of potential. Only two options are allowed; ```lj/cut``` or ```lj/gromacs```

```ofile``` : Name of output file

## Testing

Go to tests folder. Run ```martini_test.py```. 

*Bug me if you find bugs*

