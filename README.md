# Branch-and-Bound and Cutting plane algorithms for the TSP

**Author** : Théo Guyard, Nans Préjean

**Description** : Work done during my 4th school yeat at INSA Rennes in the Mathematical Departement. The aim is to solve the TSP using two different methods.

**Language** : `Python 3`

**Dependencies** : `numpy`, `networkx`, `pulp`, `tsplib95`


---

## Features

- [x] Data reader from `.tsp` files
- [x] BnB solving method
- [x] CP solving method

## Usage

Clone the repository :
```
git clone https://github.com/TheoGuyard/TSP-BnB-CP.git
```
### Option 1 :
Go into the folder :
```
cd TSP-BnB-CP
```
Run one of the following commands depending on the wanted solving method :
```
python cuttingplanes.py <dataset> <maxtime> 
python branchandbound.py <dataset> <maxtime> 
# Example : 
# $ python cuttingplanes.py pr107.tsp 600
```
where `<maxtime>` is expressed in seconds.

### Option 2 :
Open `main.py`, modify it and run it.

## Datasets

TSP datasets can be found in the [TSPLIB](http://elib.zib.de/pub/mp-testdata/tsp/tsplib/tsp/) website. Make sure that data files are placed under the `benchmarks` folder.
