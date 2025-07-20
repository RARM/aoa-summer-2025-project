# The Closest Pair Problem
_A Comparative Analysis of the Brute-Force and Divide-and-Conquer Approach_

## Getting Started
To get a local copy up and running, follow these simple steps.

### Prerequisites
- A C compiles (e.g., `gcc`).
- The `make` utility.

### Installation
1. Clone the repo.
```
git clone https://github.com/RARM/aoa-summer-2025-project.git
```

2. Navigate to the project directory.
```
cd aoa-summer-2025-project
```

3. Compile the project using the provided Makefile.
```
make
```

### Usage
Run the program to start the empirical analysis and get time measurements.
```
make run
```

Keep these considerations in mind:
- The repository comes with points stores in the `build/bin/data/` directory.
  The program uses those for replicability, but it can still produce different
  results due to hardware differences, OS, and system load. If you delete the
  points from the data directory, the program will generate new points and
  store them.
- The program will produce messages in the console, but it will save a CSV with
  all the measurements in the `build/bin/performance.csv` file.