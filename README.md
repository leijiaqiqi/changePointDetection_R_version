# changepointtool

`changepointtool` is an R package for fixed-segment change-point detection using
a mixed-integer optimization model solved by Gurobi.

## Package structure

The exported function is:

```r
change_point_detection_fixed_num(Length, X, Delta, output_flag = 1L)
```

- `Length`: number of segments. The current formulation requires `Length >= 2`.
- `X`: numeric sequence.
- `Delta`: minimum number of observations assigned to each segment.
- `output_flag`: set to `1L` to display the Gurobi log or `0L` to suppress it.

The returned vector contains the cumulative ending positions of the segments.

## Prerequisites

You need:

1. R 3.6.0 or newer.
2. The `Matrix` package.
3. Gurobi Optimizer and a valid Gurobi license.
4. The Gurobi R package.

The Gurobi R package is installed from the `R` directory inside your Gurobi
installation. For example, in R:

```r
install.packages("/path/to/gurobi_R_package.tar.gz", repos = NULL)
```

Install `Matrix` if it is not already available:

```r
install.packages("Matrix")
```

## Install from GitHub

Install `remotes` if needed:

```r
install.packages("remotes")
```

Then install this package:

```r
remotes::install_github("leijiaqiqi/changepointtool")
```

Replace `YOUR_GITHUB_USERNAME` with the GitHub username that hosts this
repository.

## Local installation

From the directory above this package folder:

```bash
R CMD INSTALL changepointtool
```

## Example

```r
library(changepointtool)

x <- c(1, 1, 2, 2, 10, 10)
result <- change_point_detection_fixed_num(
  Length = 2,
  X = x,
  Delta = 2,
  output_flag = 1L
)

print(result)
```

## Notes

This package requires Gurobi. Installing the R package does not install the
Gurobi Optimizer or create a license automatically.
