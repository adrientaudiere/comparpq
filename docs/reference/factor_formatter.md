# Format factor columns with funky colored backgrounds

[![lifecycle-experimental](https://img.shields.io/badge/lifecycle-experimental-orange)](https://adrientaudiere.github.io/MiscMetabar/articles/Rules.html#lifecycle)

## Usage

``` r
factor_formatter(x)
```

## Arguments

- x:

  A factor or character vector

## Value

A formattable formatter function that applies colored backgrounds based
on unique factor levels.

## Examples

``` r
fmt <- factor_formatter(c("low", "medium", "high"))
class(fmt)
#> [1] "formatter" "function" 
```
