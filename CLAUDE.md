# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**comparpq** is an R package extending MiscMetabar and phyloseq for comparing phyloseq objects in microbiome analysis. It focuses on taxonomic comparison, accuracy metrics computation, and interactive visualizations.

## Common Commands

You may propose to use the following commands when working with the package. Please do not run them automatically; only suggest them when appropriate.

```bash
# Run code with loaded package
Rscript -e "devtools::load_all(); code"

# Run all tests
Rscript -e "devtools::test()"

# Run tests for files starting with {name}
Rscript -e "devtools::test(filter = '^{name}')"

# Run tests for R/{name}.R
Rscript -e "devtools::test_active_file('R/{name}.R')"

# Run a single test "blah" for R/{name}.R
Rscript -e "devtools::test_active_file('R/{name}.R', desc = 'blah')"

# Generate documentation
Rscript -e "devtools::document()"

# Check pkgdown documentation
Rscript -e "pkgdown::check_pkgdown()"

# Full package check
Rscript -e "devtools::check()"
```

## Architecture

### Core Dependencies
- **phyloseq**: S4 class system for microbiome data (slots: @tax_table, @otu_table, @refseq, @sam_data)
- **MiscMetabar**: Parent package with metabarcoding utilities
- **S7**: Modern class system used for `list_phyloseq` class

### Key Modules

| Module | Purpose |
|--------|---------|
| `list_phyloseq.R` | S7 class for storing/comparing multiple phyloseq objects (6 comparison types: REPRODUCIBILITY, ROBUSTNESS, NESTED_ROBUSTNESS, REPLICABILITY, EXPLORATION, SEPARATE_ANALYSIS) |
| `compare_taxo.R` | Taxonomic accuracy metrics using confusion matrix approach (TP/FP/FN/TN, FDR, TPR, PPV, F1_score, MCC, ACC) |
| `compare_taxo_plot.R` | Visualization for taxonomic comparisons (`tc_bar`, `tc_circle`) |
| `bubbles_pq.R` | Interactive bubble plots using Observable HQ notebooks |
| `fake_creation.R` | Mock community preparation (`add_shuffle_seq_pq`, `add_external_seq_pq`) |
| `taxtab_modification.R` | Tax table utilities (rename, select, resolve conflicts) |
| `analysis_lpq.R` | Statistical analysis (PERMANOVA/ADONIS) for list_phyloseq |
| `formattable_lpq.R` | Formatted table visualizations with color bars |

## Coding Conventions

- Use base pipe (`|>`) not magrittr (`%>%`)
- Use `\() ...` for single-line anonymous functions, `function() {...}` otherwise
- Format with `air format .` after generating code
- Line length limit: 120 characters
- Tests for `R/{name}.R` go in `tests/testthat/test-{name}.R`

## Documentation Requirements

- Every user-facing function must be exported with roxygen2 documentation
- Wrap roxygen comments at 80 characters
- Parameter docs start with type: `@param x (Int, required) An integer vector of length 1`
- Add new topics to `_pkgdown.yml` reference section
- Run `pkgdown::check_pkgdown()` to verify all topics are indexed

## Test Writing Guidelines

- Use `testthat` framework
- Name test files `test-{name}.R` for `R/{name}.R`
- Each test case should have a descriptive name
- Cover edge cases and typical usage scenarios


## NEWS.md Updates

- Every user-facing change gets a bullet (no line wraps)
- Mention function name early in bullet
- Include related issue in parentheses
- Order bullets alphabetically by function name

## Commit Convention

Follow Conventional Commits: `<type>[optional scope]: <description>`
- Types: feat, fix, docs, style, refactor, test, chore
- Lowercase type and description, imperative mood
