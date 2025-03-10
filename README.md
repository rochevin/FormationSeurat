# Formation Seurat Genotoul Biostats

![](https://www.genotoul.fr/wp-content/uploads/2017/01/Biostat_logo-rvb.png)


FormationSeurat is a project focused on the formation of Seurat objects and their analysis. This project includes scripts and data files necessary for single-cell analysis using [`Seurat v5`](https://satijalab.org/seurat/).

## Project Structure

The project is organized as follows:

- `vignettes/`: Contains R scripts for basic vignettes from seurat.
  - `integration_introduction.qmd`: [Introduction to scRNA-seq integration
](https://satijalab.org/seurat/articles/integration_introduction)
  - `pbmc3k_tutorial.qmd`: [Seurat - Guided Clustering Tutorial
](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
- `data/`: Contains raw data files used in the project.
- `renv/`: Contains the `renv` lockfile and project settings (should be ignored by Git).
- `renv.lock`: Lockfile for `renv` dependencies.

## Installation

### Dependencies

To install the dependencies, follow these steps:

1. **Install `renv`**:
   ```r
   install.packages("renv")
   ```

2. **Restore the `renv` environment**:
   ```r
   renv::restore()
   ```


### Rendering the Project

To render the project using Quarto, follow these steps:

1. **Install Quarto**:
   Follow the instructions on the [Quarto website](https://quarto.org/docs/get-started/) to install Quarto on your system.

   > Quarto should be already installed in rstudio.

2. **Render the project**:
   ```sh
   quarto render 
   ```

This will generate the rendered output of your project.

By default, vignettes are not rendered. To run vignettes using quarto: 

```sh
quarto render vignettes/*
```


> Running this for the first time will download the data for the vignettes

