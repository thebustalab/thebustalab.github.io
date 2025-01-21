--- 
title: "Integrated Bioanalytics"
author: "Lucas Busta and members of the Busta lab"
date: "2025-01-19"
site: bookdown::bookdown_site
documentclass: krantz
bibliography: [book.bib, packages.bib]
# url: your book url like https://bookdown.org/yihui/bookdown
# cover-image: path to the social sharing image like images/cover.jpg
description: |
  Integrated Bioanalytics is a book describing how to perform chemical, phylogenetic, and genomic analyses.
biblio-style: apalike
csl: chicago-fullnote-bibliography.csl
output:
  bookdown::bs4_book:
    theme:
      primary: "#3860b6" #links  
      base_font: 
        google: 
          family: Lato
      heading_font:
        google:
          family: Montserrat
          wght: 600
      code_font:
        google: 
          family: Roboto Mono
      bg: "#fefefe" #backgrounds
      fg: "#000000" #fonts
    repo: 
      base: https://github.com/thebustalab/thebustalab.github.io/tree/master/integrated_bioanalytics
      branch: main
    includes:
      in_header: style/ga.html
    template: style/bs4_book.html
    css: style/style.css
---



# WELCOME {-}

<!-- start preface-->
<img src="https://thebustalab.github.io/integrated_bioanalytics/images/cover3.png" width="100%" style="display: block; margin: auto;" />

Integrated Bioanalytics documents methods for analyzing chemical and sequence data in R as well as some basics of scientific writing. It is maintained by Lucas Busta and members of the Busta lab. To run the analyses described in this book you will need to run a source script that will set up your R environment with a variety of packages, custom functions, and datasets. If you don't have R, see "installation" under "Data Analysis In R" in the table of contents. Run the source script by pasting and executing the following in your R command line (RStudio recommended). If you are in the Busta Lab (or want access to full features), define an object `bustalab = TRUE` before running the source command. If you have trouble running the source script, please reach out to Lucas Busta at: bust0037@d.umn.edu. The source script: 


```r
source("https://thebustalab.github.io/phylochemistry/phylochemistry.R")
```

<!-- Features provided by the source script: -->

<!-- --Analysis and visualization tools-- -->

<!-- * A GC-MS data analysis application with a MS reference library. -->
<!-- * A sequence alignment analysis application for trimming alignments. -->
<!-- * BLAST searches that export .fasta files of hits and store results in a .csv file. -->
<!-- * Functions for dimensionality reduction, clustering, modeling, and visualization. -->

<!-- --Useful data and data structures-- -->

<!-- * More than 12 chemical data sets for running practice analyses. -->
<!-- * A phylogeny for >30,000 plant species, including nearly 30,000 angiosperms, >500 gymnosperms, nearly 500 pteridophytes, and 100 bryophytes (https://thebustalab.github.io/data/plant_phylogeny.newick). -->
<!-- * A list of nearly 400,000 plant species as well as the families, orders, and phyla to which they belong (https://thebustalab.github.io/data/plant_species.csv). -->
<!-- * Support for multiple column and row names. -->

<!-- ## new features

1. A Shiny app for GC-FID and GC-MS data analysis, including a large MS library.
2. Open reading frame extraction from multiple fasta files.
3. 
4. Minor ticks for ggplot2 axes.
5. Phylogenetic signal for discrete traits.
6. Analyze multiple sequence alignments for sites associated with user-defined function
7. Multiple column name, multiple row name data structures (aka "polylists").
8. Draw annotated multiple sequence alignments.
9. Use image analysis to automatically get the csv of a mass spectrum from a published image.
10. Draw chemical structures in R from a csv of molecular coordinates.

## wrapped features

1. BLAST transcriptomes, via [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download).
2. Multiple sequence alignments and codon alignments of amino acid and nucleotide sequences, via [msa](https://bioconductor.org/packages/release/bioc/html/msa.html) and [orthologr](https://github.com/HajkD/orthologr).
3. Phylogenetic tree construction (including g-blocks trimming, pruning, ancestral states reconstruction), via [phangorn](https://cran.r-project.org/web/packages/phangorn/index.html).
4. Systematic read/write functions (csv, newick, wide tables, fasta, summary statistic tables, GFFs, chromatograms, mass spectra).
5. Phylogenetic signal for continuous traits, via [picante](https://cran.r-project.org/web/packages/picante/index.html). -->

<!-- end -->

# (PART) GETTING STARTED 

<!-- start overview -->

# overview {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/chemometrics.jpeg" width="100%" style="display: block; margin: auto;" />

In bioanalytical science, we separate, identify, and quantify matter - be it DNA, RNA, proteins, small molecules, or even atoms. To connect our data with the world around us and answer scientific questions, multiple chemical entities must be separated, quantified, and identified. As our ability to collect analytical data expands, so must our ability to effectively analyze that data - whether its 10 data points or 10,000.

This book first covers data analysis in R. We will first look at tools for hypothesis generation, including: (i) encoding variables in visual representations of data and (ii) summarizing and providing overviews of large data set. We will then turn to evaluating hypothesese with data by looking at statistical tests and models. Finally, we will look at how to communicate our results in a clear and effective way. These techniques will also allow us to answer common quesions we may have about our data: "Which of my samples are most closely related?", "Which analytes are driving differences among my samples?", "Do my samples fall into definable clusters?", "Are any of my variables related?", and "Are any of these distributions different?".

Let's get started!

<!-- end -->

<!-- start installation -->

# installation {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/installation.png" width="100%" style="display: block; margin: auto;" />

## R {-}

R is the computing language we will use to run our analyses and produce high quality plots. If you already have R installed (you will need at least version 4.1.1), you can go straight to installing RStudio. If not, follow these steps to install R:

1. Go to https://cran.r-project.org/

2. Click on "Download R for \<your operating system\>" (see footnote), depending on your operating system you will select "Download R for Linux", "Download R for (Mac) OS X", or "Download R for Windows".

We will use \<this notation\> quite a bit. It indicates a place where you should insert information, data, or something similar that corresponds to your particular situation. In this example it means insert "your operating system", i.e. Linux, (Mac) OS X, or Windows.

3. For Mac: download the .pkg file for the latest release. For PC: click "install R for the first time", then click "Download R \<version\> for Windows".

4. After the executable finishes downloading (in Windows, it is a file with .exe extension; for Mac, it is a .dmg file or a .dmg inside a .pkg file), open the file as an administrator, and follow the installation instructions. R should install without any problems. You can click OK for all of the windows that pop-up during installation, and choose a "regular" installation (if given the choice). 

If you have trouble installing R please google "Install R Mac" or "Install R PC" and follow one the many video tutorials out there. If you have tried this and are still having trouble, please contact me.

## RStudio {-}

Once we install R, we can install RStudio, which is essentially a convenient way of interacting with R. Some people do not like RStudio and prefer to interact with R directly. This is fine, but many beginning R users find RStudio helpful, so I recommend it. Follow these steps to install RStudio:

1. Go to https://rstudio.com/

2. Click "DOWNLOAD" at the top of the page.

3. Click the "DOWNLOAD" button that corresponds to RStudio Desktop with the free Open Source License.

4. The page may automatically detect which operating system you are using and recommend a version for you. If it does, download that file (.exe for PC or .dmg for Mac). If not, scroll down to the "All Installers" section and download the file that is right for you. Open the file as an administrator, and follow the installation instructions. RStudio should install without any problems. You can click OK for all of the windows that pop-up during installation, and choose a "regular" installation (if given the choice).

If you have trouble installing RStudio please google "Install RStudio Mac" or "Install RStudio PC" and following one the many video tutorials out there. If you have tried this and are still having trouble, please contact me.

## verification {-}

Open RStudio by clicking on the appropriate file in your applications folder, or wherever it is saved on your computer. If you are on Windows, be sure to run RStudio as administrator. You will see several windows. One is the Code Editor, one is the R Console, one is the Workspace and History, and one is the Plots and Files window.

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/rstudio_components.png" width="100%" style="display: block; margin: auto;" />

The R Console window should have a `>` in it. Type `head(Indometh)`. This should display the first six lines of a data set describing the pharmacokinets of indomethacin. This is one of the built in datasets in R - you do not need any additional files to run this test.


``` r
head(Indometh)
## Grouped Data: conc ~ time | Subject
##   Subject time conc
## 1       1 0.25 1.50
## 2       1 0.50 0.94
## 3       1 0.75 0.78
## 4       1 1.00 0.48
## 5       1 1.25 0.37
## 6       1 2.00 0.19
```

Next, type `plot(Indometh)` into the R Console. This will plot the indomethacin dataset in a basic way.


``` r
plot(Indometh)
```

<img src="index_files/figure-html/unnamed-chunk-8-1.png" width="100%" style="display: block; margin: auto;" />

If both the above commands (`head(Indometh)` and `plot(Indometh)`) worked and there were no error messages during installation, then you should be ready to proceed.

<!-- ## tidyverse

For us to run our analyses, we need to install a set of add-on functions that expand R's capabilities. These functions are collected in something called the tidyverse, a very well-known and widely-used R package. You do not need to manually download anything to complete this installation - R will do it for you. In the R Console, type `install.packages("tidyverse", repos = "https://cran.us.r-project.org")` to install the tidyverse. Let's try it:

RSudio might ask you: "Do you want to install from sources the packages which need compilation? (Yes/no/cancel)", for now, type `no` and press enter.


``` r
install.packages("tidyverse", repos = "https://cran.us.r-project.org")
```

Let's make sure your version of the tidyverse is installed correctly. To do this, we will load the `tidyverse` library/package inside of an R session. We can do this using `library(tidyverse)`. Let's try it:

``` r
library(tidyverse)
```

If the library load correctly - then you are set to go! If not, try updating your R / RStudio installations, the re installing the tidyverse. If this still fails, please contact me. -->

## TeX {-}

In this class we will generate high quality reports suitable for submission to supervisors, academic journals, etc. For this, we need the typesetting engine TeX. There are a few ways to do this. The easiest way is using the following commands:


``` r
install.packages(c('tinytex', 'rmarkdown'))
```

If you are on Mac, you may get an error about "not being able to write to a path" or something like that. In that case you probably need to open your terminal and run the following two commands: 

```
sudo chown -R \`whoami\`:admin /usr/local/bin
```

  and then 

```
~/Library/TinyTeX/bin/\*/tlmgr path add
```

Then, on both Mac and PC, you then need to do:

``` r
tinytex::install_tinytex()
```

Other options are: if you have Windows, download and install [MikTeX](https://miktex.org/download). If you have OSX, you can download and install [MacTeX](https://www.tug.org/mactex/morepackages.html).

## phylochemistry {-}

In addition to the tidyverse, there are a variety of other packages we will need, as well as some datasets and custom functions. These call all be loaded by doing the following.

First, attempt to load phylochemistry, if you are on Windows, be sure you've opened RStudio as an administrator (right click, "run as administrator"):

``` r
source("https://thebustalab.github.io/phylochemistry/phylochemistry.R")
```

The first time you try this, it will very likely say: "You need to install the following packages before proceeding […] Is it okay if phylochemistry installs them for you?" You should say "yes".

<!-- end -->

# (PART) DATA VISUALIZATION

<!-- start data visualization -->


# data visualization I {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/tufte_train.jpeg" width="100%" style="display: block; margin: auto;" />

Visualization is one of the most fun parts of working with data. In this section, we will jump into visualization as quickly as possible - after just a few prerequisites. Please note that data visualization is a whole field in and of itself (just google "data visualization" and see what happens). Data visualization is also rife with "trendy" visuals, misleading visuals, and visuals that look cool but don't actually communicate much information. We will touch on these topics briefly, but will spend most of our time practicing how to represent our data in intuitive and interpretable ways. Let's get started!

## {-}

## objects {-}

In R, data is stored in objects. You can think of objects as if they were "files" inside an R session. `phylochemistry` provides a variety of objects for us to work with. Let's look at how to create an object. For this, we can use an arrow: `<-` . The arrow will take something and store it inside an object. For example:


``` r
new_object <- 1
```

Now we've got a new object called `new_object`, and inside of it is the number 1. To look at what's inside an object, we can simply type the name of the object into the console:


``` r
new_object 
## [1] 1
```

Easy! Let's look at one of the objects that comes with our class code base. What are the dimensions of the "algae_data" data set?


``` r
algae_data
## # A tibble: 180 × 5
##    replicate algae_strain harvesting_regime chemical_species
##        <dbl> <chr>        <chr>             <chr>           
##  1         1 Tsv1         Heavy             FAs             
##  2         1 Tsv1         Heavy             saturated_Fas   
##  3         1 Tsv1         Heavy             omega_3_polyuns…
##  4         1 Tsv1         Heavy             monounsaturated…
##  5         1 Tsv1         Heavy             polyunsaturated…
##  6         1 Tsv1         Heavy             omega_6_polyuns…
##  7         1 Tsv1         Heavy             lysine          
##  8         1 Tsv1         Heavy             methionine      
##  9         1 Tsv1         Heavy             essential_Aas   
## 10         1 Tsv1         Heavy             non_essential_A…
## # ℹ 170 more rows
## # ℹ 1 more variable: abundance <dbl>
```

## functions {-}

Excellent - we've got data. Now we need to manipulate it. For this we need functions:

* A function is a command that tells R to perform an action!
* A function begins and ends with parentheses: `this_is_a_function()`
* The stuff inside the parentheses are the details of how you want the function to perform its action: `run_this_analysis(on_this_data)`

Let's illustrate this with an example. `algae_data` is a pretty big object. For our next chapter on visualization, it would be nice to have a smaller dataset object to work with. Let's use another `tidyverse` command called `filter` to filter the `algae_data` object. We will need to tell the filter command what to filter out using "logical predicates" (things like equal to: `==`, less than: `<`, greater than: `>`, greater-than-or-equal-to: `<=`, etc.). Let's filter `algae_data` so that only rows where the `chemical_species` is equal to `FAs` (fatty acids) is preserved. This will look like `chemical_species == "FAs"`. Here we go:


``` r
filter(algae_data, chemical_species == "FAs")
## # A tibble: 18 × 5
##    replicate algae_strain harvesting_regime chemical_species
##        <dbl> <chr>        <chr>             <chr>           
##  1         1 Tsv1         Heavy             FAs             
##  2         2 Tsv1         Heavy             FAs             
##  3         3 Tsv1         Heavy             FAs             
##  4         1 Tsv1         Light             FAs             
##  5         2 Tsv1         Light             FAs             
##  6         3 Tsv1         Light             FAs             
##  7         1 Tsv2         Heavy             FAs             
##  8         2 Tsv2         Heavy             FAs             
##  9         3 Tsv2         Heavy             FAs             
## 10         1 Tsv2         Light             FAs             
## 11         2 Tsv2         Light             FAs             
## 12         3 Tsv2         Light             FAs             
## 13         1 Tsv11        Heavy             FAs             
## 14         2 Tsv11        Heavy             FAs             
## 15         3 Tsv11        Heavy             FAs             
## 16         1 Tsv11        Light             FAs             
## 17         2 Tsv11        Light             FAs             
## 18         3 Tsv11        Light             FAs             
## # ℹ 1 more variable: abundance <dbl>
```

Cool! Now it's just showing us the 18 rows where the chemical_species is fatty acids (FAs). Let's write this new, smaller dataset into a new object. For that we use `<-`, remember?


``` r
algae_data_small <- filter(algae_data, chemical_species == "FAs")
algae_data_small
## # A tibble: 18 × 5
##    replicate algae_strain harvesting_regime chemical_species
##        <dbl> <chr>        <chr>             <chr>           
##  1         1 Tsv1         Heavy             FAs             
##  2         2 Tsv1         Heavy             FAs             
##  3         3 Tsv1         Heavy             FAs             
##  4         1 Tsv1         Light             FAs             
##  5         2 Tsv1         Light             FAs             
##  6         3 Tsv1         Light             FAs             
##  7         1 Tsv2         Heavy             FAs             
##  8         2 Tsv2         Heavy             FAs             
##  9         3 Tsv2         Heavy             FAs             
## 10         1 Tsv2         Light             FAs             
## 11         2 Tsv2         Light             FAs             
## 12         3 Tsv2         Light             FAs             
## 13         1 Tsv11        Heavy             FAs             
## 14         2 Tsv11        Heavy             FAs             
## 15         3 Tsv11        Heavy             FAs             
## 16         1 Tsv11        Light             FAs             
## 17         2 Tsv11        Light             FAs             
## 18         3 Tsv11        Light             FAs             
## # ℹ 1 more variable: abundance <dbl>
```

Here are a variety of ways to filter:

`filter(<data>, <variable> < 18)` ## less than 18

`filter(<data>, <variable> <= 18)` ## less than or equal to 18

`filter(<data>, <variable> > 18)` ## greater than 18

`filter(<data>, <variable> >= 18)` ## greater than or equal to 18

`filter(<data>, <variable> == 18)` ## equals than 18

`filter(<data>, <variable> != 18)` ## not equal to 18

`filter(<data>, <variable> == 18 | <variable> == 19)` ## equal to 18 or 19

`filter(<data>, <variable> %in% c(18, 19, 20)` ## equal to 18 or 19 or 20

## ggplot & geoms {-}

Now we have a nice, small table that we can use to practice data visualization. For visualization, we're going to use `ggplot2` - a powerful set of commands for plot generation. 

There are three steps to setting up a ggplot:

1. **Define the data you want to use.**

We do this using the ggplot function's data argument. When we run that line, it just shows a grey plot space. Why is this? It's because all we've done is told ggplot that (i) we want to make a plot and (ii) what data should be used. We haven't explained how to represent features of the data using ink.


``` r
ggplot(data = algae_data_small)
```

<img src="index_files/figure-html/unnamed-chunk-20-1.png" width="50%" style="display: block; margin: auto;" />

2. **Define how your variables map onto the axes.**

This is called aesthetic mapping and is done with the `aes()` function. `aes()` should be placed inside the `ggplot` command. Now when we run it, we get our axes!


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance))
```

<img src="index_files/figure-html/unnamed-chunk-21-1.png" width="50%" style="display: block; margin: auto;" />

3. **Use geometric shapes to represent other variables in your data.**

Map your variables onto the geometric features of the shapes. To define which shape should be used, use a `geom_*` command. Some options are, for example, `geom_point()`, `geom_boxplot()`, and `geom_violin()`. These functions should be added to your plot using the `+` sign. We can use a new line to keep the code from getting too wide, just make sure the `+` sign is at the end fo the top line. Let's try it:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_point()
```

<img src="index_files/figure-html/unnamed-chunk-22-1.png" width="50%" style="display: block; margin: auto;" />

In the same way that we mapped variables in our dataset to the plot axes, we can map variables in the dataset to the geometric features of the shapes we are using to represent our data. For this, again, use `aes()` to map your variables onto the geometric features of the shapes:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_point(aes(color = harvesting_regime))
```

<img src="index_files/figure-html/unnamed-chunk-23-1.png" width="50%" style="display: block; margin: auto;" />

In the plot above, the points are a bit small, how could we fix that? We can modify the features of the shapes by adding additional arguments to the `geom_*()` functions. To change the size of the points created by the `geom_point()` function, this means that we need to add the `size = ` argument. IMPORTANT! Please note that when we map a feature of a shape to a *variable* in our data(as we did with color/harvesting regime, above) then it goes *inside* aes(). In contrast, when we map a feature of a shape to a *constant*, it goes *outside* aes(). Here's an example:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_point(aes(color = harvesting_regime), size = 5)
```

<img src="index_files/figure-html/unnamed-chunk-24-1.png" width="50%" style="display: block; margin: auto;" />

One powerful aspect of `ggplot` is the ability to quickly change mappings to see if alternative plots are more effective at bringing out the trends in the data. For example, we could modify the plot above by switching how harvesting_regime is mapped:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_point(aes(size = harvesting_regime), color = "black")
```

<img src="index_files/figure-html/unnamed-chunk-25-1.png" width="50%" style="display: block; margin: auto;" />

** Important note: Inside the `aes()` function, map aesthetics (the features of the geom's shape) to a *variable*. Outside the `aes()` function, map aesthetics to *constants*. You can see this in the above two plots - in the first one, color is inside `aes()` and mapped to the variable called harvesting_regime, while size is outside the `aes()` call and is set to the constant 5. In the second plot, the situation is reversed, with size being inside the `aes()` function and mapped to the variable harvesting_regime, while color is outside the `aes()` call and is mapped to the constant "black".

We can also stack geoms on top of one another by using multiple `+` signs. We also don't have to assign the same mappings to each geom.


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_violin() +
  geom_point(aes(color = harvesting_regime), size = 5)
```

<img src="index_files/figure-html/unnamed-chunk-26-1.png" width="50%" style="display: block; margin: auto;" />

As you can probably guess right now, there are lots of mappings that can be done, and lots of different ways to look at the same data!


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_violin(aes(fill = algae_strain)) +
  geom_point(aes(color = harvesting_regime, size = replicate))
```

<img src="index_files/figure-html/unnamed-chunk-27-1.png" width="50%" style="display: block; margin: auto;" />


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_boxplot()
```

<img src="index_files/figure-html/unnamed-chunk-28-1.png" width="50%" style="display: block; margin: auto;" />

## markdown {-}

Now that we are able to filter our data and make plots, we are ready to make reports to show others the data processing and visualization that we are doing. For this, we will use R Markdown. You can open a new markdown document in RStudio by clicking: `File -> New File -> R Markdown`. You should get a template document that compiles when you press "knit".

Customize this document by modifying the title, and add `author: "your_name"` to the header. Delete all the content below the header, then compile again. You should get a page that is blank except for the title and the author name.

You can think of your markdown document as a stand-alone R Session. This means you will need to load our class code base into each new markdown doument you create. You can do this by adding a "chunk" or R code. That looks like this:

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/markdown_1.png" width="100%" style="display: block; margin: auto;" />

You can compilie this document into a pdf. We can also run R chunks right inside the document and create figures. You should notice a few things when you compile this document:

1. Headings: When you compile that code, the "# My first analysis" creates a header. You can create headers of various levels by increasing the number of hashtags you use in front of the header. For example, "## Part 1" will create a subheading, "### Part 1.1" will create a sub-subheading, and so on.

2. Plain text: Plain text in an R Markdown document creates a plan text entry in your compiled document. You can use this to explain your analyses and your figures, etc.

3. You can modify the output of a code chunk by adding arguments to its header. Useful arguments are fig.height, fig.width, and fig.cap. Dr. Busta will show you how to do this in class.

<!-- ## exercises {-}

In this set of exercises we're going to practice filtering and plotting data in R Markdown. We're going to work with two datasets: (i) algae_data and (ii) alaska_lake_data. **For these exercises, you will write your code and answers to all questions in an R Markdown report, compile it as a pdf, and submit it on Canvas. If you have any questions please let me know**

Some pointers:

- If your code goes off the page, don't be afraid to wrap it across multiple lines, as shown in some of the examples.

- Don't be afraid to put the variable with the long elements / long text on the y-axis and the continuous variable on the x-axis.

### algae {-}

1. You will have `algae_data` stored in an object called `algae_data` as soon as you run `source("https://thebustalab.github.io/phylochemistry/phylochemistry.R")`. For this question, filter the data so that only entries are shown for which the `chemical_species` is "FAs" (remember that quotes are needed around FAs here!). What are the dimensions (i.e. number of rows and columns) of the resulting dataset?



2. Now filter the original dataset (`algae_data`) so that only entries for the `algae_strain` "Tsv1" are shown. What are the dimensions of the resulting dataset?



3. Now filter the original dataset (`algae_data`) so that only entries with an abundance greater than 250 are shown. Note that `>` can be used in the filter command instead of `==`, and that numbers inside a filter command do not require quotes around them. What are the dimensions of the resulting dataset?



4. Use the original dataset (`algae_data`) to make a ggplot that has `algae_strain` on the x axis and `abundance` on the y axis. Remember about `aes()`. Use points (`geom_point()`) to represent each compound. You don't need to color the points. Which algae strain has the most abundant compound out of all the compounds in the dataset?



5. Make a ggplot that has `abundance` on the x axis and `chemical_species` on the y axis. Use points to represent each compound. You don't need to color the points. Generally speaking, what are the two most abundant classes of chemical species in these algae strains? (FAs/Fas stand for fatty acids, AAs/Aas stand for amino acids.)



6. I am going to show you an example of how you can filter and plot at the same time. To do this, we nest the filter command inside ggplot's data argument:


``` r
ggplot(
  data = filter(algae_data, chemical_species == "essential_Aas"),
  aes(x = algae_strain, y = abundance)) +
geom_point()
```

<img src="index_files/figure-html/unnamed-chunk-35-1.png" width="100%" style="display: block; margin: auto;" />

Using the above as a template, make a plot that shows just `omega_3_polyunsaturated_Fas`, with algae_strain on the x axis, and abundance on the y axis. Color the points so that they correspond to `harvesting_regime`. Remember that mapping a feature of a shape onto a variable must be done inside `aes()`. Change the plot so that all the points are size = 5. Remember that mapping features of a shape to a constant needs to be done outside `aes()`. Which harvesting regime leads to higher levels of `omega_3_polyunsaturated_Fas`?



7. Use a combination of filtering and plotting to show the abundance of the different chemical species in just the `algae_strain` called "Tsv1". Use an x and y axis, as well as points to represent the measurements. Make point size correspond to the replicate, and color the points according to harvesting regime.



8. Make a plot that checks to see which `chemical_species` were more abundant under light as opposed to heavy `harvesting_regime` in all three replicates. Use filtered data so that just one `algae_strain` is shown, an x and a y axis, and points to represent the measurements. Make the points `size = 5` and also set the point's `alpha = 0.6`. The points should be colored according to harvesting_regime. Make 3 plots, one for each strain of algae.







9. Take the code that you made for the question above. Remove the filtering. Add the following line to the end of the plot: `facet_grid(.~algae_strain)`. Remember that adding things to plots is done with the `+` sign, so your code should look something like:


``` r
ggplot(data = algae_data, aes(x = <something>, y = <something else>)) +
  geom_point(aes(<some things>), <some here too>) +
  facet_grid(.~algae_strain)
```



Also try, instead of `facet_grid(.~algae_strain)`, `facet_grid(algae_strain~.)` at the end of you plot command. (note the swap in the position of the `.~` relative to `algae_strain`). This means your code should look something like:


``` r
ggplot(data = algae_data, aes(x = <something>, y = <something else>)) +
  geom_point(aes(<some things>), <some here too>) +
  facet_grid(algae_strain~.)
```



What advantages does this one extra line (i.e. facet_grid) provide over what you had to do in question 8?

### alaska lakes {-}

1. Use R to view the first few lines of the `alaska_lake_data` dataset. Do your best to describe, in written format, the kind of data that are in this data set.



2. How many variables are in the Alaska lakes dataset?



3. Filter the data set so only meausurements of free elements (i.e. element_type is "free") are shown. Remember, it's `==`, not `=`. What are the dimensions of the resulting dataset?



4. Make a plot that shows the water temperatures of each lake. Don't worry if you get a warning message from R about "missing values". Which is the hottest lake? The coolest?



5. Make a plot that shows the water temperature of each lake. The x axis should be `park`, the y axis `water temp`. Add geom_violin() to the plot first, then geom_point(). Make the points size = 5. Color the points according to water_temp. Which park has four lakes with very similar temperatures?



6. From the plot you made for question 5, it should be apparent that there is one lake in NOAT that is much warmer than the others. Filter the data so that only entries from `park == "NOAT"` are shown (note the double equals sign and the quotes around NOAT...). Combine this filtering with plotting and use geom_point() to make a plot that shows which specific lake that is.



7. Make a plot that shows which lake has the highest abundance of sulfur.



8. Make a plot that uses geom_point(). Set the "shape" aesthetic of the points to 21, i.e. `geom_point(aes(...), shape = 21)`. This gives you access to a new aesthetics: `fill`. It also changes the behaviour of the `color` aesthetic slightly, in that it now controls border color, not the internal color. Here is an example (though it doesn't make a very nice plot):


``` r
ggplot(
  data = filter(alaska_lake_data, lake == "Lake_Narvakrak"),
  aes(x = lake, y = mg_per_L)
) +
  geom_point(
    shape = 21, size = 10,
    color = "black", fill = "green"
  )
```

<img src="index_files/figure-html/unnamed-chunk-52-1.png" width="100%" style="display: block; margin: auto;" />

Now we have lots of aesthetics we can map to: x, y, size, color, and fill (leave shape set to 21 for now). Make a plot of your own design. It should include filtering, and all the aesthetics listed above, though whether you map them to a variable or a constant is up to you.



When you are done with this plot, take a screen shot of it. Go to [THIS GOOGLE SHEET](https://docs.google.com/presentation/d/1G0BJ_qye9a_HAPLktFytj66qSj20BjoUOTKtjmCyuN0/edit?usp=sharing), make a slide for yourself (you don't have to include your name), and paste your screen shot there. Add a small caption that explains how your variables are mapped. -->

<!-- end -->

<!-- start data visualization II -->


# data visualization II {-}

## {-}

## more geoms {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/chart_suggestions.jpeg" width="100%" style="display: block; margin: auto;" />

We've looked at how to filter data and map variables in our data to geometric shapes to make plots. Let's have a look at a few more things. For these examples, we're going to use the data set called `solvents`. In these examples, I'd like to introduce you to two new geoms. The first `geom_smooth()` is used when there are two continuous variables. It is particularly nice when geom_point() is stacked on top of it.


``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point()
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-55-1.png" width="100%" style="display: block; margin: auto;" />

Also, please be aware of `geom_tile()`, which is nice for situations with two discrete variables and one continuous variable. `geom_tile()` makes what are often referred to as heat maps. Note that `geom_tile()` is somewhat similar to `geom_point(shape = 21)`, in that it has both `fill` and `color` aesthetics that control the fill color and the border color, respectively.


``` r
ggplot(
  data = filter(algae_data, harvesting_regime == "Heavy"),
  aes(x = algae_strain, y = chemical_species)
) + 
  geom_tile(aes(fill = abundance), color = "black", size = 1)
```

<img src="index_files/figure-html/unnamed-chunk-56-1.png" width="100%" style="display: block; margin: auto;" />

These examples should illustrate that there is, to some degree, correspondence between the type of data you are interested in plotting (number of discrete and continuous variables) and the types of geoms that can effectively be used to represent the data.

## facets {-}

As alluded to in Exercises 1, it is possible to map variables in your dataset to more than the geometric features of shapes (i.e. geoms). One very common way of doing this is with facets. Faceting creates small multiples of your plot, each of which shows a different subset of your data based on a categorical variable of your choice. Let's check it out.

Here, we can facet in the horizontal direction:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(.~replicate)
```

<img src="index_files/figure-html/unnamed-chunk-57-1.png" width="100%" style="display: block; margin: auto;" />

We can facet in the vertical direction:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(replicate~.)
```

<img src="index_files/figure-html/unnamed-chunk-58-1.png" width="100%" style="display: block; margin: auto;" />

And we can do both at the same time:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(harvesting_regime~replicate)
```

<img src="index_files/figure-html/unnamed-chunk-59-1.png" width="100%" style="display: block; margin: auto;" />

Faceting is a great way to describe more variation in your plot without having to make your geoms more complicated. For situations where you need to generate lots and lots of facets, consider `facet_wrap` instead of `facet_grid`:



``` r
ggplot(data = algae_data, aes(x = replicate, y = algae_strain)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_wrap(chemical_species~.)
```

<img src="index_files/figure-html/unnamed-chunk-60-1.png" width="100%" style="display: block; margin: auto;" />

## scales {-}

Every time you define an aesthetic mapping (e.g. aes(x = algae_strain)), you are defining a new scale that is added to your plot. You can control these scales using the `scale_*` family of commands. Consider our faceting example above. In it, we use `geom_tile(aes(fill = abundance))` to map the abundance variable to the fill aesthetic of the tiles. This creates a scale called `fill` that we can adjust using `scale_fill_*`. In this case, fill is mapped to a continuous variable and so the fill scale is a color gradient. Therefore, `scale_fill_gradient()` is the command we need to change it. Remember that you could always type `?scale_fill_` into the console and it will help you find relevant help topics that will provide more detail. Another option is to google: "How to modify color scale ggplot geom_tile", which will undoubtedly turn up a wealth of help.


``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(harvesting_regime~replicate) +
  scale_fill_gradient(low = "white", high = "black") +
  theme_classic()
```

<img src="index_files/figure-html/unnamed-chunk-61-1.png" width="100%" style="display: block; margin: auto;" />

One particularly useful type of scale are the color scales provided by RColorBrewer:


``` r
display.brewer.all()
```

<img src="index_files/figure-html/unnamed-chunk-62-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ggplot(mtcars) +
  geom_point(
    aes(x = mpg, y = factor(cyl), fill = factor(carb)), 
    shape = 21, size = 6
  ) +
  scale_fill_brewer(palette = "Set1")
```

<img src="index_files/figure-html/unnamed-chunk-63-1.png" width="100%" style="display: block; margin: auto;" />
  
## themes {-}
  
So far we've just looked at how to control the means by which your *data* is represented on the plot. There are also components of the plot that are, strictly speaking, not *data* per se, but rather non-data ink. These are controlled using the `theme()` family of commands. There are two ways to go about this.

`ggplot` comes with a handful of built in "complete themes". These will change the appearance of your plots with respect to the non-data ink. Compare the following plots:


``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point() +
  theme_classic()
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-64-1.png" width="100%" style="display: block; margin: auto;" />


``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point() +
  theme_dark()
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-65-1.png" width="100%" style="display: block; margin: auto;" />
  

``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point() +
  theme_void()
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-66-1.png" width="100%" style="display: block; margin: auto;" />

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/what_is_ggplot.jpeg" width="100%" style="display: block; margin: auto;" />

You can also change individual components of themes. This can be a bit tricky, but it's all explained if you run `?theme()`. Hare is an example (and google will provide many, many more).


``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point() +
  theme(
    text = element_text(size = 20, color = "black")
  )
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-68-1.png" width="100%" style="display: block; margin: auto;" />

Last, here is an example of combining `scale_*` and `theme_*` with previous commands to really get a plot looking sharp.


``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth(color = "#4daf4a") +
  scale_x_continuous(
    name = "Boiling Point", breaks = seq(0,200,25), limits = c(30,210)
  ) +
  scale_y_continuous(
    name = "Vapor Pressure", breaks = seq(0,600,50)
  ) +
  geom_point(color = "#377eb8", size = 4, alpha = 0.6) +
  theme_bw() +
  theme(
    axis.text = element_text(color = "black"),
    text = element_text(size = 16, color = "black")
  )
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-69-1.png" alt="Vapor pressure as a function of boiling point. A scatter plot with trendline showing the vapor pressure of thirty-two solvents (y-axis) a as a function of their boiling points (x-axis). Each point represents the boiling point and vapor pressure of one solvent. Data are from the 'solvents' dataset used in UMD CHEM5725." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-69)Vapor pressure as a function of boiling point. A scatter plot with trendline showing the vapor pressure of thirty-two solvents (y-axis) a as a function of their boiling points (x-axis). Each point represents the boiling point and vapor pressure of one solvent. Data are from the 'solvents' dataset used in UMD CHEM5725.</p>
</div>

## subplots {-}

We can make subplots using the `cowplot` package, which comes with the `source()` command. Let's see:


``` r
library(patchwork)
plot1 <-  ggplot(
            filter(alaska_lake_data, element_type == "free")
          ) +
          geom_violin(aes(x = park, y = mg_per_L)) + theme_classic() +
          ggtitle("A")

plot2 <-  ggplot(
            filter(alaska_lake_data, element_type == "bound")
          ) +
          geom_violin(aes(x = park, y = mg_per_L)) + theme_classic() +
          ggtitle("B")

plot3 <-  ggplot(
            filter(alaska_lake_data, element == "C")
          ) +
          geom_violin(aes(x = park, y = mg_per_L)) + theme_classic() +
          coord_flip() + ggtitle("C")

plot_grid(plot_grid(plot1, plot2), plot3, ncol = 1)
```

<img src="index_files/figure-html/unnamed-chunk-70-1.png" width="100%" style="display: block; margin: auto;" />

<!-- ## exercises {-}

In this set of exercises we're going to practice making more plots using the dataset `solvents`. Well, you don't have to use `solvents`, you could use something else if you want, but `solvents` is a fun one to explore. Since you are now familiar with filtering and plotting data, the prompts in this assignment are going to be relatively open ended - I do not care what variables you map to x, y, fill, color, etc. Rather, I expect your submission to demonstrate to me that you have explored each of the new topics covered in the previous chapter. This includes geoms beyond `geom_point()` and `geom_violin()`, facets, scale modifications, and theme adjustments. Be creative! Explore the solvents dataset. Find something interesting! **Show me that you have mastered this material.** Don't forget about the ggplot cheat sheet (see the "Links" section in this book).

As before, for these exercises, you will write your code and answers to any questions in the Script Editor window of your RStudio as an R Markdown document. You will compile that file as a pdf and submit it on Canvas. If you have any questions please let me know.

Some pointers:

- If your code goes off the page, don't be afraid to wrap it across multiple lines, as shown in some of the examples in the previous set of exercises.

- Don't be afraid to put the variable with the long elements / long text on the y-axis and the continuous variable on the x-axis.

1. Create a plot that has x and y axes that are continuous variables. Add to this plot `facet_grid`, and specify that the facets should be based on a categorical variable (ideally a categorical variable with a small number of total categories). Now make two versions of that plot, one that uses the `scales = "free"` feature of `facet_grid` and a second the other does not (i.e. one should use `facet_grid(<things>)`, while the other uses `facet_grid(<things>, scales = "free")`). Write a single caption that describes *both* plots, highlighting the advantages provided by each plot over the other. For additional tips on writing captions, please see the "Writing" chapter in this book.

2. Using a continuous variable on one axis and a discrete (categorical) variable on the other, create two plots that are identical except that one uses `geom_point()`, while the other uses `geom_jitter()`. Write a single caption that describes *both* plots. The caption should highlight the differences between these two plots and it should describe case(s) in which you think it would be appropriate to use `geom_jitter()` over `geom_point()`.

3. Make a plot that has four aesthetic mappings (x and y mappings count). Use the `scales_*` family of commands to modify some aspect of each scale create by the four mappings. Hint: some scales are somewhat tricky to modify (alpha, linetype, ...), and some scales are easier to modify (x, y, color, fill, shape). You may need to use some google searches to help you. Queries along the lines of "how to modify point color in ggplot" should direct you to a useful resource.

4. Make a duplicate of the plot you created in the previous question and modify its theme.

5. Identify a relationship between two variables in the dataset. Create a plot that is optimized (see note) to highlight the features of this relationship. Write a short caption that describes the plot *and* the trend you've identified and highlighted. Note: I realize that the word "optimize" is not clearly defined here. That's ok! You are the judge of what is optimized and what is not. Use your caption to make a case for *why* your plot is optimized.

6. Watch [this video](https://www.youtube.com/watch?v=LFDbqw2xPbQ) on bar plots. Add a section to the end of the R Markdown document you made for Part 2 that describes the problem outlined in the video and one potential solution to the problem. -->

## {-}

## further reading {-}

There is a [handy cheat sheet](https://www.maths.usyd.edu.au/u/UG/SM/STAT3022/r/current/Misc/data-visualization-2.1.pdf) that can help you identify the right geom for your situation. Please keep this cheat sheet in mind for your future plotting needs...

For additional explanations of ggplot2: [ggplot2-book](https://ggplot2-book.org/).

Check out some of the incredible geoms that are easy to access using R and ggplot2: [R Graph Gallery](https://r-graph-gallery.com/). Use these to make your figures attractive and easy to interpret!

For a challenge, try implementing these awesome color scales: [Famous R Color Palettes](https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/). Note that some of these are optimized for colorblind individuals and that other are optimized for continuous hue gradients, etc.

For a list of data visualization sins: [Friends Don't Let Friends](https://github.com/cxli233/FriendsDontLetFriends). Some interesting things in here!

For more information on data visualization and graphics theory, check out the works by Edward Tufte: [Edward Tufte](https://www.edwardtufte.com/tufte/). A digital text that covers similar topics is here: [Look At Data] (https://socviz.co/lookatdata.html).

Some examples of award winning data visualization: [Information Is Beautiful Awards](https://www.informationisbeautifulawards.com/showcase?award=2019&type=awards) and [Data Vis Inspiration](https://www.dataviz-inspiration.com/).

Additional color palettes: [MetBrewer](https://github.com/BlakeRMills/MetBrewer) and [Paletteer](https://github.com/EmilHvitfeldt/paletteer).

<!-- end -->

<!-- start data visualization III -->

# data visualization III {-}

## {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/datavis3.png" width="100%" style="display: block; margin: auto;" />

## advanced plots {-}

### 3D scatter plots {-}

`phylochemistry` contains a function to help you make somewhat decent 3D scatter plots. Let's look at an example (see below). For this, we use the function `points3D`. Se give it a `data` argument that gives it vectors of data that should be on the x, y, and z axes, along with a vector that uniquely identifies each observation. We also tell it the angle of the z axis that we want, the integer to which ticks should be rounded, and the tick intervals. The function returns data that we can pass to ggplot to make a 3D plot.


``` r
pivot_wider(hawaii_aquifers, names_from = "analyte", values_from = "abundance") %>%
  mutate(sample_unique_ID = paste0(aquifer_code, "_", well_name)) -> aquifers

output <- points3D(
  data = data.frame(
    x = aquifers$SiO2,
    y = aquifers$Cl,
    z = aquifers$Mg,
    sample_unique_ID = aquifers$sample_unique_ID
  ),
  angle = pi/2.4,
  tick_round = 10,
  x_tick_interval = 10,
  y_tick_interval = 20,
  z_tick_interval = 20
)

str(output)
## List of 6
##  $ grid          :'data.frame':	14 obs. of  4 variables:
##   ..$ y   : num [1:14] 0 0 0 0 0 ...
##   ..$ yend: num [1:14] 96.6 96.6 96.6 96.6 96.6 ...
##   ..$ x   : num [1:14] 10 20 30 40 50 ...
##   ..$ xend: num [1:14] 35.9 45.9 55.9 65.9 75.9 ...
##  $ ticks         :'data.frame':	37 obs. of  4 variables:
##   ..$ y   : num [1:37] 0 0 0 0 0 0 0 0 0 0 ...
##   ..$ yend: num [1:37] 1.93 1.93 1.93 1.93 1.93 ...
##   ..$ x   : num [1:37] 10 20 30 40 50 60 70 80 10 20 ...
##   ..$ xend: num [1:37] 10.5 20.5 30.5 40.5 50.5 ...
##  $ labels        :'data.frame':	29 obs. of  3 variables:
##   ..$ y    : num [1:29] -11.2 -11.2 -11.2 -11.2 -11.2 -11.2 -11.2 -11.2 0 20 ...
##   ..$ x    : num [1:29] 7.2 17.2 27.2 37.2 47.2 57.2 67.2 77.2 4.4 4.4 ...
##   ..$ label: num [1:29] 10 20 30 40 50 60 70 80 0 20 ...
##  $ axes          :'data.frame':	3 obs. of  4 variables:
##   ..$ x   : num [1:3] 10 10 80
##   ..$ xend: num [1:3] 80 10 106
##   ..$ y   : num [1:3] 0 0 0
##   ..$ yend: num [1:3] 0 280 96.6
##  $ point_segments:'data.frame':	106 obs. of  4 variables:
##   ..$ x   : num [1:106] 13 33.1 39.4 53.1 22.5 ...
##   ..$ xend: num [1:106] 13 33.1 39.4 53.1 22.5 ...
##   ..$ y   : num [1:106] 27.3 81.6 82.6 109.5 19.5 ...
##   ..$ yend: num [1:106] 7.34 11.59 12.56 15.45 5.51 ...
##  $ points        :'data.frame':	106 obs. of  3 variables:
##   ..$ x               : num [1:106] 13 33.1 39.4 53.1 22.5 ...
##   ..$ y               : num [1:106] 27.3 81.6 82.6 109.5 19.5 ...
##   ..$ sample_unique_ID: chr [1:106] "aquifer_1_Alewa_Heights_Spring" "aquifer_1_Beretania_High_Service" "aquifer_1_Beretania_Low_Service" "aquifer_1_Kuliouou_Well" ...
```

The output from points3D contains a grid, axes, and ticks, which should all be plotted using geom_segment. It also contains points that should be plotted with geom_point, and point segments that should be plotted with geom_segement. We can take the output from points3D and join it with the original data, which will occurr according to our sample_unique_ID column. Then, we can also plot point metadata:


``` r
output$points <- left_join(output$points, aquifers)
## Joining with `by = join_by(sample_unique_ID)`
```

``` r
  
ggplot() +
  geom_segment(
    data = output$grid, aes(x = x, xend = xend, y = y, yend = yend),
    color = "grey80"
  ) +
  geom_segment(data = output$axes, aes(x = x, xend = xend, y = y, yend = yend)) +
  geom_segment(data = output$ticks, aes(x = x, xend = xend, y = y, yend = yend)) +
  geom_text(
    data = output$labels, aes(x = x, y = y, label = label),
    hjust = 0.5
  ) +
  geom_segment(
    data = output$point_segments,
    aes(x = x, xend = xend, y = y, yend = yend),
    linetype = "dotted", color = "black"
  ) +
  geom_point(
    data = output$points, aes(x = x, y = y, fill = aquifer_code),
    size = 3, shape = 21
  ) +
  theme_void() +
  scale_fill_manual(values = discrete_palette)
```

<img src="index_files/figure-html/unnamed-chunk-73-1.png" width="100%" style="display: block; margin: auto;" />

### marginal summaries {-}


``` r
i2 <- iris %>%
  mutate(Species2 = rep(c("A","B"), 75))
p <- ggplot(i2, aes(Sepal.Width, Sepal.Length, color = Species)) +
  geom_point()

p + geom_xsidedensity(aes(y=stat(density), xfill = Species), position = "stack")+
  geom_ysidedensity(aes(x=stat(density), yfill = Species2), position = "stack") +
  theme_bw() + 
  facet_grid(Species~Species2, space = "free", scales = "free") +
  labs(title = "FacetGrid", subtitle = "Collapsing All Side Panels") +
  ggside(collapse = "all") +
  scale_xfill_manual(values = c("darkred","darkgreen","darkblue")) +
  scale_yfill_manual(values = c("black","gold"))
```

<img src="index_files/figure-html/unnamed-chunk-74-1.png" width="100%" style="display: block; margin: auto;" />

### representing distributions {-}

You can also combine geoms to create more detailed representations of distributions:


``` r
mpg %>% filter(cyl %in% c(4,6,8)) %>%
  ggplot(aes(x = factor(cyl), y = hwy, fill = factor(cyl))) +
  ggdist::stat_halfeye(
    adjust = 0.5, justification = -0.2, .width = 0, point_colour = NA
  ) +
  geom_boxplot(width = 0.12, outlier.color = NA, alpha = 0.5) +
  ggdist::stat_dots(side = "left", justification = 1.1, binwidth = .25)
```

<img src="index_files/figure-html/unnamed-chunk-75-1.png" width="100%" style="display: block; margin: auto;" />

### venn digrams {-}


``` r
df <- data.frame(
  plant1 = sample(c(TRUE, FALSE), 24, replace = TRUE),
  plant2 = sample(c(TRUE, FALSE), 24, replace = TRUE),
  plant3 = sample(c(TRUE, FALSE), 24, replace = TRUE),
  attribute_name = sample(letters, 24, replace = FALSE)
)

vennAnalysis(df[,1:3]) %>%
  ggplot() +
    geom_circle(
      aes(x0 = x, y0 = y, r = r, fill = category),
      alpha = 0.4
    ) +
  scale_fill_brewer(palette = "Set1") +
  theme_void()
```

<img src="index_files/figure-html/unnamed-chunk-76-1.png" width="100%" style="display: block; margin: auto;" />


### ternary plots {-}


``` r
library(ggplot2)
library(ggtern)
alaska_lake_data %>%
  pivot_wider(names_from = "element", values_from = "mg_per_L") %>%
  ggtern(aes(
    x = Ca,
    y = S,
    z = Na,
    color = park,
    size = pH
    )) +
  geom_point() 
```

<img src="index_files/figure-html/unnamed-chunk-77-1.png" width="100%" style="display: block; margin: auto;" />


## map data {-}

### plotting boundaries {-}

There is a simple way to plot maps with ggplot. The map data comes with `ggplot2`! Let's have a look. See below some of the data sets included. Options included with ggplot are: `world`, `world2`, `usa`, `state` (US), `county` (US), `nz`, `italy`, and `france`. `geom_polygon()` is useful for plotting these, at (at least to me) seems more intuitive than `geom_map()`.


``` r
head(map_data("world"))
##        long      lat group order region subregion
## 1 -69.89912 12.45200     1     1  Aruba      <NA>
## 2 -69.89571 12.42300     1     2  Aruba      <NA>
## 3 -69.94219 12.43853     1     3  Aruba      <NA>
## 4 -70.00415 12.50049     1     4  Aruba      <NA>
## 5 -70.06612 12.54697     1     5  Aruba      <NA>
## 6 -70.05088 12.59707     1     6  Aruba      <NA>
```

``` r
head(map_data("state"))
##        long      lat group order  region subregion
## 1 -87.46201 30.38968     1     1 alabama      <NA>
## 2 -87.48493 30.37249     1     2 alabama      <NA>
## 3 -87.52503 30.37249     1     3 alabama      <NA>
## 4 -87.53076 30.33239     1     4 alabama      <NA>
## 5 -87.57087 30.32665     1     5 alabama      <NA>
## 6 -87.58806 30.32665     1     6 alabama      <NA>
```


``` r
head(map_data("county"))
##        long      lat group order  region subregion
## 1 -86.50517 32.34920     1     1 alabama   autauga
## 2 -86.53382 32.35493     1     2 alabama   autauga
## 3 -86.54527 32.36639     1     3 alabama   autauga
## 4 -86.55673 32.37785     1     4 alabama   autauga
## 5 -86.57966 32.38357     1     5 alabama   autauga
## 6 -86.59111 32.37785     1     6 alabama   autauga
```


``` r
head(map_data("france"))
##       long      lat group order region subregion
## 1 2.557093 51.09752     1     1   Nord      <NA>
## 2 2.579995 51.00298     1     2   Nord      <NA>
## 3 2.609101 50.98545     1     3   Nord      <NA>
## 4 2.630782 50.95073     1     4   Nord      <NA>
## 5 2.625894 50.94116     1     5   Nord      <NA>
## 6 2.597699 50.91967     1     6   Nord      <NA>
```

Cool! We can see that lat, lon, group, order, region, and subregion are included. That makes plotting easy. Note that `coord_map()` can help preserve aspect ratios:


``` r
ggplot(map_data("world")) +
  geom_point(aes(x = long, y = lat, color = group), size = 0.5) +
  theme_void() +
  coord_map()
```

<img src="index_files/figure-html/unnamed-chunk-82-1.png" width="100%" style="display: block; margin: auto;" />

Note that we can use `coord_map()` to do some pretty cool things!


``` r
ggplot(map_data("world")) +
  geom_point(aes(x = long, y = lat, color = group), size = 0.5) +
  theme_void() +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)
```

<img src="index_files/figure-html/unnamed-chunk-83-1.png" width="100%" style="display: block; margin: auto;" />

We can use filtering to produce maps of specific regions.


``` r
ggplot() +
  geom_polygon(
    data = filter(map_data("county"), region == "minnesota"),
    aes(x = long, y = lat, group = subregion, fill = subregion),
    color = "black"
  ) +
  theme_void() +
  coord_map()
```

<img src="index_files/figure-html/unnamed-chunk-84-1.png" width="100%" style="display: block; margin: auto;" />

### maps with plots {-}

Please note that the Great Lakes are in map_data()!


``` r
filter(map_data("lakes"), region == "Great Lakes", subregion == "Superior") %>%
    ggplot() +
      geom_path(aes(x = long, y = lat)) +
      coord_map() +
      theme_minimal()
```

<img src="index_files/figure-html/unnamed-chunk-85-1.png" width="100%" style="display: block; margin: auto;" />

We can clean up the map by making different groups for geom_path() whenever two consecutive points are far apart:


``` r
# Step 1: Filter and prepare your data for Lake Superior (though yes, that includes Michigan and Huron)
lake_superior <- map_data("lakes") %>%
  filter(region == "Great Lakes", subregion == "Superior") %>%
  arrange(order)

# Step 2: Calculate distances between consecutive points
lake_superior <- lake_superior %>%
  mutate(lag_long = lag(long),
         lag_lat = lag(lat),
         dist_to_prev = geosphere::distHaversine(cbind(long, lat), cbind(lag_long, lag_lat))) # distHaversine calculates distances

# Step 3: Define a threshold (e.g., 50 km) and create the "distance_group"
threshold <- 50000  # 50 km
lake_superior <- lake_superior %>%
  mutate(distance_group = cumsum(ifelse(dist_to_prev > threshold | is.na(dist_to_prev), 1, 0)))

# Step 4: Plot the map with `distance_group`
ggplot(lake_superior, aes(x = long, y = lat, group = distance_group)) +
  geom_path() +
  coord_map() +
  theme_minimal()
```

<img src="index_files/figure-html/unnamed-chunk-86-1.png" width="100%" style="display: block; margin: auto;" />

Now we could add some data. We could do something simple like plot total abundances as the size of a point:


``` r
lake_superior_PFAS <- readMonolist("/Users/bust0037/Documents/Science/Websites/pfas_data_private.csv")
lake_superior_PFAS %>%
  group_by(site, lon, lat) %>%
  summarize(total = sum(abundance)) -> lake_superior_PFAS_summarized

ggplot() +
  geom_path(
    data = filter(lake_superior, lat > 46, long < -84),
    aes(x = long, y = lat, group = distance_group)
  ) +
  geom_point(
    data =  lake_superior_PFAS_summarized,
    aes(x = lon, y = lat, size = total),
    color = "black"
  ) +
  coord_map() +
  theme_cowplot()
```

<img src="index_files/figure-html/unnamed-chunk-87-1.png" width="100%" style="display: block; margin: auto;" />

Or we could do something more sophisticated like add pie charts at each point:



``` r
lake_superior_PFAS <- readMonolist("/Users/bust0037/Documents/Science/Websites/pfas_data_private.csv")

grouped_by_site <- filter(lake_superior_PFAS, component == "PFBA")
site_less_than_90lon <- filter(grouped_by_site, lon <= -90)
site_more_than_90lon <- filter(grouped_by_site, lon >= -90)

unique_sites <- unique(lake_superior_PFAS$site)
dataframe_of_pies <- list()
for (i in 1:length(unique_sites)) { #i=1
  this_site <- filter(lake_superior_PFAS, site == unique_sites[i])
  this_site %>%
    ggplot()+
    geom_col(aes(x = 1, y = abundance, fill = class_name), color = "black") +
    coord_polar(theta = "y") +
    theme_void() +
    scale_fill_brewer(palette = "Set1", guide = "none") -> this_sites_pie
  dataframe_of_pies[[i]] <- tibble(x = this_site$lon[1], y = this_site$lat[1], plot = list(this_sites_pie))
}
dataframe_of_pies <- do.call(rbind, dataframe_of_pies)

ggplot() +
  geom_path(
    data = filter(lake_superior, lat > 46, long < -84),
    aes(x = long, y = lat, group = distance_group)
  ) +
  geom_point(
    data =  lake_superior_PFAS,
    aes(x = lon, y = lat),
    color = "black"
  ) +
  geom_plot(
    data = dataframe_of_pies, aes(x = x, y = y, label = plot),
    vp.width = 1/20, hjust = 0.5, vjust = 0.5, alpha = 0.5
  ) +
  geom_label_repel(data = site_less_than_90lon, aes(x = lon, y = lat, label = site), size = 2.5,min.segment.length = 0.01) +
  geom_label(data = site_more_than_90lon, aes(x = lon, y = lat, label = site), size = 2.5) +
  coord_map() +
  theme_cowplot()
```

<img src="index_files/figure-html/unnamed-chunk-88-1.png" width="100%" style="display: block; margin: auto;" />

You can also access a high resolution shoreline dataset for Lake Superior directly from the source() command as `lake_superior_shoreline`:


``` r
shore <- readMonolist("/Users/bust0037/Documents/Science/Websites/thebustalab.github.io/phylochemistry/sample_data/lake_superior_shoreline.csv")

wide_view <- ggplot(shore) +
    geom_point(aes(y = lat, x = lon), size = 0.01) +
    coord_map() +
    theme_minimal()

zoom_view <- ggplot(filter(shore, lat < 47.2, lat > 46.6, lon < -90)) +
    geom_point(aes(y = lat, x = lon), size = 0.01) +
    coord_map() +
    theme_minimal()
    
plot_grid(wide_view, zoom_view, nrow = 1, rel_widths = c(1,2))
```

<img src="index_files/figure-html/unnamed-chunk-89-1.png" width="100%" style="display: block; margin: auto;" />

## {-}

## further reading {-}

For more on plotting maps in R: [datavizplyr](https://datavizpyr.com/how-to-make-us-state-and-county-level-maps-in-r/)

For more advanced map plotting: [R Spatial](https://r-spatial.org/r/2018/10/25/ggplot2-sf.html)

For more on ternary plots: [ggtern](https://www.jstatsoft.org/article/view/v087c03)

<!-- ## exercises {-} -->

<!-- Using the `hawaii_aquifers` data set, please complete the following:

1. Choose one analyte and filter the data so only the rows for that analyte are shown.

2. Choose two of the aquifers. Are the mean abundances for your chosen analyte different in these two aquifers? Don't forget to test your data for normality and homogeneity of variance before selecting a statistical test. Use a plot to illustrate whether the means are similar or different.

3. Choose a second analyte, different from the first one you chose. Considering all the aquifers in the dataset, do any of them have the same abundance of this analyte? Again, don't forget about normality and homogeneity of variance tests. Use a plot to illustrate your answer.

4. Repeat #3 above, but switch the type of test used (i.e. use non-parametric if you used parametric for #3 and vice-versa). Compare the *p* values and *p* groups obtained by the two methods. Use a graphic to illustrate this. Why are they different? -->




<!-- end -->

<!-- start data wrangling -->

# (PART) STATISTICAL METHODS

# wrangling and summaries {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/wrangling.png" width="100%" style="display: block; margin: auto;" />

Data wrangling refers to the process of organizing, cleaning up, and making a "raw" data set more ready for downstream analysis. It is a key piece of any data analysis process. Here we will look at a few different aspects of wrangling, including data import, subsetting, pivoting, and summarizing data.

## {-}

## data import {-}

To analyze data that is stored on your own computer you can indeed import it into RStudio.

The easiest way to do this is to use the interactive command `readCSV()`, a function that comes with the phylochemistry source command. You run `readCSV()` in your console, then navigate to the data on your hard drive.

Another option is to read the data in from a path. For this, you will need to know the "path" to your data file. This is essentially the street address of your data on your computer's hard drive. Paths look different on Mac and PC.

* On Mac: `/Users/lucasbusta/Documents/sample_data_set.csv` (note the forward slashes!)
* On PC: `C:\\My Computer\\Documents\\sample_data_set.csv` (note double backward slashes!)

You can quickly find paths to files via the following:

* On Mac: Locate the file in Finder. Right-click on the file, hold the Option key, then click "Copy <file> as Pathname"
* On PC: Locate the file in Windows Explorer. Hold down the Shift key then right-click on the file. Click "Copy As Path"

With these paths, we can read in data using the `read_csv` command. We'll run `read_csv("<path_to_your_data>")`. Note the use of QUOTES `""`! Those are necessary. Also make sure your path uses the appropriate direction of slashes for your operating system.

## subsetting {-}

So far, we have always been passing whole data sets to ggplot to do our plotting. However, suppose we wanted to get at just certain portions of our dataset, say, specific columns, or specific rows? Here are a few ways to do that:


``` r
# To look at a single column (the third column)
head(alaska_lake_data[,3])
## # A tibble: 6 × 1
##   water_temp
##        <dbl>
## 1       6.46
## 2       6.46
## 3       6.46
## 4       6.46
## 5       6.46
## 6       6.46
```

``` r

# To look at select columns:
head(alaska_lake_data[,2:5])
## # A tibble: 6 × 4
##   park  water_temp    pH element
##   <chr>      <dbl> <dbl> <chr>  
## 1 BELA        6.46  7.69 C      
## 2 BELA        6.46  7.69 N      
## 3 BELA        6.46  7.69 P      
## 4 BELA        6.46  7.69 Cl     
## 5 BELA        6.46  7.69 S      
## 6 BELA        6.46  7.69 F
```

``` r

# To look at a single row (the second row)
head(alaska_lake_data[2,])
## # A tibble: 1 × 7
##   lake  park  water_temp    pH element mg_per_L element_type
##   <chr> <chr>      <dbl> <dbl> <chr>      <dbl> <chr>       
## 1 Devi… BELA        6.46  7.69 N          0.028 bound
```

``` r

# To look at select rows:
head(alaska_lake_data[2:5,])
## # A tibble: 4 × 7
##   lake  park  water_temp    pH element mg_per_L element_type
##   <chr> <chr>      <dbl> <dbl> <chr>      <dbl> <chr>       
## 1 Devi… BELA        6.46  7.69 N          0.028 bound       
## 2 Devi… BELA        6.46  7.69 P          0     bound       
## 3 Devi… BELA        6.46  7.69 Cl        10.4   free        
## 4 Devi… BELA        6.46  7.69 S          0.62  free
```

``` r

# To look at just a single column, by name
head(alaska_lake_data$pH)
## [1] 7.69 7.69 7.69 7.69 7.69 7.69
```

``` r

# To look at select columns by name
head(select(alaska_lake_data, park, water_temp))
## # A tibble: 6 × 2
##   park  water_temp
##   <chr>      <dbl>
## 1 BELA        6.46
## 2 BELA        6.46
## 3 BELA        6.46
## 4 BELA        6.46
## 5 BELA        6.46
## 6 BELA        6.46
```

## wide and long data {-}

When we make data tables by hand, it's often easy to make a **wide-style table** like the following. In it, the abundances of 7 different fatty acids in 10 different species are tabulated. Each fatty acid gets its own row, each species, its own column.


``` r
head(fadb_sample)
## # A tibble: 6 × 11
##   fatty_acid      Agonandra_brasiliensis Agonandra_silvatica
##   <chr>                            <dbl>               <dbl>
## 1 Hexadecanoic a…                    3.4                 1  
## 2 Octadecanoic a…                    6.2                 0.1
## 3 Eicosanoic acid                    4.7                 3.5
## 4 Docosanoic acid                   77.4                 0.4
## 5 Tetracosanoic …                    1.4                 1  
## 6 Hexacosanoic a…                    1.9                12.6
## # ℹ 8 more variables: Agonandra_excelsa <dbl>,
## #   Heisteria_silvianii <dbl>, Malania_oleifera <dbl>,
## #   Ximenia_americana <dbl>, Ongokea_gore <dbl>,
## #   Comandra_pallida <dbl>, Buckleya_distichophylla <dbl>,
## #   Nuytsia_floribunda <dbl>
```

While this format is very nice for filling in my hand (such as in a lab notebook or similar), it does not groove with ggplot and other `tidyverse` functions very well. We need to convert it into a **long-style table**. This is done using `pivot_longer()`. You can think of this function as transforming both your data's column names (or some of the column names) and your data matrix's values (in this case, the measurements) each into their own variables (i.e. columns). We can do this for our fatty acid dataset using the command below. In it, we specify what data we want to transform (`data = fadb_sample`), we need to tell it what columns we want to transform (`cols = 2:11`), what we want the new variable that contains column names to be called (`names_to = "plant_species"`) and what we want the new variable that contains matrix values to be called (`values_to = "relative_abundance"`). All together now:


``` r
pivot_longer(data = fadb_sample, cols = 2:11, names_to = "plant_species", values_to = "relative_abundance")
## # A tibble: 70 × 3
##    fatty_acid        plant_species        relative_abundance
##    <chr>             <chr>                             <dbl>
##  1 Hexadecanoic acid Agonandra_brasilien…                3.4
##  2 Hexadecanoic acid Agonandra_silvatica                 1  
##  3 Hexadecanoic acid Agonandra_excelsa                   1.2
##  4 Hexadecanoic acid Heisteria_silvianii                 2.9
##  5 Hexadecanoic acid Malania_oleifera                    0.7
##  6 Hexadecanoic acid Ximenia_americana                   3.3
##  7 Hexadecanoic acid Ongokea_gore                        1  
##  8 Hexadecanoic acid Comandra_pallida                    2.3
##  9 Hexadecanoic acid Buckleya_distichoph…                1.6
## 10 Hexadecanoic acid Nuytsia_floribunda                  3.8
## # ℹ 60 more rows
```

Brilliant! Now we have a long-style table that can be used with ggplot.

## the pipe (%>%) {-}

We have seen how to create new objects using `<-`, and we have been filtering and plotting data using, for example:


``` r
ggplot(filter(alaska_lake_data, park == "BELA"), aes(x = pH, y = lake)) + geom_col()
```

<img src="index_files/figure-html/unnamed-chunk-94-1.png" width="100%" style="display: block; margin: auto;" />

However, as our analyses get more complex, the code can get long and hard to read. We're going to use the pipe `%>%` to help us with this. Check it out:


``` r
alaska_lake_data %>%
  filter(park == "BELA") %>%
  ggplot(aes(x = pH, y = lake)) + geom_col()
```

<img src="index_files/figure-html/unnamed-chunk-95-1.png" width="100%" style="display: block; margin: auto;" />

Neat! Another way to think about the pipe:

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/pipe.jpg" width="100%" style="display: block; margin: auto;" />

The pipe will become more important as our analyses become more sophisticated, which happens very quickly when we start working with summary statistics, as we shall now see...

## summary statistics {-}

So far, we have been plotting raw data. This is well and good, but it is not always suitable. Often we have scientific questions that cannot be answered by looking at raw data alone, or sometimes there is too much raw data to plot. For this, we need summary statistics - things like averages, standard deviations, and so on. While these metrics can be computed in Excel, programming such can be time consuming, especially for group statistics. Consider the example below, which uses the `ny_trees` dataset. The NY Trees dataset contains information on nearly half a million trees in New York City (this is after considerable filtering and simplification):


``` r
head(ny_trees)
## # A tibble: 6 × 14
##   tree_height tree_diameter address        tree_loc pit_type
##         <dbl>         <dbl> <chr>          <chr>    <chr>   
## 1        21.1             6 1139 57 STREET Front    Sidewal…
## 2        59.0             6 2220 BERGEN A… Across   Sidewal…
## 3        92.4            13 2254 BERGEN A… Across   Sidewal…
## 4        50.2            15 2332 BERGEN A… Across   Sidewal…
## 5        95.0            21 2361 EAST   7… Front    Sidewal…
## 6        67.5            19 2409 EAST   7… Front    Continu…
## # ℹ 9 more variables: soil_lvl <chr>, status <chr>,
## #   spc_latin <chr>, spc_common <chr>, trunk_dmg <chr>,
## #   zipcode <dbl>, boroname <chr>, latitude <dbl>,
## #   longitude <dbl>
```

More than 300,000 observations of 14 variables! That's 4.2M data points! Now, what is the average and standard deviation of the height and diameter of each tree species within each NY borough? Do those values change for trees that are in parks versus sidewalk pits?? I don't even know how one would begin to approach such questions using traditional spreadsheets. Here, we will answer these questions with ease using two new commands: `group_by()` and `summarize()`. Let's get to it.

Say that we want to know (and of course, visualize) the mean and standard deviation of the heights of each tree species in NYC. We can see that data in first few columns of the NY trees dataset above, but how to calculate these statistics? In R, mean can be computed with `mean()` and standard deviation can be calculated with `sd()`. We will use the function `summarize()` to calculate summary statistics. So, we can calculate the average and standard deviation of all the trees in the data set as follows:


``` r
ny_trees %>%
  summarize(mean_height = mean(tree_height))
## # A tibble: 1 × 1
##   mean_height
##         <dbl>
## 1        72.6
```

``` r

ny_trees %>%
  summarize(stdev_height = sd(tree_height))
## # A tibble: 1 × 1
##   stdev_height
##          <dbl>
## 1         28.7
```

Great! But how to do this for each species? We need to subdivide the data by species, then compute the mean and standard deviation, then recombine the results into a new table. First, we use `group_by()`. Note that in ny_trees, species are indicated in the column called `spc_latin`. Once the data is grouped, we can use `summarize()` to compute statistics.


``` r
ny_trees %>%
  group_by(spc_latin) %>%
  summarize(mean_height = mean(tree_height))
## # A tibble: 12 × 2
##    spc_latin              mean_height
##    <chr>                        <dbl>
##  1 ACER PLATANOIDES              82.6
##  2 ACER RUBRUM                  106. 
##  3 ACER SACCHARINUM              65.6
##  4 FRAXINUS PENNSYLVANICA        60.6
##  5 GINKGO BILOBA                 90.4
##  6 GLEDITSIA TRIACANTHOS         53.0
##  7 PLATANUS ACERIFOLIA           82.0
##  8 PYRUS CALLERYANA              21.0
##  9 QUERCUS PALUSTRIS             65.5
## 10 QUERCUS RUBRA                111. 
## 11 TILIA CORDATA                 98.8
## 12 ZELKOVA SERRATA              101.
```

Bam. Mean height of each tree species. We can also count the number of observations using `n()`:


``` r
ny_trees %>%
  group_by(spc_latin) %>%
  summarize(number_of_individuals = n())
## # A tibble: 12 × 2
##    spc_latin              number_of_individuals
##    <chr>                                  <int>
##  1 ACER PLATANOIDES                       67260
##  2 ACER RUBRUM                            11506
##  3 ACER SACCHARINUM                       13161
##  4 FRAXINUS PENNSYLVANICA                 16987
##  5 GINKGO BILOBA                          15672
##  6 GLEDITSIA TRIACANTHOS                  48707
##  7 PLATANUS ACERIFOLIA                    80075
##  8 PYRUS CALLERYANA                       39125
##  9 QUERCUS PALUSTRIS                      37058
## 10 QUERCUS RUBRA                          10020
## 11 TILIA CORDATA                          25970
## 12 ZELKOVA SERRATA                        13221
```

Cool! `summarize()` is more powerful though, we can do many summary statistics at once:


``` r
ny_trees %>%
  group_by(spc_latin) %>%
  summarize(
    mean_height = mean(tree_height),
    stdev_height = sd(tree_height)
  ) -> ny_trees_by_spc_summ
ny_trees_by_spc_summ
## # A tibble: 12 × 3
##    spc_latin              mean_height stdev_height
##    <chr>                        <dbl>        <dbl>
##  1 ACER PLATANOIDES              82.6        17.6 
##  2 ACER RUBRUM                  106.         15.7 
##  3 ACER SACCHARINUM              65.6        16.6 
##  4 FRAXINUS PENNSYLVANICA        60.6        21.3 
##  5 GINKGO BILOBA                 90.4        24.5 
##  6 GLEDITSIA TRIACANTHOS         53.0        13.0 
##  7 PLATANUS ACERIFOLIA           82.0        16.0 
##  8 PYRUS CALLERYANA              21.0         5.00
##  9 QUERCUS PALUSTRIS             65.5         6.48
## 10 QUERCUS RUBRA                111.         20.7 
## 11 TILIA CORDATA                 98.8        32.6 
## 12 ZELKOVA SERRATA              101.         10.7
```

Now we can use this data in plotting. For this, we will use a new geom, `geom_pointrange`, which takes `x` and `y` aesthetics, as usual, but also requires two additional y-ish aesthetics `ymin` and `ymax` (or `xmin` and `xmax` if you want them to vary along x). Also, note that in the aesthetic mappings for `xmin` and `xmax`, we can use a mathematical expression: `mean-stdev` and `mean+stdev`, respectivey. In our case, these are `mean_height - stdev_height` and `mean_height + stdev_height`. Let's see it in action:


``` r
ny_trees_by_spc_summ %>%
ggplot() +
  geom_pointrange(
      aes(
        y = spc_latin,
        x = mean_height,
        xmin = mean_height - stdev_height,
        xmax = mean_height + stdev_height
      )
    )
```

<img src="index_files/figure-html/unnamed-chunk-102-1.png" width="100%" style="display: block; margin: auto;" />

Cool! Just like that, we've found (and visualized) the average and standard deviation of tree heights, by species, in NYC. But it doesn't stop there. We can use `group_by()` and `summarize()` on multiple variables (i.e. more groups). We can do this to examine the properties of each tree species in each NYC borough. Let's check it out:


``` r
ny_trees %>%
  group_by(spc_latin, boroname) %>%
  summarize(
    mean_diam = mean(tree_diameter),
    stdev_diam = sd(tree_diameter)
  ) -> ny_trees_by_spc_boro_summ
ny_trees_by_spc_boro_summ
## # A tibble: 48 × 4
## # Groups:   spc_latin [12]
##    spc_latin        boroname  mean_diam stdev_diam
##    <chr>            <chr>         <dbl>      <dbl>
##  1 ACER PLATANOIDES Bronx         13.9        6.74
##  2 ACER PLATANOIDES Brooklyn      15.4       14.9 
##  3 ACER PLATANOIDES Manhattan     11.6        8.45
##  4 ACER PLATANOIDES Queens        15.1       12.9 
##  5 ACER RUBRUM      Bronx         11.4        7.88
##  6 ACER RUBRUM      Brooklyn      10.5        7.41
##  7 ACER RUBRUM      Manhattan      6.63       4.23
##  8 ACER RUBRUM      Queens        14.1        8.36
##  9 ACER SACCHARINUM Bronx         19.7       10.5 
## 10 ACER SACCHARINUM Brooklyn      22.2       10.1 
## # ℹ 38 more rows
```

Now we have summary statistics for each tree species within each borough. This is different from the previous plot in that we now have an additional variable (boroname) in our summarized dataset. This additional variable needs to be encoded in our plot. Let's map boroname to x and facet over tree species, which used to be on x. We'll also manually modify the theme element `strip.text.y` to get the species names in a readable position.


``` r
ny_trees_by_spc_boro_summ %>%
ggplot() +
  geom_pointrange(
    aes(
      y = boroname,
      x = mean_diam,
      xmin = mean_diam-stdev_diam,
      xmax = mean_diam+stdev_diam
    )
  ) +
  facet_grid(spc_latin~.) +
  theme(
    strip.text.y = element_text(angle = 0)
  )
```

<img src="index_files/figure-html/unnamed-chunk-104-1.png" width="100%" style="display: block; margin: auto;" />

Excellent! And if we really want to go for something pretty:


``` r
ny_trees_by_spc_boro_summ %>%
ggplot() +
  geom_pointrange(
    aes(
      y = boroname,
      x = mean_diam,
      xmin = mean_diam-stdev_diam,
      xmax = mean_diam+stdev_diam,
      fill = spc_latin
    ), color = "black", shape = 21
  ) +
  labs(
    y = "Borough", 
    x = "Trunk diameter"
    # caption = str_wrap("Figure 1: Diameters of trees in New York City. Points correspond to average diameters of each tree species in each borough. Horizontal lines indicate the standard deviation of tree diameters. Points are colored according to tree species.", width = 80)
  ) +
  facet_grid(spc_latin~.) +
  guides(fill = "none") +
  scale_fill_brewer(palette = "Paired") +
  theme_bw() +
  theme(
    strip.text.y = element_text(angle = 0),
    plot.caption = element_text(hjust = 0.5)
  )
```

<img src="index_files/figure-html/unnamed-chunk-105-1.png" width="100%" style="display: block; margin: auto;" />

*Now* we are getting somewhere. It looks like there are some really big maple trees (Acer) in Queens.

## ordering {-}

We can also sort or order a data frame based on a specific column with the command `arrange()`. Let's have a quick look. Suppose we wanted to know which lake was the coldest:


``` r
arrange(alaska_lake_data, water_temp)
## # A tibble: 220 × 7
##    lake             park  water_temp    pH element mg_per_L
##    <chr>            <chr>      <dbl> <dbl> <chr>      <dbl>
##  1 Desperation_Lake NOAT        2.95  6.34 C          2.1  
##  2 Desperation_Lake NOAT        2.95  6.34 N          0.005
##  3 Desperation_Lake NOAT        2.95  6.34 P          0    
##  4 Desperation_Lake NOAT        2.95  6.34 Cl         0.2  
##  5 Desperation_Lake NOAT        2.95  6.34 S          2.73 
##  6 Desperation_Lake NOAT        2.95  6.34 F          0.01 
##  7 Desperation_Lake NOAT        2.95  6.34 Br         0    
##  8 Desperation_Lake NOAT        2.95  6.34 Na         1.11 
##  9 Desperation_Lake NOAT        2.95  6.34 K          0.16 
## 10 Desperation_Lake NOAT        2.95  6.34 Ca         5.87 
## # ℹ 210 more rows
## # ℹ 1 more variable: element_type <chr>
```

Or suppose we wanted to know which was the warmest?


``` r
arrange(alaska_lake_data, desc(water_temp))
## # A tibble: 220 × 7
##    lake      park  water_temp    pH element mg_per_L
##    <chr>     <chr>      <dbl> <dbl> <chr>      <dbl>
##  1 Lava_Lake BELA        20.2  7.42 C          8.3  
##  2 Lava_Lake BELA        20.2  7.42 N          0.017
##  3 Lava_Lake BELA        20.2  7.42 P          0.001
##  4 Lava_Lake BELA        20.2  7.42 Cl         2.53 
##  5 Lava_Lake BELA        20.2  7.42 S          0.59 
##  6 Lava_Lake BELA        20.2  7.42 F          0.04 
##  7 Lava_Lake BELA        20.2  7.42 Br         0.01 
##  8 Lava_Lake BELA        20.2  7.42 Na         2.93 
##  9 Lava_Lake BELA        20.2  7.42 K          0.57 
## 10 Lava_Lake BELA        20.2  7.42 Ca        11.8  
## # ℹ 210 more rows
## # ℹ 1 more variable: element_type <chr>
```

`arrange()` will work on grouped data, which is particularly useful in combination with `slice()`, which can show us the first n elements in each group:


``` r
alaska_lake_data %>%
  group_by(park) %>%
  arrange(water_temp) %>%
  slice(1)
## # A tibble: 3 × 7
## # Groups:   park [3]
##   lake  park  water_temp    pH element mg_per_L element_type
##   <chr> <chr>      <dbl> <dbl> <chr>      <dbl> <chr>       
## 1 Devi… BELA        6.46  7.69 C            3.4 bound       
## 2 Wild… GAAR        5.5   6.98 C            6.5 bound       
## 3 Desp… NOAT        2.95  6.34 C            2.1 bound
```

It looks like the coldest lakes in the three parks are Devil Mountain Lake, Wild Lake, and Desperation Lake!

## mutate {-}

One last thing before our exercises... there is another command called `mutate()`. It is like summarize it calculates user-defined statistics, but it creates output on a per-observation level instead of for each group. This means that it doesn't make the data set smaller, in fact it makes it bigger, by creating a new row for the new variables defined inside `mutate()`. It can also take grouped data. This is really useful for calculating percentages within groups. For example: within each park, what percent of the park's total dissolved sulfur does each lake have?


``` r
alaska_lake_data %>%
  filter(element == "S") %>%
  group_by(park) %>%
  select(lake, park, element, mg_per_L) %>%
  mutate(percent_S = mg_per_L/sum(mg_per_L)*100)
## # A tibble: 20 × 5
## # Groups:   park [3]
##    lake                park  element mg_per_L percent_S
##    <chr>               <chr> <chr>      <dbl>     <dbl>
##  1 Devil_Mountain_Lake BELA  S           0.62     32.3 
##  2 Imuruk_Lake         BELA  S           0.2      10.4 
##  3 Kuzitrin_Lake       BELA  S           0.29     15.1 
##  4 Lava_Lake           BELA  S           0.59     30.7 
##  5 North_Killeak_Lake  BELA  S           0.04      2.08
##  6 White_Fish_Lake     BELA  S           0.18      9.38
##  7 Iniakuk_Lake        GAAR  S          12.1      13.2 
##  8 Kurupa_Lake         GAAR  S          12.4      13.6 
##  9 Lake_Matcharak      GAAR  S          13.3      14.5 
## 10 Lake_Selby          GAAR  S           7.92      8.64
## 11 Nutavukti_Lake      GAAR  S           2.72      2.97
## 12 Summit_Lake         GAAR  S           3.21      3.50
## 13 Takahula_Lake       GAAR  S           5.53      6.03
## 14 Walker_Lake         GAAR  S           5.77      6.30
## 15 Wild_Lake           GAAR  S          28.7      31.3 
## 16 Desperation_Lake    NOAT  S           2.73     25.8 
## 17 Feniak_Lake         NOAT  S           4.93     46.5 
## 18 Lake_Kangilipak     NOAT  S           0.55      5.19
## 19 Lake_Narvakrak      NOAT  S           1.38     13.0 
## 20 Okoklik_Lake        NOAT  S           1.01      9.53
```

The percent columns for each park add to 100%, so, for example, Devil Mountain Lake has 32.3% of BELA's dissolved sulfur.

<!-- ## exercises {-}

Isn’t seven the most powerfully magical number? *Isn’t seven the most powerfully magical number?*  Yes... I think the idea of a seven-part assignment would greatly appeal to an alchemist.

In this set of exercises we are going to use the periodic table. After you run source() you can load that data set using `periodic_table`. Please use that dataset to run analyses and answer the following questions/prompts. Compile the answers in an R Markdown document, compile it as a pdf, and upload it to the Canvas assignment. Please let me know if you have any questions. Good luck, and have fun!

Some pointers:

- If your code goes off the page, please wrap it across multiple lines, as shown in some of the examples in the previous set of exercises.

- Don't be afraid to put the variable with the long elements / long text on the y-axis and the continuous variable on the x-axis.

- If your axis tick labels are overlapping or not visible, do something to fix that. Some solutions could be: move the legend to the top of the plot (`theme(legend.position = "top")`), rotate the text (`theme(axis.text.x = element_text(angle = 90)`), or make the text smaller (`theme(axis.text.x = element_text(size = 8)`).

1. Make a plot using `geom_point()` that shows the average atomic weight of the elements discovered in each year spanned by the dataset (i.e. what was the average weight of the elements discovered in 1900? 1901? 1902? etc.). You should see a trend, particularly after 1950. What do you think has caused this trend?



2. The column `state_at_RT` indicates the state of each element at room temperate. Make a plot that shows the average first ionization potential of all the elements belonging to each state group indicated in `state_at_RT` (i.e. what is the average 1st ionization potential of all elements that are solid at room temp? liquid? etc.). Which is the highest?



3. Filter the dataset so that only elements with atomic number less than 85 are included. Considering only these elements, what is the average and standard deviation of boiling points for each type of `crystal_structure`? Make a plot using `geom_pointrange()` that shows the mean and standard deviation of each of these groups. What's up with elements that have a cubic crystal structure?



4. Now filter the original dataset so that only elements with atomic number less than 37 are considered. The elements in this dataset belong to the first four periods. What is the average abundance of each of these four *periods* in seawater? i.e. what is the average abundance of all elements from period 1? period 2? etc. Which period is the most abundant? In this context what does "CHON" mean? (not the rock band, though they are also excellent, especially that song that features GoYama)



5. Now filter the original dataset so that only elements with atomic number less than 103 are considered. Filter it further so that elements from group number 18 are excluded. Using this twice-filtered dataset, compute the average, minimum, and maximum values for electronegativiy for each `group_number`. Use `geom_point()` and `geom_errorbar()` to illustrate the average, minimum, and maximum values for each group number.



6. Filter the dataset so that only elements with atomic number less than 85 are considered. Group these by `color`. Now filter out those that have `color == "colorless"`. Of the remaining elements, which has the widest range of specific heats? Use `geom_point()` and `geom_errorbar()` to illustrate the mean and standard deviation of each color's specific heats.



7. You have learned many things in this course so far. `read_csv()`, `filter()`, `ggplot()`, and now `group_by()`, `summarize()`, `mutate()`, `arrange()`, and `slice()`. Using **all** these commands, create one or more graphics to illustrate what you consider to be one or more interesting trends in a data set of your own choosing. Use theme elements and scales to enhance your plot. Give your plot a nice caption based on the caption guide in this book. -->

## {-}

## further reading {-}

Be sure to check out the Tidy Data Tutor: https://tidydatatutor.com/vis.html. An easy way to visualize what's going on during data wrangling!

<!-- end -->

<!-- start dimensionality reduction -->

# dimensional reduction {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/dimensionality.png" width="100%" style="display: block; margin: auto;" />

In the previous chapters, we looked at how to explore our data sets by visualizing many variables and manually identifying trends. Sometimes, we encounter data sets with so many variables, that it is not reasonable to manually select certain variables with which to create plots and manually search for trends. In these cases, we need dimensionality reduction - a set of techniques that helps us identify which variables are driving differences among our samples. In this course, we will conduct dimensionality reduction useing `runMatrixAnalysis()`, a function that is loaded into your R Session when you run the source() command.

Matrix analyses can be a bit tricky to set up. There are two things that we can do to help us with this: (i) we will use a template for `runMatrixAnalysis()` (see below) and (ii) it is *critical* that we think about our data in terms of **samples** and **analytes**. Let's consider our Alaska lakes data set:


``` r
alaska_lake_data
## # A tibble: 220 × 7
##    lake              park  water_temp    pH element mg_per_L
##    <chr>             <chr>      <dbl> <dbl> <chr>      <dbl>
##  1 Devil_Mountain_L… BELA        6.46  7.69 C          3.4  
##  2 Devil_Mountain_L… BELA        6.46  7.69 N          0.028
##  3 Devil_Mountain_L… BELA        6.46  7.69 P          0    
##  4 Devil_Mountain_L… BELA        6.46  7.69 Cl        10.4  
##  5 Devil_Mountain_L… BELA        6.46  7.69 S          0.62 
##  6 Devil_Mountain_L… BELA        6.46  7.69 F          0.04 
##  7 Devil_Mountain_L… BELA        6.46  7.69 Br         0.02 
##  8 Devil_Mountain_L… BELA        6.46  7.69 Na         8.92 
##  9 Devil_Mountain_L… BELA        6.46  7.69 K          1.2  
## 10 Devil_Mountain_L… BELA        6.46  7.69 Ca         5.73 
## # ℹ 210 more rows
## # ℹ 1 more variable: element_type <chr>
```

We can see that this dataset is comprised of measurements of various *analytes* (i.e. several chemical elements, as well as water_temp, and pH), in different *samples* (i.e. lakes). We need to tell the `runMatrixAnalysis()` function how each column relates to this samples and analytes structure. See the image below for an explanation.

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/runMatrixAnalysis1.png" width="100%" style="display: block; margin: auto;" />

## {-}

## pca {-}

"Which analytes are driving differences among my samples?"
"Which analytes in my data set are correlated?"

### theory {-}

PCA looks at all the variance in a high dimensional data set and chooses new axes within that data set that align with the directions containing highest variance. These new axes are called principal components. Let's look at an example:

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/PCA.png" width="100%" style="display: block; margin: auto;" />

In the example above, the three dimensional space can be reduced to a two dimensional space with the principal components analysis. New axes (principal components) are selected (bold arrows on left) that become the x and y axes in the principal components space (right).

We can run and visualize principal components analyses using the `runMatrixAnalysis()` function as in the example below. As you can see in the output, the command provides the sample_IDs, sample information, then the coordinates for each sample in the 2D projection (the "PCA plot") and the raw data, in case you wish to do further processing.


``` r
AK_lakes_pca <- runMatrixAnalysis(
  data = alaska_lake_data,
  analysis = c("pca"),
  column_w_names_of_multiple_analytes = "element",
  column_w_values_for_multiple_analytes = "mg_per_L",
  columns_w_values_for_single_analyte = c("water_temp", "pH"),
  columns_w_additional_analyte_info = "element_type",
  columns_w_sample_ID_info = c("lake", "park"),
  scale_variance = TRUE
)
head(AK_lakes_pca)
## # A tibble: 6 × 18
##   sample_unique_ID      lake  park   Dim.1  Dim.2 water_temp
##   <chr>                 <chr> <chr>  <dbl>  <dbl>      <dbl>
## 1 Devil_Mountain_Lake_… Devi… BELA   0.229 -0.861       6.46
## 2 Imuruk_Lake_BELA      Imur… BELA  -1.17  -1.62       17.4 
## 3 Kuzitrin_Lake_BELA    Kuzi… BELA  -0.918 -1.15        8.06
## 4 Lava_Lake_BELA        Lava… BELA   0.219 -1.60       20.2 
## 5 North_Killeak_Lake_B… Nort… BELA   9.46   0.450      11.3 
## 6 White_Fish_Lake_BELA  Whit… BELA   4.17  -0.972      12.0 
## # ℹ 12 more variables: pH <dbl>, C <dbl>, N <dbl>, P <dbl>,
## #   Cl <dbl>, S <dbl>, F <dbl>, Br <dbl>, Na <dbl>,
## #   K <dbl>, Ca <dbl>, Mg <dbl>
```

Let's plot the 2D projection of the Alaska lakes data:


``` r
ggplot(data = AK_lakes_pca, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(fill = park), shape = 21, size = 4, alpha = 0.8) +
  geom_label_repel(aes(label = lake), alpha = 0.5) +
  theme_classic()
```

<img src="index_files/figure-html/unnamed-chunk-121-1.png" width="100%" style="display: block; margin: auto;" />

Great! In this plot we can see that White Fish Lake and North Killeak Lake, both in BELA park, are quite different from the other parks (they are separated from the others along dimension 1, i.e. the first principal component). At the same time, Wild Lake, Iniakuk Lake, Walker Lake, and several other lakes in GAAR park are different from all the others (they are separated from the others along dimension 2, i.e. the second principal component).

Important question: what makes the lakes listed above different from the others? Certainly some aspect of their chemistry, since that's the data that this analysis is built upon, but how do we determine which analyte(s) are driving the differences among the lakes that we see in the PCA plot?

### ordination plots {-}

Let's look at how to access the information about which analytes are major contributors to each principal component. This is important because it will tell you which analytes are associated with particular dimensions, and by extension, which analytes are associated with (and are markers for) particular groups in the PCA plot. This can be determined using an ordination plot. Let's look at an example. We can obtain the ordination plot information using `runMatrixAnalysis()` with `analysis = "pca_ord"`:


```
## # A tibble: 6 × 3
##   analyte      Dim.1   Dim.2
##   <chr>        <dbl>   <dbl>
## 1 water_temp 0.0750  -0.261 
## 2 pH         0.686    0.0185
## 3 C          0.290   -0.242 
## 4 N          0.00435  0.714 
## 5 P          0.473   -0.0796
## 6 Cl         0.953    0.0148
```

We can now visualize the ordination plot using our standard ggplot plotting techniques. Note the use of `geom_label_repel()` and `filter()` to label certain segments in the ordination plot. You do not need to use `geom_label_repel()`, you could use the built in `geom_label()`, but `geom_label_repel()` can make labelling your segments easier.


``` r
AK_lakes_pca_ord <- runMatrixAnalysis(
  data = alaska_lake_data,
  analysis = c("pca_ord"),
  column_w_names_of_multiple_analytes = "element",
  column_w_values_for_multiple_analytes = "mg_per_L",
  columns_w_values_for_single_analyte = c("water_temp", "pH"),
  columns_w_additional_analyte_info = "element_type",
  columns_w_sample_ID_info = c("lake", "park")
)
head(AK_lakes_pca_ord)
## # A tibble: 6 × 3
##   analyte      Dim.1   Dim.2
##   <chr>        <dbl>   <dbl>
## 1 water_temp 0.0750  -0.261 
## 2 pH         0.686    0.0185
## 3 C          0.290   -0.242 
## 4 N          0.00435  0.714 
## 5 P          0.473   -0.0796
## 6 Cl         0.953    0.0148
```

``` r

ggplot(AK_lakes_pca_ord) +
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, color = analyte), size = 1) +
  geom_circle(aes(x0 = 0, y0 = 0, r = 1)) +
  geom_label_repel(
    data = filter(AK_lakes_pca_ord, Dim.1 > 0.9, Dim.2 < 0.1, Dim.2 > -0.1),
    aes(x = Dim.1, y = Dim.2, label = analyte), xlim = c(1,1.5)
  ) +
  geom_label_repel(
    data = filter(AK_lakes_pca_ord, Dim.2 > 0.5),
    aes(x = Dim.1, y = Dim.2, label = analyte), direction = "y", ylim = c(1,1.5)
  ) +
  coord_cartesian(xlim = c(-1,1.5), ylim = c(-1,1.5)) +
  theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-123-1.png" width="100%" style="display: block; margin: auto;" />

Great! Here is how to read the ordination plot:

1. When considering one analyte's vector: the vector's projected value on an axis shows how much its variance is aligned with that principal component.

2. When considering two analyte vectors: the angle between two vectors indicates how correlated those two variables are. If they point in the same direction, they are highly correlated. If they meet each other at 90 degrees, they are not very correlated. If they meet at ~180 degrees, they are negatively correlated. If say that one analyte is "1.9" with respect to dimension 2 and another is "-1.9" with respect to dimension 2. Let's also say that these vectors are ~"0" with respect to dimension 1.

With the ordination plot above, we can now see that the abundances of K, Cl, Br, and Na are the major contributors of variance to the first principal component (or the first dimension). The abundances of these elements are what make White Fish Lake and North Killeak Lake different from the other lakes. We can also see that the abundances of N, S, and Ca are the major contributors to variance in the second dimension, which means that these elements ar what set Wild Lake, Iniakuk Lake, Walker Lake, and several other lakes in GAAR park apart from the rest of the lakes in the data set. It slightly easier to understand this if we look at an overlay of the two plots, which is often called a "biplot":


``` r
AK_lakes_pca <- runMatrixAnalysis(
  data = alaska_lake_data,
  analysis = c("pca"),
  column_w_names_of_multiple_analytes = "element",
  column_w_values_for_multiple_analytes = "mg_per_L",
  columns_w_values_for_single_analyte = c("water_temp", "pH"),
  columns_w_additional_analyte_info = "element_type",
  columns_w_sample_ID_info = c("lake", "park"),
  scale_variance = TRUE
)
## Replacing NAs in your data with mean
```

``` r

AK_lakes_pca_ord <- runMatrixAnalysis(
  data = alaska_lake_data,
  analysis = c("pca_ord"),
  column_w_names_of_multiple_analytes = "element",
  column_w_values_for_multiple_analytes = "mg_per_L",
  columns_w_values_for_single_analyte = c("water_temp", "pH"),
  columns_w_additional_analyte_info = "element_type",
  columns_w_sample_ID_info = c("lake", "park")
)
## Replacing NAs in your data with mean
```

``` r

ggplot() +
  geom_point(
    data = AK_lakes_pca, 
    aes(x = Dim.1, y = Dim.2, fill = park), shape = 21, size = 4, alpha = 0.8
  ) +
  # geom_label_repel(aes(label = lake), alpha = 0.5) +
  geom_segment(
    data = AK_lakes_pca_ord,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2, color = analyte),
    size = 1
  ) +
  scale_color_manual(values = discrete_palette) +
  theme_classic()
```

<img src="index_files/figure-html/unnamed-chunk-124-1.png" width="100%" style="display: block; margin: auto;" />

Note that you do not have to plot ordination data as a circular layout of segments. Sometimes it is much easier to plot (and interpret!) alternatives:


``` r
AK_lakes_pca_ord %>%
  ggplot(aes(x = Dim.1, y = analyte)) +
    geom_point(aes(fill = analyte), shape = 22, size = 3) +
    scale_fill_manual(values = discrete_palette) +
    theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-125-1.png" width="100%" style="display: block; margin: auto;" />

### principal components {-}

We also can access information about the how much of the variance in the data set is explained by each principal component, and we can plot that using ggplot:


``` r
AK_lakes_pca_dim <- runMatrixAnalysis(
  data = alaska_lake_data,
  analysis = c("pca_dim"),
  column_w_names_of_multiple_analytes = "element",
  column_w_values_for_multiple_analytes = "mg_per_L",
  columns_w_values_for_single_analyte = c("water_temp", "pH"),
  columns_w_additional_analyte_info = "element_type",
  columns_w_sample_ID_info = c("lake", "park")
)
head(AK_lakes_pca_dim)
## # A tibble: 6 × 2
##   principal_component percent_variance_explained
##                 <dbl>                      <dbl>
## 1                   1                      48.8 
## 2                   2                      18.6 
## 3                   3                      11.6 
## 4                   4                       7.88
## 5                   5                       4.68
## 6                   6                       3.33
```

``` r

ggplot(
  data = AK_lakes_pca_dim, 
  aes(x = principal_component, y = percent_variance_explained)
) +
  geom_line() +
  geom_point() +
  theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-126-1.png" width="100%" style="display: block; margin: auto;" />

Cool! We can see that the first principal component retains nearly 50% of the variance in the original dataset, while the second dimension contains only about 20%. We can derive an important notion about PCA visualization from this: the scales on the two axes need to be the same for distances between points in the x and y directions to be comparable. This can be accomplished using `coord_fixed()` as an addition to your ggplots.

<!-- ### exercises {-}

In this set of exercises, as you are filling out the `runMatrixAnalysis()` template, you can use the `colnames()` function to help you specify a long list of column names rather than typing them out by hand. For example, in the periodic table data set, we can refer to a set of columns (columns 10 through 20) with the following command:


``` r
colnames(periodic_table_subset)[10:20]
##  [1] "electronegativity_pauling"         
##  [2] "first_ionization_poten_eV"         
##  [3] "second_ionization_poten_eV"        
##  [4] "third_ionization_poten_eV"         
##  [5] "electron_affinity_eV"              
##  [6] "atomic_radius_ang"                 
##  [7] "ionic_radius_ang"                  
##  [8] "covalent_radius_ang"               
##  [9] "atomic_volume_cm3_per_mol"         
## [10] "electrical_conductivity_mho_per_cm"
## [11] "specific_heat_J_per_g_K"
```

``` r
colnames(periodic_table_subset)[c(18:20, 23:25)]
## [1] "atomic_volume_cm3_per_mol"         
## [2] "electrical_conductivity_mho_per_cm"
## [3] "specific_heat_J_per_g_K"           
## [4] "thermal_conductivity_W_per_m_K"    
## [5] "polarizability_A_cubed"            
## [6] "heat_atomization_kJ_per_mol"
```
We can use that command in the template, as in the example below. With the notation `colnames(periodic_table_subset)[c(5:7,9:25)]`, we can mark columns 5 - 7 and 9 - 25 as columns_w_values_for_single_analyte (note what happens when you run `c(5:7,9:25)` in the console, and what happens when you run `colnames(periodic_table_subset)[c(5:7,9:25)]` in the console). With the notation `colnames(periodic_table_subset)[c(1:4, 8)]` we can mark columns 1 - 4 and column 8 as columns_w_sample_ID_info (note what happens when you run `c(1:4, 8)` in the console, and what happens when you run `colnames(periodic_table_subset)[c(1:4, 8)]` in the console).

#### human metabolomics {-}

For these questions, work with a dataset describing metabolomics data (i.e. abundances of > 100 different biochemicals) from each of 93 human patients, some of which have Chronic Kidney Disease. Your task is to discover a biomarker for Chronic Kidney Disease. This means that you will need to determine a metabolite whose abundance is strongly associated with the disease. To do this you should complete the following:

1. Run a PCA analysis on `metabolomics_data` (i.e. `runMatrixAnalysis()` with `analysis = "pca"`)
2. Plot the results of the analysis to determine which principal component (i.e. dimension) separates the healthy and kidney_disease samples.
3. Obtain the ordination plot coordinates for the analytes in the PCA analysis (i.e. `runMatrixAnalysis()` with `analysis = "pca_ord"`). In your own words, how does this plot correspond to the original data set?
4. Visualize the ordination plot and determine which of the analytes are strongly associated with the principal component (i.e. dimension) separates the healthy and kidney_disease samples.
5. Bingo! These analytes are associated with Chronic Kidney Disease and could be biomarkers for such.
6. Complete this PCA analysis by creating a scree plot (i.e. use `analysis = "pca_dim"`). In your own words, what does this plot mean?



#### grape vine chemistry {-}

For this set of quesions, work with a dataset describing metabolomics data (i.e. abundances of > 100 different biochemicals) from 5 different wine grape varieties. Your task is to discover a biomarker for Chardonnay and a biomarker for Cabernet Sauvignon. This means that you will need to identify two metabolites, each of which are associated with one of those two grape varieties. To do this you should complete the following:

1. Run a PCA analysis on `wine_grape_data` (i.e. `runMatrixAnalysis()` with `analysis = "pca"`)
2. Plot the results of the analysis to determine which principal component (i.e. dimension) separates the Chardonnay samples from the other varieties and the Cabernet Sauvignon samples from the other varieties.
3. Obtain the ordination plot coordinates for the analytes in the PCA analysis (i.e. `runMatrixAnalysis()` with `analysis = "pca_ord"`). In your own words, how does this plot correspond to the original data set?
4. Visualize the ordination plot and determine which of the analytes are strongly associated with the principal component (i.e. dimension) separates the Chardonnay samples from the other varieties and the Cabernet Sauvignon samples from the other varieties.
5. Bingo! These analytes are associated with those varieites and could be biomarkers for such.
6. Complete this PCA analysis by creating a scree plot (i.e. use `analysis = "pca_dim"`). In your own words, what does this plot mean?


 -->


## tsne and umap {-}


``` r
set.seed(235)
runMatrixAnalysis(
  data = metabolomics_data,
  analysis = "pca",
  column_w_names_of_multiple_analytes = NULL,
  column_w_values_for_multiple_analytes = NULL,
  columns_w_values_for_single_analyte = colnames(metabolomics_data)[c(3:126)],
  columns_w_additional_analyte_info = NULL,
  columns_w_sample_ID_info = colnames(metabolomics_data)[c(1:2)],
  na_replacement = "mean"
) -> pca_data
## Replacing NAs in your data with mean
```

``` r
pca_data$technique <- "pca_data"
colnames(pca_data) <- gsub("\\.", "_", colnames(pca_data))
pca_data$Dim_1 <- as.numeric(scale(pca_data$Dim_1))
pca_data$Dim_2 <- as.numeric(scale(pca_data$Dim_2))


runMatrixAnalysis(
  data = metabolomics_data,
  analysis = "umap",
  column_w_names_of_multiple_analytes = NULL,
  column_w_values_for_multiple_analytes = NULL,
  columns_w_values_for_single_analyte = colnames(metabolomics_data)[c(3:126)],
  columns_w_additional_analyte_info = NULL,
  columns_w_sample_ID_info = colnames(metabolomics_data)[c(1:2)],
  na_replacement = "mean"
) -> umap_data
## Replacing NAs in your data with mean
```

``` r
umap_data$technique <- "umap_data"
umap_data$Dim_1 <- as.numeric(scale(umap_data$Dim_1))
umap_data$Dim_2 <- as.numeric(scale(umap_data$Dim_2))


runMatrixAnalysis(
  data = metabolomics_data,
  analysis = "tsne",
  column_w_names_of_multiple_analytes = NULL,
  column_w_values_for_multiple_analytes = NULL,
  columns_w_values_for_single_analyte = colnames(metabolomics_data)[c(3:126)],
  columns_w_additional_analyte_info = NULL,
  columns_w_sample_ID_info = colnames(metabolomics_data)[c(1:2)],
  na_replacement = "mean"
) -> tsne_data
## Replacing NAs in your data with mean
```

``` r
tsne_data$technique <- "tsne_data"
tsne_data$Dim_1 <- as.numeric(scale(tsne_data$Dim_1))
tsne_data$Dim_2 <- as.numeric(scale(tsne_data$Dim_2))


data <- rbind(pca_data, umap_data, tsne_data)

p1 <- ggplot(data) +
  geom_point(aes(x = Dim_1, y = Dim_2, fill = patient_status), shape = 21, size= 4) +
  facet_grid(technique~., scales = "free") +
  scale_fill_brewer(palette = "Set1")

p2 <- ggplot(data) +
  geom_point(aes(x = Dim_1, y = Dim_2, fill = patient_status), shape = 21, size= 4) +
  facet_grid(technique~., scales = "free") +
  scale_fill_brewer(palette = "Set1")

p1 + p2
```

<img src="index_files/figure-html/unnamed-chunk-130-1.png" width="100%" style="display: block; margin: auto;" />

## {-}

## further reading {-}

Here is a video that nicely explains PCA: https://www.youtube.com/watch?v=FgakZw6K1QQ&list=PLblh5JKOoLUIcdlgu78MnlATeyx4cEVeR

https://datavizpyr.com/how-to-make-umap-plot-in-r/

https://datavizpyr.com/how-to-make-tsne-plot-in-r/

https://pair-code.github.io/understanding-umap/

https://www.youtube.com/watch?v=jth4kEvJ3P8

<!-- end -->

<!-- start clustering -->

# flat clustering {-}

"Do my samples fall into definable clusters?"

## {-}

## kmeans {-}

One of the questions we've been asking is "which of my samples are most closely related?". We've been answering that question using clustering. However, now that we know how to run principal components analyses, we can use another approach. This alternative approach is called k-means, and can help us decide how to assign our data into clusters. It is generally desirable to have a small number of clusters, however, this must be balanced by not having the variance within each cluster be too big. To strike this balance point, the elbow method is used. For it, we must first determine the maximum within-group variance at each possible number of clusters. An illustration of this is shown in **A** below:

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/kmeans.png" width="100%" style="display: block; margin: auto;" />

One we know within-group variances, we find the "elbow" point - the point with minimum angle theta - thus picking the outcome with a good balance of cluster number and within-cluster variance (illustrated above in **B** and **C**.)

Let's try k-means using `runMatrixAnalysis`. For this example, let's run it on the PCA projection of the alaska lakes data set. We can set `analysis = "kmeans"`. When we do this, an application will load that will show us the threshold value for the number of clusters we want. We set the number of clusters and then close the app. In the context of markdown document, simply provide the number of clusters to the `parameters` argument:



``` r
alaska_lake_data_pca <- runMatrixAnalysis(
    data = alaska_lake_data,
    analysis = c("pca"),
    column_w_names_of_multiple_analytes = "element",
    column_w_values_for_multiple_analytes = "mg_per_L",
    columns_w_values_for_single_analyte = c("water_temp", "pH"),
    columns_w_additional_analyte_info = "element_type",
    columns_w_sample_ID_info = c("lake", "park")
)

alaska_lake_data_pca_clusters <- runMatrixAnalysis(
    data = alaska_lake_data_pca,
    analysis = c("kmeans"),
    parameters = c(5),
    column_w_names_of_multiple_analytes = NULL,
    column_w_values_for_multiple_analytes = NULL,
    columns_w_values_for_single_analyte = c("Dim.1", "Dim.2"),
    columns_w_sample_ID_info = "sample_unique_ID"
)

alaska_lake_data_pca_clusters <- left_join(alaska_lake_data_pca_clusters, alaska_lake_data_pca) 
```

We can plot the results and color them according to the group that kmeans suggested. We can also highlight groups using `geom_mark_ellipse`. Note that it is recommended to specify both `fill` and `label` for geom_mark_ellipse:


``` r
alaska_lake_data_pca_clusters$cluster <- factor(alaska_lake_data_pca_clusters$cluster)
ggplot() +
  geom_point(
    data = alaska_lake_data_pca_clusters,
    aes(x = Dim.1, y = Dim.2, fill = cluster), shape = 21, size = 5, alpha = 0.6
  ) +
  geom_mark_ellipse(
    data = alaska_lake_data_pca_clusters,
    aes(x = Dim.1, y = Dim.2, label = cluster, fill = cluster), alpha = 0.2
  ) +
  theme_classic() +
  coord_cartesian(xlim = c(-7,12), ylim = c(-4,5)) +
  scale_fill_manual(values = discrete_palette) 
```

<img src="index_files/figure-html/unnamed-chunk-133-1.png" width="100%" style="display: block; margin: auto;" />

## dbscan {-}

There is another method to define clusters that we call dbscan. In this method, not all points are necessarily assigned to a cluster, and we define clusters according to a set of parameters, instead of simply defining the number of clusteres, as in kmeans. In interactive mode, `runMatrixAnalysis()` will again load an interactive means of selecting parameters for defining dbscan clusters ("k", and "threshold"). In the context of markdown document, simply provide "k" and "threshold" to the `parameters` argument:


``` r
alaska_lake_data_pca <- runMatrixAnalysis(
    data = alaska_lake_data,
    analysis = c("pca"),
    column_w_names_of_multiple_analytes = "element",
    column_w_values_for_multiple_analytes = "mg_per_L",
    columns_w_values_for_single_analyte = c("water_temp", "pH"),
    columns_w_additional_analyte_info = "element_type",
    columns_w_sample_ID_info = c("lake", "park")
) 

alaska_lake_data_pca_clusters <- runMatrixAnalysis(
    data = alaska_lake_data_pca,
    analysis = c("dbscan"),
    parameters = c(4, 0.45),
    column_w_names_of_multiple_analytes = NULL,
    column_w_values_for_multiple_analytes = NULL,
    columns_w_values_for_single_analyte = c("Dim.1", "Dim.2"),
    columns_w_sample_ID_info = "sample_unique_ID"
)

alaska_lake_data_pca_clusters <- left_join(alaska_lake_data_pca_clusters, alaska_lake_data_pca)
```
We can make the plot in the same way, but please note that to get `geom_mark_ellipse` to omit the ellipse for NAs you need to feed it data without NAs:


``` r
alaska_lake_data_pca_clusters$cluster <- factor(alaska_lake_data_pca_clusters$cluster)
ggplot() +
  geom_point(
    data = alaska_lake_data_pca_clusters,
    aes(x = Dim.1, y = Dim.2, fill = cluster), shape = 21, size = 5, alpha = 0.6
  ) +
  geom_mark_ellipse(
    data = drop_na(alaska_lake_data_pca_clusters),
    aes(x = Dim.1, y = Dim.2, label = cluster, fill = cluster), alpha = 0.2
  ) +
  theme_classic() +
  coord_cartesian(xlim = c(-7,12), ylim = c(-4,5)) +
  scale_fill_manual(values = discrete_palette) 
```

<img src="index_files/figure-html/unnamed-chunk-135-1.png" width="100%" style="display: block; margin: auto;" />

## summarize by cluster {-}

One more important point: when using kmeans or dbscan, we can use the clusters as groupings for summary statistics. For example, suppose we want to see the differences in abundances of certain chemicals among the clusters:


``` r
alaska_lake_data_pca <- runMatrixAnalysis(
  data = alaska_lake_data,
  analysis = c("pca"),
  column_w_names_of_multiple_analytes = "element",
  column_w_values_for_multiple_analytes = "mg_per_L",
  columns_w_values_for_single_analyte = c("water_temp", "pH"),
  columns_w_additional_analyte_info = "element_type",
  columns_w_sample_ID_info = c("lake", "park")
)

alaska_lake_data_pca_clusters <- runMatrixAnalysis(
  data = alaska_lake_data_pca,
  analysis = c("dbscan"),
  parameters = c(4, 0.45),
  column_w_names_of_multiple_analytes = NULL,
  column_w_values_for_multiple_analytes = NULL,
  columns_w_values_for_single_analyte = c("Dim.1", "Dim.2"),
  columns_w_sample_ID_info = "sample_unique_ID",
  columns_w_additional_analyte_info = colnames(alaska_lake_data_pca)[6:18]
) 

alaska_lake_data_pca_clusters <- left_join(alaska_lake_data_pca_clusters, alaska_lake_data_pca)

alaska_lake_data_pca_clusters %>%
  select(cluster, S, Ca) %>%
  pivot_longer(cols = c(2,3), names_to = "analyte", values_to = "mg_per_L") %>%
  drop_na() %>%
  group_by(cluster, analyte) -> alaska_lake_data_pca_clusters_clean

plot_2 <- ggplot() + 
  geom_col(
    data = summarize(
      alaska_lake_data_pca_clusters_clean,
      mean = mean(mg_per_L), sd = sd(mg_per_L)
    ),
    aes(x = cluster, y = mean, fill = cluster),
    color = "black", alpha = 0.6
  ) +
  geom_errorbar(
    data = summarize(
      alaska_lake_data_pca_clusters_clean,
      mean = mean(mg_per_L), sd = sd(mg_per_L)
    ),
    aes(x = cluster, ymin = mean-sd, ymax = mean+sd, fill = cluster),
    color = "black", alpha = 0.6, width = 0.5, size = 1
  ) +
  facet_grid(.~analyte) + theme_bw() +
  geom_jitter(
    data = alaska_lake_data_pca_clusters_clean,
    aes(x = cluster, y = mg_per_L, fill = cluster), width = 0.05,
    shape = 21
  ) +
  scale_fill_manual(values = discrete_palette) 

plot_1<- ggplot() +
  geom_point(
    data = alaska_lake_data_pca_clusters,
    aes(x = Dim.1, y = Dim.2, fill = cluster), shape = 21, size = 5, alpha = 0.6
  ) +
  geom_mark_ellipse(
    data = drop_na(alaska_lake_data_pca_clusters),
    aes(x = Dim.1, y = Dim.2, label = cluster, fill = cluster), alpha = 0.2
  ) +
  theme_classic() + coord_cartesian(xlim = c(-7,12), ylim = c(-4,5)) +
  scale_fill_manual(values = discrete_palette)

plot_1 + plot_2
```

<img src="index_files/figure-html/unnamed-chunk-136-1.png" width="100%" style="display: block; margin: auto;" />
 
## {-}

## further reading {-}

http://www.sthda.com/english/wiki/wiki.php?id_contents=7940

https://ryanwingate.com/intro-to-machine-learning/unsupervised/hierarchical-and-density-based-clustering/

https://ryanwingate.com/intro-to-machine-learning/unsupervised/hierarchical-and-density-based-clustering/hierarchical-4.png

https://www.geeksforgeeks.org/dbscan-clustering-in-r-programming/


<!-- ## exercises {-}

For this set of exercises, please use the dataset `hawaii_aquifers` or `tequila_chemistry`, available after you run the `source()` command. Do the following:

1. Run a PCA analysis on the data set and plot the results. 

2. Create an ordination plot and identify one analyte that varies with Dim.1 and one analyte that varies with Dim.2 (these are your "variables of interest").

3. Run kmeans clustering on your PCA output. Create a set of clusters that seems to appropriately subdivide the data set.

4. Use the clusters defined by kmeans as groupings on which to run summary statistics for your two variables of interest.

5. Create a plot with three subpanels that shows: (i) the PCA analysis (colored by kmeans clusters), (ii) the ordination analysis, and (iii) the summary statistics for your two variables of interest within the kmeans groups. Please note that subpanel plots can be created by sending ggplots to their own objects and then adding those objects together. Please see the subsection in the data visualization chapter on subplots.

6. Run dbscan clustering on your PCA output. Create a set of clusters that seems to appropriately subdivide the data set.

7. Use the clusters defined by dbscan as groupings on which to run summary statistics for your two variables of interest.

8. Create a plot with three subpanels that shows: (i) the PCA analysis (colored by dbscan clusters), (ii) the ordination analysis, and (iii) the summary statistics for your two variables of interest within the dbscan groups. -->


<!-- 



* Question 1

Run a principal components analysis on the dataset. Use `na_replacement = "drop"` (so that variables with NA values are not included in the analysis) and generate clusters automatically using kmeans by setting `kmeans = "auto"`. Make scatter plot of the results. How many clusters does kmeans recommend? 



* Question 2

Modify your code from Question 1 so that only two clusters are generated. Plot the results. Use `geom_mark_ellipse` to highlight each cluster in your plot (note that the `fill` aesthetic is required to mark groups). Which varieties are put into each of the two clusters?



* Question 3

Use an ordination plot to determine what chemicals makes Chardonnay so different from the other varieties. To what class of compounds do these chemical belong?



* Question 4

Modify your code from Question 2 so that five clusters are generated. Plot the results. Use `geom_mark_ellipse` to highlight each cluster in your plot (note that the `fill` aesthetic is required to mark groups). Based on this plot, which grape variety undergoes the least amount of change, chemically speaking, between dry and well-watered conditions?



* Question 5

Run a heirarchical clustering analysis on the wine grapes data set, using kmeans to create five groups, and also continue using `na_replacement = "drop"`. Plot the results. Which grape variety undergoes the most change in terms of its chemistry between well-watered and dry conditions? (hint: remember that the x-axis shows the distances between nodes and tips, the y-axis is meaningless). Compare the method you used to compare sample shifts in question 4 (i.e. pca+kmeans) versus the method you used in this question (i.e. hclust+kmeans). Which do you is better? Would this change depending on the circumstances?



* Question 6

Google "Quercetin". What kind of compound is it? Use the clusters created by the heirarchical clustering analysis in question 5 as groups for which to calculate summary statistics. Calculate the mean and standard deviation of the concentration of Quercetin in each group. Plot the result using `geom_pointrange` and adjust axis font sizes so that they are in good proportion with the size of the plot. Also specify a theme (for example, `theme_classic()`).

Does one cluster have a large amount of variation in Quercetin abundance? Why do you think this might be?



* Question 7

Use `cowplot::plot_grid` to display your plots from questions 4 and 5 next to each other.



* Challenge (optional)

Use cowplot to display your plots from questions 4, 5, and 6 alongside each other. **Make your combined plot as attractive as possible!** Use each of the following:

`align = TRUE` inside `geom_tiplab()`

`nrow = 1` inside `plot_grid()`

`rel_widths = <your_choice>` inside `plot_grid()`

`name = <your_choice>` inside `scale_*_*`

`label = cluster` inside `geom_mark_ellipse()`

`breaks = <your_choice>` inside `scale_x_continuous()` or `scale_y_continuous()` (as an example, `breaks = seq(0,10,1)`)

Also, consider using:

`guides(fill = "none", color = "none")`

Install the RColorBrewer package, and use one of its color schemes. As an example with the color scheme `Set1`:

`scale_fill_brewer(palette = "Set1", na.value = "grey")`

`scale_color_brewer(palette = "Set1", na.value = "grey")`

Save your plot as a png using `ggsave()`.



Maybe something like this:

{r fig.align='center', echo=FALSE, include=identical(knitr:::pandoc_to(), 'html'), results="markup"}
knitr:::include_graphics('https://thebustalab.github.io/integrated_bioanalytics/images/Ex6Challenge.png', dpi = NA)
 -->

<!-- end -->

<!-- start comparing means -->


# hierarchical clustering {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/clustering.png" width="100%" style="display: block; margin: auto;" />

"Which of my samples are most closely related?"

## {-}

So far we have been looking at how to plot raw data, summarize data, and reduce a data set's dimensionality. It's time to look at how to identify relationships between the samples in our data sets. For example: in the Alaska lakes dataset, which lake is most similar, chemically speaking, to Lake Narvakrak? Answering this requires calculating numeric distances between samples based on their chemical properties. For this, the first thing we need is a distance matrix:

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/dist_matrix.jpg" width="100%" style="display: block; margin: auto;" />

Please note that we can get distance matrices directly from `runMatrixAnalysis` by specifying `analysis = "dist"`:


``` r
dist <- runMatrixAnalysis(
    data = alaska_lake_data,
    analysis = c("dist"),
    column_w_names_of_multiple_analytes = "element",
    column_w_values_for_multiple_analytes = "mg_per_L",
    columns_w_values_for_single_analyte = c("water_temp", "pH"),
    columns_w_additional_analyte_info = "element_type",
    columns_w_sample_ID_info = c("lake", "park")
)
## Replacing NAs in your data with mean
```

``` r

as.matrix(dist)[1:3,1:3]
##                          Devil_Mountain_Lake_BELA
## Devil_Mountain_Lake_BELA                 0.000000
## Imuruk_Lake_BELA                         3.672034
## Kuzitrin_Lake_BELA                       1.663147
##                          Imuruk_Lake_BELA
## Devil_Mountain_Lake_BELA         3.672034
## Imuruk_Lake_BELA                 0.000000
## Kuzitrin_Lake_BELA               3.062381
##                          Kuzitrin_Lake_BELA
## Devil_Mountain_Lake_BELA           1.663147
## Imuruk_Lake_BELA                   3.062381
## Kuzitrin_Lake_BELA                 0.000000
```

There is more that we can do with distance matrices though, lots more. Let's start by looking at an example of hierarchical clustering. For this, we just need to tell `runMatrixAnalysis()` to use `analysis = "hclust"`: 


``` r
AK_lakes_clustered <- runMatrixAnalysis(
    data = alaska_lake_data,
    analysis = "hclust",
    column_w_names_of_multiple_analytes = "element",
    column_w_values_for_multiple_analytes = "mg_per_L",
    columns_w_values_for_single_analyte = c("water_temp", "pH"),
    columns_w_additional_analyte_info = "element_type",
    columns_w_sample_ID_info = c("lake", "park"),
    na_replacement = "mean"
)
AK_lakes_clustered
## # A tibble: 39 × 26
##    sample_unique_ID   lake  park  parent  node branch.length
##    <chr>              <chr> <chr>  <int> <int>         <dbl>
##  1 Devil_Mountain_La… Devi… BELA      33     1          8.12
##  2 Imuruk_Lake_BELA   Imur… BELA      32     2          4.81
##  3 Kuzitrin_Lake_BELA Kuzi… BELA      37     3          3.01
##  4 Lava_Lake_BELA     Lava… BELA      38     4          2.97
##  5 North_Killeak_Lak… Nort… BELA      21     5        254.  
##  6 White_Fish_Lake_B… Whit… BELA      22     6         80.9 
##  7 Iniakuk_Lake_GAAR  Inia… GAAR      29     7          3.60
##  8 Kurupa_Lake_GAAR   Kuru… GAAR      31     8          8.57
##  9 Lake_Matcharak_GA… Lake… GAAR      29     9          3.60
## 10 Lake_Selby_GAAR    Lake… GAAR      30    10          4.80
## # ℹ 29 more rows
## # ℹ 20 more variables: label <chr>, isTip <lgl>, x <dbl>,
## #   y <dbl>, branch <dbl>, angle <dbl>, bootstrap <dbl>,
## #   water_temp <dbl>, pH <dbl>, C <dbl>, N <dbl>, P <dbl>,
## #   Cl <dbl>, S <dbl>, F <dbl>, Br <dbl>, Na <dbl>,
## #   K <dbl>, Ca <dbl>, Mg <dbl>
```

It works! Now we can plot our cluster diagram with a ggplot add-on called ggtree. We've seen that ggplot takes a "data" argument (i.e. `ggplot(data = <some_data>) + geom_*()` etc.). In contrast, ggtree takes an argument called `tr`, though if you're using the `runMatrixAnalysis()` function, you can treat these two (`data` and `tr`) the same, so, use: `ggtree(tr = <output_from_runMatrixAnalysis>) + geom_*()` etc.

Note that `ggtree` also comes with several great new geoms: `geom_tiplab()` and `geom_tippoint()`. Let's try those out:


``` r
library(ggtree)
AK_lakes_clustered %>%
ggtree() +
  geom_tiplab() +
  geom_tippoint() +
  theme_classic()
```

<img src="index_files/figure-html/unnamed-chunk-149-1.png" width="100%" style="display: block; margin: auto;" />

Cool! Though that plot could use some tweaking... let's try:


``` r
AK_lakes_clustered %>%
ggtree() +
    geom_tiplab(aes(label = lake), offset = 10, align = TRUE) +
    geom_tippoint(shape = 21, aes(fill = park), size = 4) +
    scale_x_continuous(limits = c(0,375)) +
    scale_fill_brewer(palette = "Set1") +
    # theme_classic() +
    theme(
      legend.position = c(0.2,0.8)
    )
```

<img src="index_files/figure-html/unnamed-chunk-150-1.png" width="100%" style="display: block; margin: auto;" />

Very nice!

## further reading {-}
 
For more information on plotting annotated trees, see: https://yulab-smu.top/treedata-book/chapter10.html.

For more on clustering, see: https://ryanwingate.com/intro-to-machine-learning/unsupervised/hierarchical-and-density-based-clustering/.

## exercises {-}

For this set of exercises, please use `runMatrixAnalysis()` to run and visualize a hierarchical cluster analysis with each of the main datasets that we have worked with so far, except for NY_trees. This means: `algae_data` (which algae strains are most similar to each other?), `alaska_lake_data` (which lakes are most similar to each other?). and `solvents` (which solvents are most similar to each other?). It also means you should use the periodic table (which elements are most similar to each other?), though please don't use the whole periodic table, rather, use `periodic_table_subset`. Please also conduct a heirarchical clustering analysis for a dataset of your own choice that is not provided by the `source()` code. For each of these, create (i) a tree diagram that shows how the "samples" in each data set are related to each other based on the numerical data associated with them, (ii) a caption for each diagram, and (iii) describe, in two or so sentences, an interesting trend you see in the diagram. You can ignore columns that contain categorical data, or you can list those columns as "additional_analyte_info".

For this assignment, you may again find the `colnames()` function and square bracket-subsetting useful. It will list all or a subset of the column names in a dataset for you. For example:


``` r
colnames(solvents)
##  [1] "solvent"             "formula"            
##  [3] "boiling_point"       "melting_point"      
##  [5] "density"             "miscible_with_water"
##  [7] "solubility_in_water" "relative_polarity"  
##  [9] "vapor_pressure"      "CAS_number"         
## [11] "formula_weight"      "refractive_index"   
## [13] "specific_gravity"    "category"
```

``` r

colnames(solvents)[1:3]
## [1] "solvent"       "formula"       "boiling_point"
```

``` r

colnames(solvents)[c(1,5,7)]
## [1] "solvent"             "density"            
## [3] "solubility_in_water"
```

# comparing means {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/hawaii_aquifers.jpeg" width="100%" style="display: block; margin: auto;" />

**"Are these two things the same?"**

## {-}

Often, we want to know if our study subjects contain different amounts of certain analytes. For example, "Does this lake over here contain more potassium than that lake over there?" For this, we need statistical tests. Here, we will have a look at comparing mean values for analyte abundance in situations with two samples and in situations with more than two samples.

I find many of the concepts discussed in this chapter easier to think about with an example in mind. For that, suppose that you are an analytical chemist on Hawaii that is studying the chemistry of the island's aquifers. you have the data set `hawaii_aquifers`. You can see in the output below the structure of the data set - we have 990 measurements of a 9 different analytes in multiple wells that draw on a set of 10 aquifers.


``` r
hawaii_aquifers
## # A tibble: 954 × 6
##    aquifer_code well_name         longitude latitude analyte
##    <chr>        <chr>                 <dbl>    <dbl> <chr>  
##  1 aquifer_1    Alewa_Heights_Sp…        NA       NA SiO2   
##  2 aquifer_1    Alewa_Heights_Sp…        NA       NA Cl     
##  3 aquifer_1    Alewa_Heights_Sp…        NA       NA Mg     
##  4 aquifer_1    Alewa_Heights_Sp…        NA       NA Na     
##  5 aquifer_1    Alewa_Heights_Sp…        NA       NA K      
##  6 aquifer_1    Alewa_Heights_Sp…        NA       NA SO4    
##  7 aquifer_1    Alewa_Heights_Sp…        NA       NA HCO3   
##  8 aquifer_1    Alewa_Heights_Sp…        NA       NA dissol…
##  9 aquifer_1    Alewa_Heights_Sp…        NA       NA Ca     
## 10 aquifer_1    Beretania_High_S…        NA       NA SiO2   
## # ℹ 944 more rows
## # ℹ 1 more variable: abundance <dbl>
```

``` r
unique(hawaii_aquifers$aquifer_code)
##  [1] "aquifer_1"  "aquifer_2"  "aquifer_3"  "aquifer_4" 
##  [5] "aquifer_5"  "aquifer_6"  "aquifer_7"  "aquifer_8" 
##  [9] "aquifer_9"  "aquifer_10"
```
Importantly, there are many wells that draw on each aquifer, as shown in the graph below.


``` r
hawaii_aquifers %>%
  select(aquifer_code, well_name) %>%
  group_by(aquifer_code) %>%
  summarize(n_wells = length(unique(well_name))) -> aquifers_summarized

aquifers_summarized
## # A tibble: 10 × 2
##    aquifer_code n_wells
##    <chr>          <int>
##  1 aquifer_1         12
##  2 aquifer_10         7
##  3 aquifer_2          5
##  4 aquifer_3          3
##  5 aquifer_4         16
##  6 aquifer_5          4
##  7 aquifer_6         12
##  8 aquifer_7          9
##  9 aquifer_8          3
## 10 aquifer_9         30
```

``` r

ggplot(aquifers_summarized) + geom_col(aes(x = n_wells, y = aquifer_code))
```

<img src="index_files/figure-html/unnamed-chunk-154-1.png" width="100%" style="display: block; margin: auto;" />

<!-- To run these statistical analyses, we will need several new R packages: `rstatix`, `agricolae`, and `multcompView`. Please install these with `install.packages("rstatix")`, `install.packages("agricolae")`, and `install.packages("multcompView")`. Load them into your R session using `library(rstatix)`, `library(agricolae)`, and `library(multcompView)`.
 -->
<!-- ```{r, echo = FALSE, message = FALSE}
source("https://thebustalab.github.io/phylochemistry/phylochemistry.R")
library(rstatix)
library(agricolae)
library(multcompView)
``` -->

## definitions {-}

1. **populations and independent measurements**: When we are comparing means, we are comparing two *sets* of values. It is important to consider where these values came from in the first place. In particular, it is usually useful to think of these values as representatives of larger populations. In the example of our aquifer data set, we can think of the measurements from different wells drawing on the same aquifer as independent measurements of the "population" (i.e. the aquifer).

2. **the null hypothesis**: When we conduct a statistical test, we are testing the null hypothesis. The null (think "default") hypothesis is that there is no difference bewteen the means (hence the name "null"). In the example of our aquifers, let's say that we're interested in whether two aquifers have different abundances of potassium - in this case the null hypothesis is that they do not differ, in other words, that they have the same amount of potassium.

3. **the *p* value**: The *p* value represents the probability of getting data as extreme as our results if the null hypothesis is true. In other words - the *p* value is the probability that we would observe the differences we did, if in fact there were no differences in the means at all. To continue with our example: suppose we measure potassium levels in 10% of the wells that access each aquifer and find that aquifer_1 has potassium levels of 14 +/- 2 and aquifer_2 has potassium levels of 12 +/- 1. Suppose that we then conduct a statistical test and get a p value of 0.04. This means that, assuming the aquifers have the same magneisum levels (i.e. assuming the null hypothesis is true), there is a 4% chance that we would get the measured values that we did. In other words, IF the aquifers have the same potassium abundance, it is pretty unlikely that we would have obtained the measurements that we did.

Please note that the the *p* value is not the probability of a detected difference being a false positive. The probability of a false positive requires additional information in order to be calculated. For further discussion please see the end of this chapter.

## test selection {-}

There are many different types of statistical tests. Below is a flow chart illustrating how it is recommended that statistical tests be used in this course. You can see that there are three regimes of tests: variance and normality tests (blue), parametric tests (green), and non-parametric tests (orange):

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/stats.png" width="100%" style="display: block; margin: auto;" />

When we are comparing means, we need to first determine what kind of statistical tests we can use with our data. If (i) our data can be reasonably modelled by a normal distribution and (ii) the variances about the two means are similar, then we can use the more powerful "parametric" tests (i.e. tests that will be more likely to detect a difference in means, assuming one exists). If one of these criteria are not met, then we need to use less powerful "non-parametric" tests.

We can check our data for normality and similar variances using the Shapiro test and the Levene test. Let's use the `hawaii_aquifers` data as an example, and let's consider only the element potassium:


``` r
K_data <- hawaii_aquifers %>%
  filter(analyte == "K")
  K_data
## # A tibble: 106 × 6
##    aquifer_code well_name         longitude latitude analyte
##    <chr>        <chr>                 <dbl>    <dbl> <chr>  
##  1 aquifer_1    Alewa_Heights_Sp…       NA      NA   K      
##  2 aquifer_1    Beretania_High_S…       NA      NA   K      
##  3 aquifer_1    Beretania_Low_Se…       NA      NA   K      
##  4 aquifer_1    Kuliouou_Well         -158.     21.3 K      
##  5 aquifer_1    Manoa_Well_II         -158.     21.3 K      
##  6 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  7 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  8 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  9 aquifer_1    Nuuanu_Aerator_W…     -158.     21.4 K      
## 10 aquifer_1    Palolo_Tunnel         -158.     21.3 K      
## # ℹ 96 more rows
## # ℹ 1 more variable: abundance <dbl>
```

To work with two means, let's just look at aquifers 1 and 6:


``` r
K_data_1_6 <- K_data %>%
    filter(aquifer_code %in% c("aquifer_1", "aquifer_6"))

K_data_1_6
## # A tibble: 24 × 6
##    aquifer_code well_name         longitude latitude analyte
##    <chr>        <chr>                 <dbl>    <dbl> <chr>  
##  1 aquifer_1    Alewa_Heights_Sp…       NA      NA   K      
##  2 aquifer_1    Beretania_High_S…       NA      NA   K      
##  3 aquifer_1    Beretania_Low_Se…       NA      NA   K      
##  4 aquifer_1    Kuliouou_Well         -158.     21.3 K      
##  5 aquifer_1    Manoa_Well_II         -158.     21.3 K      
##  6 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  7 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  8 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  9 aquifer_1    Nuuanu_Aerator_W…     -158.     21.4 K      
## 10 aquifer_1    Palolo_Tunnel         -158.     21.3 K      
## # ℹ 14 more rows
## # ℹ 1 more variable: abundance <dbl>
```

``` r

ggplot(K_data_1_6, aes(x = aquifer_code, y = abundance)) +
    geom_boxplot() +
    geom_point()
```

<img src="index_files/figure-html/unnamed-chunk-157-1.png" width="100%" style="display: block; margin: auto;" />

Are these data normally distributed? Do they have similar variance? Let's get a first approximation by looking at a plot:


``` r
K_data_1_6 %>%
  ggplot(aes(x = abundance)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~aquifer_code) +
    geom_density(aes(y = ..density..*10), color = "blue")
```

<img src="index_files/figure-html/unnamed-chunk-158-1.png" width="100%" style="display: block; margin: auto;" />

Based on this graphic, it's hard to say! Let's use a statistical test to help. When we want to run the Shaprio test, we are looking to see if each group has normally distributed here (here group is "aquifer_code", i.e. aquifer_1 and aquifer_6). This means we need to `group_by(aquifer_code)` before we run the test:


``` r
K_data_1_6 %>%
  group_by(aquifer_code) %>% 
  shapiroTest(abundance)
## # A tibble: 2 × 4
##   aquifer_code variable  statistic     p
##   <chr>        <chr>         <dbl> <dbl>
## 1 aquifer_1    abundance     0.885 0.102
## 2 aquifer_6    abundance     0.914 0.239
```

Both p-values are above 0.05! This means that the distributions are not significantly different from a normal distribution. What about the variances about the two means? Are they similar? For this we need a Levene test. With that test, we are not looking within each group, but rather across groups - this means we do NOT need to `group_by(aquifer_code)` and should specify a `y ~ x` formula instead:


``` r
K_data_1_6 %>%
  leveneTest(abundance ~ aquifer_code)
## # A tibble: 1 × 4
##     df1   df2 statistic     p
##   <int> <int>     <dbl> <dbl>
## 1     1    22     0.289 0.596
```

The p-value from this test is 0.596! This means that their variances are not significantly different. If their shape is the same (normality) and their variances are the same (equal variances), then the only thing that can be different is their means. This allows us to set up our null hypothesis and calculate the probability of obtaining these different means if the two sets of observations come from an identical source.

## two means {-}

Now, since our data passed both test, this means we can use a normal t-test. A t-test is a parametric test. This means that it relies on modelling the data using a normal distribution in order to make comparisons. It is also a powerful test. This means that it is likely to detect a difference in means, assuming one is present. Let's try it out:


``` r
K_data_1_6 %>%
  tTest(abundance ~ aquifer_code)
## # A tibble: 1 × 8
##   .y.       group1 group2    n1    n2 statistic    df      p
## * <chr>     <chr>  <chr>  <int> <int>     <dbl> <dbl>  <dbl>
## 1 abundance aquif… aquif…    12    12     -2.75  20.5 0.0121
```

A p-value of 0.012! This is below 0.05, meaning that there is a 95% chance that the two means are different. Suppose that our data had not passed the Shapiro and/or Levene tests. We would then need to use a Wilcox test. The Wilcox test is a non-parametric test, which means that it does not use a normal distribution to model the data in order to make comparisons. This means that is a less powerful test than the t-test, which means that it is less likely to detect a difference in the means, assuming there is one. For fun, let's try that one out and compare the p-values from the two methods:


``` r
K_data_1_6 %>%
  wilcoxTest(abundance ~ aquifer_code)
## # A tibble: 1 × 7
##   .y.       group1    group2       n1    n2 statistic      p
## * <chr>     <chr>     <chr>     <int> <int>     <dbl>  <dbl>
## 1 abundance aquifer_1 aquifer_6    12    12      33.5 0.0282
```

A p-value of 0.028! This is higher than the value given by the t-test (0.012). That is because the Wilcox test is a less powerful test: it is less likely to detect differences in means, assuming they exist.

## more than two means {-}

In the previous section we compared two means. What if we want to compare means from more than two study subjects? The first step is again to determine which tests to use. Let's consider our hawaii aquifer data again, though this time let's use all the aquifers, not just two:


``` r
K_data <- hawaii_aquifers %>%
  filter(analyte == "K")

K_data
## # A tibble: 106 × 6
##    aquifer_code well_name         longitude latitude analyte
##    <chr>        <chr>                 <dbl>    <dbl> <chr>  
##  1 aquifer_1    Alewa_Heights_Sp…       NA      NA   K      
##  2 aquifer_1    Beretania_High_S…       NA      NA   K      
##  3 aquifer_1    Beretania_Low_Se…       NA      NA   K      
##  4 aquifer_1    Kuliouou_Well         -158.     21.3 K      
##  5 aquifer_1    Manoa_Well_II         -158.     21.3 K      
##  6 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  7 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  8 aquifer_1    Moanalua_Wells_P…     -158.     21.4 K      
##  9 aquifer_1    Nuuanu_Aerator_W…     -158.     21.4 K      
## 10 aquifer_1    Palolo_Tunnel         -158.     21.3 K      
## # ℹ 96 more rows
## # ℹ 1 more variable: abundance <dbl>
```

``` r

ggplot(data = K_data, aes(y = aquifer_code, x = abundance)) +
  geom_boxplot() +
  geom_point(color = "maroon", alpha = 0.6, size = 3)
```

<img src="index_files/figure-html/unnamed-chunk-163-1.png" width="100%" style="display: block; margin: auto;" />

Let's check visually to see if each group is normally distributed and to see if they have roughly equal variance:


``` r
K_data %>%
  group_by(aquifer_code) %>%
  ggplot(aes(x = abundance)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~aquifer_code) +
    geom_density(aes(y = ..density..*10), colour = "blue")
```

<img src="index_files/figure-html/unnamed-chunk-164-1.png" width="100%" style="display: block; margin: auto;" />

Again, it is somewhat hard to tell visually if these data are normally distributed. It seems pretty likely that they have different variances about the means, but let's check using the Shapiro and Levene tests. Don't forget: with the Shaprio test, we are looking within each group and so need to `group_by()`, with the Levene test, we are looking across groups, and so need to provide a `y~x` formula:


``` r
K_data %>%
  group_by(aquifer_code) %>% 
  shapiroTest(abundance)
## # A tibble: 10 × 4
##    aquifer_code variable  statistic         p
##    <chr>        <chr>         <dbl>     <dbl>
##  1 aquifer_1    abundance     0.885 0.102    
##  2 aquifer_10   abundance     0.864 0.163    
##  3 aquifer_2    abundance     0.913 0.459    
##  4 aquifer_3    abundance     0.893 0.363    
##  5 aquifer_4    abundance     0.948 0.421    
##  6 aquifer_5    abundance     0.993 0.972    
##  7 aquifer_6    abundance     0.914 0.239    
##  8 aquifer_7    abundance     0.915 0.355    
##  9 aquifer_8    abundance     0.842 0.220    
## 10 aquifer_9    abundance     0.790 0.0000214
```


``` r
K_data %>%
  leveneTest(abundance ~ aquifer_code)
## # A tibble: 1 × 4
##     df1   df2 statistic       p
##   <int> <int>     <dbl>   <dbl>
## 1     9    96      2.95 0.00387
```

Based on these tests, it looks like the data for aquifer 9 is significantly different from a normal distribution (Shaprio test p < 0.05), and the variances are certainly different from one another (Levene test p > 0.05).

Let's assume for a second that our data passed these tests. This means that we could reasonably model our data with normal distributions and use a parametric test to compare means. This means that we can use an ANOVA to test for differences in means.

### ANOVA, Tukey tests {-}

We will use the `anovaTest` function from the package `rstatix`. It will tell us if any of the means in the data are statistically different from one another. However, if there are differences between the means, it will not tell us which of them are different.


``` r
K_data %>%
  anovaTest(abundance ~ aquifer_code)
## ANOVA Table (type II tests)
## 
##         Effect DFn DFd     F        p p<.05   ges
## 1 aquifer_code   9  96 9.486 3.28e-10     * 0.471
```

A pretty small p-value! There are definitely some significant differences among this group. But, WHICH are different from one another though? For this, we need to run Tukey's Honest Significant Difference test (implemented using `tukey_hsd`). This will essentially run t-test on all the pairs of study subjects that we can derive from our data set (in this example, aquifer_1 vs. aquifer_2, aquifer_1 vs. aquifer_3, etc.). After that, it will correct the p-values according to the number of comparisons that it performed. This controls the rate of type I error that we can expect from the test. These corrected values are provided to us in the `p.adj` column.


``` r
K_data %>%
  tukey_hsd(abundance ~ aquifer_code)
## # A tibble: 45 × 9
##    term         group1   group2 null.value estimate conf.low
##  * <chr>        <chr>    <chr>       <dbl>    <dbl>    <dbl>
##  1 aquifer_code aquifer… aquif…          0  0.00357   -2.04 
##  2 aquifer_code aquifer… aquif…          0  1.44      -0.708
##  3 aquifer_code aquifer… aquif…          0  0.375     -2.40 
##  4 aquifer_code aquifer… aquif…          0 -1.15      -2.78 
##  5 aquifer_code aquifer… aquif…          0 -0.875     -3.36 
##  6 aquifer_code aquifer… aquif…          0  1.98       0.228
##  7 aquifer_code aquifer… aquif…          0  2.70       0.801
##  8 aquifer_code aquifer… aquif…          0 -0.125     -2.90 
##  9 aquifer_code aquifer… aquif…          0 -0.349     -1.80 
## 10 aquifer_code aquifer… aquif…          0  1.44      -0.954
## # ℹ 35 more rows
## # ℹ 3 more variables: conf.high <dbl>, p.adj <dbl>,
## #   p.adj.signif <chr>
```

``` r

# K_data %>%
#   tukeyHSD(abundance ~ aquifer_code)
```

Using the output from our tukey test, we can determine which means are similar. We can do this using the pGroups function:


``` r
groups_based_on_tukey <- K_data %>%
  tukey_hsd(abundance ~ aquifer_code) %>%
  pGroups()
groups_based_on_tukey
##             treatment group spaced_group
## aquifer_1   aquifer_1    ab         ab  
## aquifer_10 aquifer_10   abc         abc 
## aquifer_2   aquifer_2   acd         a cd
## aquifer_3   aquifer_3  abcd         abcd
## aquifer_4   aquifer_4     b          b  
## aquifer_5   aquifer_5    ab         ab  
## aquifer_6   aquifer_6    cd           cd
## aquifer_7   aquifer_7     d            d
## aquifer_8   aquifer_8  abcd         abcd
## aquifer_9   aquifer_9    ab         ab
```

We can use the output from `pGroups` to annotate our plot:


``` r
ggplot(data = K_data, aes(y = aquifer_code, x = abundance)) +
  geom_boxplot() +
  geom_point(color = "maroon", alpha = 0.6, size = 3) +
  geom_text(data = groups_based_on_tukey, aes(y = treatment, x = 9, label = group))
```

<img src="index_files/figure-html/unnamed-chunk-170-1.png" width="100%" style="display: block; margin: auto;" />

Excellent! This plot shows us, using the letters on the same line with each aquifer, which means are the same and which are different. If a letter is shared among the labels in line with two aquifers, it means that their means do not differ significantly. For example, aquifer 2 and aquifer 6 both have "b" in their labels, so their means are not different - and are the same as those of aquifers 3 and 10.

### Kruskal, Dunn tests {-}

The above ANOVA example is great, but remember - our data did not pass the Shapiro or Levene tests. This means not all our data can be modelled by a normal distribution and taht we need to use a non-parametric test. The non-parametric alternative to the ANOVA is called the Kruskal test. Like the Wilcox test, it is less powerful that its parametric relative, meaning that it is less likely to detected differences, should they exist. However, since our data do not pass the Shapiro/Levene tests, we have to resort to the Kruskal test. Let's try it out:


``` r
K_data %>%
  kruskalTest(abundance ~ aquifer_code)
## # A tibble: 1 × 6
##   .y.           n statistic    df             p method      
## * <chr>     <int>     <dbl> <int>         <dbl> <chr>       
## 1 abundance   106      55.9     9 0.00000000807 Kruskal-Wal…
```

A pretty small p-value! This is higher than the p-value from running ANOVA on the same data (remember, the Kruskal test is less powerful). Never the less, the value is still well below 0.05, meaning that some of the means are different. So, how do we determine WHICH are different from one another? When we ran ANOVA the follow-up test (the post hoc test) was Tukey's HSD. After the Kruskal test, the post hoc test we use is the Dunn test. Let's try:


``` r
K_data %>%
  dunnTest(abundance ~ aquifer_code)
## # A tibble: 45 × 9
##    .y.     group1 group2    n1    n2 statistic       p p.adj
##  * <chr>   <chr>  <chr>  <int> <int>     <dbl>   <dbl> <dbl>
##  1 abunda… aquif… aquif…    12     7    -0.194 0.846   1    
##  2 abunda… aquif… aquif…    12     6     2.24  0.0254  0.736
##  3 abunda… aquif… aquif…    12     3     0.866 0.387   1    
##  4 abunda… aquif… aquif…    12    17    -2.65  0.00806 0.266
##  5 abunda… aquif… aquif…    12     4    -1.12  0.263   1    
##  6 abunda… aquif… aquif…    12    12     2.51  0.0121  0.388
##  7 abunda… aquif… aquif…    12     9     3.01  0.00257 0.100
##  8 abunda… aquif… aquif…    12     3     0.143 0.886   1    
##  9 abunda… aquif… aquif…    12    33    -0.470 0.639   1    
## 10 abunda… aquif… aquif…     7     6     2.17  0.0296  0.830
## # ℹ 35 more rows
## # ℹ 1 more variable: p.adj.signif <chr>
```

This gives us adjusted p-values for all pairwise comparisons. Once again, we can use `pGroups()` to give us a compact letter display for each group, which can then be used to annotate the plot:


``` r
groups_based_on_dunn <- K_data %>%
  dunnTest(abundance ~ aquifer_code) %>%
  pGroups()
groups_based_on_dunn
##             treatment group spaced_group
## aquifer_1   aquifer_1  abcd         abcd
## aquifer_10 aquifer_10  abcd         abcd
## aquifer_2   aquifer_2   abc         abc 
## aquifer_3   aquifer_3  abcd         abcd
## aquifer_4   aquifer_4     d            d
## aquifer_5   aquifer_5   acd         a cd
## aquifer_6   aquifer_6    ab         ab  
## aquifer_7   aquifer_7     b          b  
## aquifer_8   aquifer_8  abcd         abcd
## aquifer_9   aquifer_9    cd           cd
```

``` r

ggplot(data = K_data, aes(y = aquifer_code, x = abundance)) +
  geom_boxplot() +
  geom_point(color = "black", alpha = 0.4, size = 2) +
  scale_x_continuous(name = "Potassium abundance", breaks = seq(0,10,1)) +
  scale_y_discrete(name = "Aquifer code") +
  geom_text(data = groups_based_on_dunn, aes(y = treatment, x = 9, label = group)) +
  theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-173-1.png" width="100%" style="display: block; margin: auto;" />

Note that these groupings are different from those generated by ANOVA/Tukey.

## pairs of means {-}

Oftentimes we have more than two means to compare, but rather than wanting to compare all means at once, we want to compare them in a pairwise fashion. For example, suppose we want to know if any of the aquifers contain different amounts of Na and Cl. We are not interested in testing for differences among *all* values of Na and Cl, rather, we want to test all *pairs* of Na and Cl values arising from each aquifer. That is to say, we want to compare the means in each facet of the plot below:


``` r
hawaii_aquifers %>%
  filter(analyte %in% c("Na", "Cl")) %>%
  ggplot(aes(x = analyte, y = abundance)) + geom_violin() + geom_point() + facet_grid(.~aquifer_code)
```

<img src="index_files/figure-html/unnamed-chunk-174-1.png" width="100%" style="display: block; margin: auto;" />

Fortunately, we can use an approach that is very similar to the what we've learned in the earlier portions of this chapter, just with minor modifications. Let's have a look! We start with the Shapiro and Levene tests, as usual (note that we group using two variables when using the Shapiro test so that each analyte within each aquifer is considered as an individual distribution):


``` r
hawaii_aquifers %>%
  filter(analyte %in% c("Na", "Cl")) %>%
  group_by(analyte, aquifer_code) %>%
  shapiroTest(abundance)
## # A tibble: 20 × 5
##    aquifer_code analyte variable  statistic        p
##    <chr>        <chr>   <chr>         <dbl>    <dbl>
##  1 aquifer_1    Cl      abundance     0.900 1.59e- 1
##  2 aquifer_10   Cl      abundance     0.486 1.09e- 5
##  3 aquifer_2    Cl      abundance     0.869 2.24e- 1
##  4 aquifer_3    Cl      abundance     0.75  0       
##  5 aquifer_4    Cl      abundance     0.903 7.49e- 2
##  6 aquifer_5    Cl      abundance     0.849 2.24e- 1
##  7 aquifer_6    Cl      abundance     0.741 2.15e- 3
##  8 aquifer_7    Cl      abundance     0.893 2.12e- 1
##  9 aquifer_8    Cl      abundance     0.878 3.17e- 1
## 10 aquifer_9    Cl      abundance     0.420 2.68e-10
## 11 aquifer_1    Na      abundance     0.886 1.06e- 1
## 12 aquifer_10   Na      abundance     0.593 2.26e- 4
## 13 aquifer_2    Na      abundance     0.884 2.88e- 1
## 14 aquifer_3    Na      abundance     0.822 1.69e- 1
## 15 aquifer_4    Na      abundance     0.933 2.41e- 1
## 16 aquifer_5    Na      abundance     0.827 1.61e- 1
## 17 aquifer_6    Na      abundance     0.764 3.80e- 3
## 18 aquifer_7    Na      abundance     0.915 3.51e- 1
## 19 aquifer_8    Na      abundance     0.855 2.53e- 1
## 20 aquifer_9    Na      abundance     0.531 3.97e- 9
```

Looks like some of those distributions are significantly different from normal! Let's run the levene test anyway. Note that for this particular case of the Levene test, we are interested in testing whether each pair of distributions has similar variances. For that we need to feed the Levene test data that is grouped by aquifer_code (so that it tests each pair as a group), then we need to specify the y ~ x formula (which in this case is abundance ~ analyte):


``` r
hawaii_aquifers %>%
  filter(analyte %in% c("Na", "Cl")) %>%
  group_by(aquifer_code) %>%
  leveneTest(abundance ~ analyte)
## # A tibble: 10 × 5
##    aquifer_code   df1   df2 statistic       p
##    <chr>        <int> <int>     <dbl>   <dbl>
##  1 aquifer_1        1    22   10.5    0.00375
##  2 aquifer_10       1    12    0.0535 0.821  
##  3 aquifer_2        1    10    0.0243 0.879  
##  4 aquifer_3        1     4    0.320  0.602  
##  5 aquifer_4        1    32    1.57   0.219  
##  6 aquifer_5        1     6    2      0.207  
##  7 aquifer_6        1    22    1.03   0.322  
##  8 aquifer_7        1    16    1.54   0.232  
##  9 aquifer_8        1     4    0.515  0.512  
## 10 aquifer_9        1    64    1.10   0.298
```

It looks like the variances of the pair in aquifer 1 have significantly different variances. So - we for sure need to be using non-parametric testing. If this were a simple case of two means we would use the `wilcox_test`, but we have may pairs, so we will use `pairwise_wilcox_test` (note that with this test there are options for various styles of controlling for multiple comparisons, see: `?pairwise_wilcox_test`):


``` r
hawaii_aquifers %>%
  filter(analyte %in% c("Na", "Cl")) %>%
  group_by(aquifer_code) %>%
  pairwiseWilcoxTest(abundance~analyte)
## # A tibble: 10 × 10
##    aquifer_code .y.      group1 group2    n1    n2 statistic
##  * <chr>        <chr>    <chr>  <chr>  <int> <int>     <dbl>
##  1 aquifer_1    abundan… Cl     Na        12    12      99.5
##  2 aquifer_10   abundan… Cl     Na         7     7      14  
##  3 aquifer_2    abundan… Cl     Na         6     6      36  
##  4 aquifer_3    abundan… Cl     Na         3     3       3  
##  5 aquifer_4    abundan… Cl     Na        17    17     189  
##  6 aquifer_5    abundan… Cl     Na         4     4      13  
##  7 aquifer_6    abundan… Cl     Na        12    12      53  
##  8 aquifer_7    abundan… Cl     Na         9     9      42  
##  9 aquifer_8    abundan… Cl     Na         3     3       6  
## 10 aquifer_9    abundan… Cl     Na        33    33     195  
## # ℹ 3 more variables: p <dbl>, p.adj <dbl>,
## #   p.adj.signif <chr>
```

Excellent! It looks like there is a statistically significant difference between the means of the abundances of Cl and Na in aquifer_2 and (surprisingly?) in aquifer_9 (perhaps due to the large number of observations?).

What would we have done if our Shaprio and Levene tests had revealed no significant differences? Well, a `pairwise_TTest` of course!


``` r
hawaii_aquifers %>%
  filter(analyte %in% c("Na", "Cl")) %>%
  group_by(aquifer_code) %>%
  pairwiseTTest(abundance~analyte) -> test_output
  test_output
## # A tibble: 10 × 10
##    aquifer_code .y.       group1 group2    n1    n2        p
##  * <chr>        <chr>     <chr>  <chr>  <int> <int>    <dbl>
##  1 aquifer_1    abundance Cl     Na        12    12  4.69e-2
##  2 aquifer_10   abundance Cl     Na         7     7  8.82e-1
##  3 aquifer_2    abundance Cl     Na         6     6  3.75e-5
##  4 aquifer_3    abundance Cl     Na         3     3  6.83e-1
##  5 aquifer_4    abundance Cl     Na        17    17  1.03e-1
##  6 aquifer_5    abundance Cl     Na         4     4  9.75e-2
##  7 aquifer_6    abundance Cl     Na        12    12  5.66e-1
##  8 aquifer_7    abundance Cl     Na         9     9  5.21e-1
##  9 aquifer_8    abundance Cl     Na         3     3  4.28e-1
## 10 aquifer_9    abundance Cl     Na        33    33  8.96e-1
## # ℹ 3 more variables: p.signif <chr>, p.adj <dbl>,
## #   p.adj.signif <chr>
```

Excellent, now we see how to run parametric and non-parametric pairwise comparisons. How do we annotate plots with the output of these tests? Here is an example:


``` r

anno <- data.frame(
  xmin = test_output$group1,
  xmax = test_output$group2,
  y_position = c(150, 150, 150, 175, 80, 50, 300, 150, 50, 125),
  text = test_output$p.signif,
  text_size = 10,
  text_vert_offset = 10,
  text_horiz_offset = 1.5,
  tip_length_xmin = 5,
  tip_length_xmax = 5,
  aquifer_code = test_output$aquifer_code,
  hjust = 0.5,
  vjust = 0.5
)

hawaii_aquifers %>%
  filter(analyte %in% c("Na", "Cl")) %>%
  ggplot(aes(x = analyte, y = abundance)) +
  geom_violin(fill = "gold", color = "black") +
  geom_point(shape = 21, fill = "maroon", color = "black") +
  facet_grid(.~aquifer_code) +
  geomSignif(data = anno, orientation = "horizontal") +
  scale_x_discrete(name = "Analyte") +
  scale_y_continuous(name = "Abundance") +
  theme_bw() +
  theme(
    text = element_text(size = 16)
    )
```

<img src="index_files/figure-html/unnamed-chunk-179-1.png" width="100%" style="display: block; margin: auto;" />

## {-}

## further reading {-}

For more on comparing multiple means in R: [www.datanovia.com](https://www.datanovia.com/en/courses/comparing-multiple-means-in-r/)

For more on parametric versus non-parametric tests: [Statistics by Jim](https://statisticsbyjim.com/hypothesis-testing/nonparametric-parametric-tests/)

For more on interpreting *p* values: [[The p value wars (again) by Ulrich Dirnagl](https://link.springer.com/article/10.1007/s00259-019-04467-5)]

Ten common statistical mistakes and their solutions: [Science Forum: Ten common statistical mistakes to watch out for when writing or reviewing a manuscript](https://elifesciences.org/articles/48175)

How to think about very small p-values: [Reporting p Values, by Wolfgang Huber](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30071-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471219300717%3Fshowall%3Dtrue)


<!-- ## exercises {-}

Using the `hawaii_aquifers` data set or the `tequila_chemistry` data set, please complete the following:

1. Choose one analyte and filter the data so only the rows for that analyte are shown.

2. Choose two of the aquifers (or bottles of tequila). Are the mean abundances for your chosen analyte different in these two aquifers? Don't forget to test your data for normality and homogeneity of variance before selecting a statistical test. Use a plot to illustrate whether the means are similar or different. If you wish, you can use the figure caption to help illustrate whether means are significantly different or not. Either way, be sure to describe the tests you ran in your figure caption, and how you interpreted the results of the tests. For example, you should include a statment like "Groups sharing a letter label (a, b) are not significantly different (ANOVA/Tukey HSD tests, p > 0.05). Groups with different letters indicate a significant difference (p < 0.05).", or, if you are using asterisks then perhaps something like "Significant differences between groups are indicated by asterisks: \* = p < 0.05, \*\* = p < 0.01 (Kruskal / Dunn tests)."

3. Choose a second analyte, different from the first one you chose. Considering all the aquifers (or all bottles of tequila) in the dataset, do any of them have the same abundance of this analyte? Again, don't forget about normality and homogeneity of variance tests. Use a plot to illustrate your answer. Be sure to describe the tests you ran in your figure caption, and how you interpreted the results of the tests.

4. Repeat #3 above, but switch the type of test used (i.e. use non-parametric if you used parametric for #3 and vice-versa). Compare the *p* values and *p* groups obtained by the two methods. Use a graphic to illustrate this. Why are they different? -->

<!-- end -->

# (PART) MODELS

<!-- start numerical models -->

# numerical models {-}


<!-- explain cross validation and rmse in some sort of metrics section -->


<img src="https://thebustalab.github.io/integrated_bioanalytics/images/models.png" width="100%" style="display: block; margin: auto;" />

## model use {-}

Next on our quest to develop our abilities in analytical data exploration is modeling. We will discuss two main ways to use numerical models: for inferential uses and predictive uses.

- Inferential uses aim to understand and quantify the relationships between variables, focusing on the significance, direction, and strength of these relationships to draw conclusions about the data. When using inferential models we care a lot about the exact inner workings of the model because those inner workings are how we understand relationships between variables.

- Predictive uses are designed with the primary goal of forecasting future outcomes or behaviors based on historical data. When using predictive models we often care much less about the exact inner workings of the model and instead care about accuracy. In other words, we don't really care how the model works as long as it is accurate.

In the preceding section, we explored inferential models, highlighting the significance of grasping the quantitative aspects of the model's mechanisms in order to understand relationships between variables. Now we will look at predictive models. Unlike inferential models, predictive models are typically more complex, often making it challenging to fully comprehend their internal processes. That is okay though, because when using predictive models we usually care most about their predictive accuracy, and are willing to sacrifice our quantitative understanding of the model’s inner workings to achieve higher accuracy.

Interestingly, increasingly complex predictive models do not always have higher accuracy. If a model is too complex we say that the model is 'overfitted' which means that the model is capturing noise and random fluctuations in the input data and using those erroneous patterns in its predictions. On the other hand, if a model is not complex enough then it will not be able to capture important true patterns in the data that are required to make accurate predictions. This means that we have to build models with the right level of complexity.

To build a model with the appropriate level of complexity we usually use this process: (i) separate out 80% of our data and call it the training data, (ii) build a series of models with varying complexity using the training data, (iii) use each of the models to make predictions about the remaining 20% of the data (the testing data), (iv) whichever model has the best predictive accuracy on the remaining 20% is the model with the appropriate level of complexity.

In this course, we will build both types of models using a function called `buildModel`. To use it, we simply give it our data, and tell it which to sets of values we want to compare. To tell it what we want to compare, we need to specify (at least) two things:

- input_variables: the input variables (sometimes called "features" or "predictors") the model should use as inputs for the prediction. Depending on the model, these could be continuous or categorical variables.

- output_variable: the variable that the model should predict. Depending on the model, it could be a continuous value (regression) or a category/class (classification).

For model building, we also need to talk about handling missing data. If we have missing data in our data set, we need to  one way forward is to impute it. This means that we need to fill in the missing values with something. There are many ways to do this, but we will use the median value of each column. We can do this using the `impute` function from the `rstatix` package. Let's do that now:


``` r
any(is.na(metabolomics_data))
## [1] TRUE
```

``` r
metabolomics_data[] <- lapply(metabolomics_data, function(x) Hmisc::impute(x, median(x, na.rm = TRUE)))
any(is.na(metabolomics_data))
## [1] FALSE
```

## single linear regression {-}

We will start with some of the simplest models - linear models. There are a variety of ways to build linear models in R, but we will use `buildModel`, as mentioned above. First, we will use least squares regression to model the relationship between input and output variables. Suppose we want to know if the abundances of ADP and AMP are related in our metabolomics dataset:


``` r
ggplot(metabolomics_data) +
  geom_point(aes(x = ADP, y = AMP))
## Don't know how to automatically pick scale for object of
## type <impute>. Defaulting to continuous.
```

<img src="index_files/figure-html/unnamed-chunk-182-1.png" width="100%" style="display: block; margin: auto;" />

It looks like there might be a relationship! Let's build an inferential model to examine the details of that that relationship:


``` r
basic_regression_model <- buildModel2(
  data = metabolomics_data,
  model_type = "linear_regression",
  input_variables = "ADP",
  output_variable = "AMP"
)
names(basic_regression_model)
## [1] "model_type" "model"      "metrics"
```
The output of `buildModel` consists of three things: the model_type, the model itself, and the metric describing certain aspects of the model and/or its performance. Let's look at the model:


``` r
basic_regression_model$model
## 
## Call:
## lm(formula = formula, data = data, model = TRUE, x = TRUE, y = TRUE, 
##     qr = TRUE)
## 
## Coefficients:
## (Intercept)          ADP  
##      4.5764       0.6705
```

This is a linear model stored inside a special object type inside R called an `lm`. They can be a bit tricky to work with, but we have a way to make it easier - we'll look at that in a second. Before that, let's look at the metrics.


``` r
basic_regression_model$metrics
##               variable   value std_err        type  p_value
## 1            r_squared  0.4693      NA   statistic       NA
## 2    total_sum_squares 38.6346      NA   statistic       NA
## 3 residual_sum_squares 20.5024      NA   statistic       NA
## 4          (Intercept)  4.5764  0.9667 coefficient 8.04e-06
## 5                  ADP  0.6705  0.0747 coefficient 0.00e+00
##   p_value_adj
## 1          NA
## 2          NA
## 3          NA
## 4    8.04e-06
## 5    0.00e+00
```

It shows us the r-squared, the total and residual sum of squares, the intercept (b in y = mx + b), and the coefficient for AMP (i.e. the slope, m), as well some other things (we will talk about them in a second).

We can also use a function called `predictWithModel` to make some predictions using the model. Let's try that for ADP and AMP. What we do is give it the model, and then tell it what values we want to predict for. In this case, we want to predict the abundance of ADP for each value of AMP in our data set. We can do that like this:


``` r
AMP_values_predicted_from_ADP_values <- predictWithModel(
  data = metabolomics_data,
  model_type = "linear_regression",
  model = basic_regression_model$model
)
head(AMP_values_predicted_from_ADP_values)
##        1        2        3        4        5        6 
## 13.18461 13.35352 13.48067 13.42760 12.57818 13.18366
```

So, `predictWithModel` is using the model to predict AMP values from ADP. However, note that we have the measured AMP values in our data set. We can compare the predicted values to the measured values to see how well our model is doing. We can do that like this:


``` r
ADP_values <- metabolomics_data$ADP

predictions_from_basic_linear_model <- data.frame(
    ADP_values = ADP_values,
    AMP_values_predicted_from_ADP_values = AMP_values_predicted_from_ADP_values,
    AMP_values_measured = metabolomics_data$AMP
)

plot1 <- ggplot() +
    geom_line(
        data = predictions_from_basic_linear_model,
        aes(x = ADP_values, y = AMP_values_predicted_from_ADP_values), color = "red"
    ) +
    geom_point(
        data = predictions_from_basic_linear_model,
        aes(x = ADP_values, y = AMP_values_predicted_from_ADP_values), color = "red"
    ) +
    geom_point(
        data = metabolomics_data,
        aes(x = ADP, y = AMP), color = "blue"
    )
plot1
## Don't know how to automatically pick scale for object of
## type <impute>. Defaulting to continuous.
```

<img src="index_files/figure-html/unnamed-chunk-187-1.png" width="100%" style="display: block; margin: auto;" />

Very good. Now let's talk about evaluating the quality of our model. For this we need some means of assessing how well our line fits our data. We will use residuals - the distance between each of our points and our line.


``` r
ggplot(predictions_from_basic_linear_model) +
  geom_point(aes(x = ADP_values, y = AMP_values_measured)) +
  geom_line(aes(x = ADP_values, y = AMP_values_predicted_from_ADP_values)) +
  geom_segment(aes(x = ADP_values, y = AMP_values_measured, xend = ADP_values, yend = AMP_values_predicted_from_ADP_values))
## Don't know how to automatically pick scale for object of
## type <impute>. Defaulting to continuous.
```

<img src="index_files/figure-html/unnamed-chunk-188-1.png" width="100%" style="display: block; margin: auto;" />

We can calculate the sum of the squared residuals:


``` r
sum(
  (predictions_from_basic_linear_model$AMP_values_measured - predictions_from_basic_linear_model$AMP_values_predicted_from_ADP_values)^2
, na.rm = TRUE)
## [1] 20.50235
```

Cool! Let's call that the "residual sum of the squares". So... does that mean our model is good? I don't know. We have to compare that number to something. Let's compare it to a super simple model that is just defined by the mean y value of the input data.


``` r
ggplot(metabolomics_data) +
  geom_point(aes(x = ADP, y = AMP)) +
  geom_hline(aes(yintercept = mean(AMP, na.rm = TRUE)))
## Don't know how to automatically pick scale for object of
## type <impute>. Defaulting to continuous.
```

<img src="index_files/figure-html/unnamed-chunk-190-1.png" width="100%" style="display: block; margin: auto;" />

A pretty bad model, I agree. How much better is our linear model that the flat line model? Let's create a measure of the distance between each point and the point predicted for that same x value on the model:


``` r
ggplot(metabolomics_data) +
  geom_point(aes(x = ADP, y = AMP)) +
  geom_hline(aes(yintercept = mean(ADP, na.rm = TRUE))) +
  geom_segment(aes(x = ADP, y = AMP, xend = ADP, yend = mean(ADP, na.rm = TRUE)))
## Don't know how to automatically pick scale for object of
## type <impute>. Defaulting to continuous.
```

<img src="index_files/figure-html/unnamed-chunk-191-1.png" width="100%" style="display: block; margin: auto;" />

``` r

sum(
  (metabolomics_data$AMP - mean(metabolomics_data$AMP, na.rm = TRUE))^2
, na.rm = TRUE)
## [1] 38.63462
```

Cool. Let's call that the "total sum of the squares", and now we can compare that to our "residual sum of the squares": 


``` r
1-(20.1904/38.63)
## [1] 0.4773389
```

Alright. That is our R squared value. It is equal to 1 minus the ratio of the "residual sum of the squares" to the "total sum of the squares". You can think of the R squared value as:
- The amount of variance in the response explained by the dependent variable.
- How much better the line of best fit describes the data than the flat line.
Now, let's put it all together and make it pretty:


``` r
top <- ggplot() +
    geom_line(
        data = predictions_from_basic_linear_model,
        aes(x = ADP_values, y = AMP_values_predicted_from_ADP_values), color = "red"
    ) +
    geom_point(
        data = predictions_from_basic_linear_model,
        aes(x = ADP_values, y = AMP_values_predicted_from_ADP_values), color = "red"
    ) +
    geom_point(
        data = metabolomics_data,
        aes(x = ADP, y = AMP), color = "blue"
    ) +
    annotate(geom = "table",
      x = 10,
      y = 15,
      label = list(select(basic_regression_model$metrics, variable, value))
    ) +
    coord_cartesian(ylim = c(10,16)) +
    theme_bw()

bottom <- ggplot(predictions_from_basic_linear_model) +
  geom_col(
    aes(x = ADP_values, y = AMP_values_measured-AMP_values_predicted_from_ADP_values),
    width = 0.03, color = "black", position = "dodge", alpha = 0.5
  ) +
  theme_bw()

cowplot::plot_grid(top, bottom, ncol = 1, labels = "AUTO", rel_heights = c(2,1))
## Don't know how to automatically pick scale for object of
## type <impute>. Defaulting to continuous.
## Don't know how to automatically pick scale for object of
## type <impute>. Defaulting to continuous.
```

<img src="index_files/figure-html/unnamed-chunk-193-1.png" width="100%" style="display: block; margin: auto;" />

## multiple linear regression {-}

Cool! Now let's try a multiple linear regression model. This is the same as a simple linear regression model, but with more than one predictor variable. Simple and multiple linear regression are both statistical methods used to explore the relationship between one or more independent variables (predictor variables) and a dependent variable (outcome variable). Simple linear regression involves one independent variable to predict the value of one dependent variable, utilizing a linear equation of the form y = mx + b. Multiple linear regression extends this concept to include two or more independent variables, with a typical form of  y = m1x1 + m2x2 + ... + b, allowing for a more complex representation of relationships among variables. While simple linear regression provides a straight-line relationship between the independent and dependent variables, multiple linear regression can model a multi-dimensional plane in the variable space, providing a more nuanced understanding of how the independent variables collectively influence the dependent variable. The complexity of multiple linear regression can offer more accurate predictions and insights, especially in scenarios where variables interact or are interdependent, although it also requires a more careful consideration of assumptions and potential multicollinearity among the independent variables. Let's try it with the first 30 metabolites in our data set:


``` r

basic_regression_model <- buildModel2(
  data = metabolomics_data,
  model_type = "linear_regression",
  input_variables = "ADP",
  output_variable = "AMP"
)

multiple_regression_model <- buildModel2(
  data = metabolomics_data,
  model_type = "linear_regression",
  input_variables = colnames(metabolomics_data)[3:33],
  output_variable = "AMP"
)

ggplot() +
  geom_point(
    data = metabolomics_data,
    aes(x = ADP, y = AMP), fill = "gold", shape = 21, color = "black"
  ) +
  geom_line(aes(
    x = metabolomics_data$ADP,
    y = mean(metabolomics_data$AMP)
  ), color = "grey") +
  geom_line(aes(
    x = metabolomics_data$ADP,
    y = predictWithModel(
      data = metabolomics_data,
      model_type = "linear_regression",
      model = basic_regression_model$model
    )),
    color = "maroon", size = 1
  ) +
  geom_line(aes(
    x = metabolomics_data$ADP,
    y = predictWithModel(
      data = metabolomics_data,
      model_type = "linear_regression",
      model = multiple_regression_model$model
    )),
    color = "black", size = 1, alpha = 0.6
  ) +
  theme_bw()
## Don't know how to automatically pick scale for object of
## type <impute>. Defaulting to continuous.
```

<img src="index_files/figure-html/unnamed-chunk-194-1.png" width="100%" style="display: block; margin: auto;" />

## assessing regression models {-}

There are three aspects of a regression model that we should check on, to make sure the model isn't violating and assumptions that we make when declaring the model to be valid:

1. Residual-Based Assumptions. These checks are specifically about the behavior of residuals, which are the differences between observed values and model predictions. This group includes:

    - Linearity: Examines if residuals show a random scatter around zero, indicating a linear relationship between predictors and the response variable. Patterns or curves in this plot may indicate that the model does not adequately capture the true relationship, suggesting potential non-linearity.
    
    - Homoscedasticity (Homogeneity of Variance): Looks at the spread of residuals to confirm that variance is consistent across fitted values. In other words, the spread of the residuals should be uniform regardless of the fitted values. To assess homoscedasticity, residuals are plotted against fitted values. A horizontal band of residuals indicates that the variance is consistent, supporting the homoscedasticity assumption. Conversely, patterns such as a funnel shape, where residuals spread out or contract as fitted values increase, suggest heteroscedasticity, indicating that the variance of errors changes with the level of the independent variables.
    
    - Normality of Residuals: Assesses whether residuals follow a normal distribution, crucial for statistical inference in regression. To check this characteristic, a Quantile-Quantile (Q-Q) plot is used, where the residuals are plotted against a theoretical normal distribution. If the residuals are normally distributed, the points will align closely along a straight line. Significant deviations from this line indicate departures from normality, which may affect the reliability of statistical inferences drawn from the model.

2. Predictor Relationships: This check pertains to relationships among the predictor/input variables themselves, rather than their relationship with the response/output variable. In this case:
    
    - Collinearity: Assesses multicollinearity, which occurs when predictors are highly correlated with each other. This can lead to inflated variances of regression coefficients, making it challenging to attribute effects to individual predictors. This check helps ensure that predictors are independent enough to provide clear, interpretable results for each variable’s influence on the response.

3. Model Fit and Influence: These checks look at the overall fit of the model and assess if specific data points have undue influence on the model’s results. This group includes:

    - Posterior Predictive Check: This checks if model predictions align well with observed data, indicating a good overall fit. While not directly a residual analysis, it’s a comprehensive check for how well the model captures data patterns.

    - Influential Observations: Identifies data points that might disproportionately affect the model. High-leverage points can distort model estimates, so it’s important to verify that no single observation is overly influential.

We can assess all of these using:


``` r
multiple_regression_model <- buildModel2(
  data = metabolomics_data,
  model_type = "linear_regression",
  input_variables = colnames(metabolomics_data)[3:10],
  output_variable = "AMP"
)

check_model(multiple_regression_model$model)
```

<img src="index_files/figure-html/unnamed-chunk-195-1.png" width="100%" style="display: block; margin: auto;" />

## random forests {-}

Random forests are collections of decision trees that can be used for predicting categorical variables (i.e. a 'classification' task) and for predicting numerical variables (a 'regression' task). Random forests are built by constructing multiple decision trees, each using a randomly selected subset of the training data, to ensure diversity among the trees. At each node of each tree, a number of input variables are randomly chosen as candidates for splitting the data, introducing further randomness beyond the data sampling. Among the variables randomly selected as candidates for splitting at each node, one is chosen for that node based on a criterion, such as maximizing purity in the tree's output, or minimizing mean squared error for regression tasks, guiding the construction of a robust ensemble model. The forest's final prediction is derived either through averaging the outputs (for regression) or a majority vote (for classification).

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/random_forests.jpeg" width="100%" style="display: block; margin: auto;" />

We can use the `buildModel()` function to make a random forest model. We need to specify data, a model type, input and output variables, and in the case of a random forest model, we also need to provide a list of optimization parameters: `n_vars_tried_at_split`, `n_trees`, and `min_n_leaves`. Here is more information on those parameters:

- n_vars_at_split (often called "mtry" in other implementations): this parameter specifies the number of variables that are randomly sampled as candidate features at each split point in the construction of a tree. The main idea behind selecting a subset of features (variables) is to introduce randomness into the model, which helps in making the model more robust and less likely to overfit to the training data. By trying out different numbers of features, the model can find the right level of complexity, leading to more generalized predictions. A smaller value of n_vars_at_split increases the randomness of the forest, potentially increasing bias but decreasing variance. Conversely, a larger mtry value makes the model resemble a bagged ensemble of decision trees, potentially reducing bias but increasing variance.

- n_trees (often referred to as "num.trees" or "n_estimators" in other implementations): this parameter defines the number of trees that will be grown in the random forest. Each individual tree predicts the outcome based on the subset of features it considers, and the final prediction is typically the mode (for classification tasks) or average (for regression tasks) of all individual tree predictions. Increasing the number of trees generally improves the model's performance because it averages more predictions, which tends to reduce overfitting and makes the model more stable. However, beyond a certain point, adding more trees offers diminishing returns in terms of performance improvement and can significantly increase computational cost and memory usage without substantial gains.

- min_n_leaves (often referred to as "min_n" in other implementations, default value is 1): This parameter sets the minimum number of samples that must be present in a node for it to be split further. Increasing this value makes each tree in the random forest less complex by reducing the depth of the trees, leading to larger, more generalized leaf nodes. This can help prevent overfitting by ensuring that the trees do not grow too deep or too specific to the training data. By carefully tuning this parameter, you can strike a balance between the model's ability to capture the underlying patterns in the data and its generalization to unseen data.

`buildModel()` is configured to allow you to explore a number of settings for both n_vars_at_split and n_trees, then pick the combination with the highest predictive accuracy. In this function:

- `data` specifies the dataset to be used for model training, here metabolomics_data.
- `model_type` defines the type of model to build, with "random_forest_regression" indicating a random forest model for regression tasks.
- `input_variables` selects the features or predictors for the model, here using columns 3 to 33 from metabolomics_data as predictors.
- `output_variable` is the target variable for prediction, in this case, "AMP".

The optimization_parameters argument takes a list to define the grid of parameters for optimization, including n_vars_tried_at_split, n_trees, and min_leaf_size. The seq() function generates sequences of numbers and is used here to create ranges for each parameter:

- `n_vars_tried_at_split` = seq(1,24,3) generates a sequence for the number of variables tried at each split, starting at 1, ending at 24, in steps of 3 (e.g., 1, 4, 7, ..., 24).
- `n_trees` = seq(1,40,2) creates a sequence for the number of trees in the forest, from 1 to 40 in steps of 2.
- `min_leaf_size` = seq(1,3,1) defines the minimal size of leaf nodes, ranging from 1 to 3 in steps of 1.

This setup creates a grid of parameter combinations where each combination of n_vars_tried_at_split, n_trees, and min_leaf_size defines a unique random forest model. The function will test each combination within this grid to identify the model that performs best according to a given evaluation criterion, effectively searching through a defined parameter space to optimize the random forest's performance. This approach allows for a systematic exploration of how different configurations affect the model's ability to predict the output variable, enabling the selection of the most effective model configuration based on the dataset and task at hand.
































































































































































































