--- 
title: "Integrated Bioanalytics"
author: "Lucas Busta and members of the Busta lab"
date: "2026-03-02"
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

________________________________________________________________________________________________
________________________________________________________________________________________________
________________________________________________________________________________________________

# (PART) GETTING STARTED 


# overview {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/chemometrics.jpeg" width="100%" style="display: block; margin: auto;" />

In bioanalytical science, we separate, identify, and quantify matter - be it DNA, RNA, proteins, small molecules, or even atoms. To connect our data with the world around us and answer scientific questions, multiple chemical entities must be separated, quantified, and identified. As our ability to collect analytical data expands, so must our ability to effectively analyze that data - whether its 10 data points or 10,000.

This book first covers data analysis in R. We will first look at tools for hypothesis generation, including: (i) encoding variables in visual representations of data and (ii) summarizing and providing overviews of large data set. We will then turn to evaluating hypothesese with data by looking at statistical tests and models. Finally, we will look at how to communicate our results in a clear and effective way. These techniques will also allow us to answer common quesions we may have about our data: "Which of my samples are most closely related?", "Which analytes are driving differences among my samples?", "Do my samples fall into definable clusters?", "Are any of my variables related?", and "Are any of these distributions different?".

Let's get started!


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

<img src="index_files/figure-html/unnamed-chunk-62-1.png" width="100%" style="display: block; margin: auto;" />

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
________________________________________________________________________________________________
________________________________________________________________________________________________
________________________________________________________________________________________________

# (PART) DATA VISUALIZATION


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

`filter(<data>, <variable> %in% c(18, 19, 20))` ## equal to 18 or 19 or 20

## ggplot & geoms {-}

Now we have a nice, small table that we can use to practice data visualization. For visualization, we're going to use `ggplot2` - a powerful set of commands for plot generation. 

There are three steps to setting up a ggplot:

1. **Define the data you want to use.**

We do this using the ggplot function's data argument. When we run that line, it just shows a grey plot space. Why is this? It's because all we've done is told ggplot that (i) we want to make a plot and (ii) what data should be used. We haven't explained how to represent features of the data using ink.


``` r
ggplot(data = algae_data_small)
```

<img src="index_files/figure-html/unnamed-chunk-85-1.png" width="50%" style="display: block; margin: auto;" />

2. **Define how your variables map onto the axes.**

This is called aesthetic mapping and is done with the `aes()` function. `aes()` should be placed inside the `ggplot` command. Now when we run it, we get our axes!


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance))
```

<img src="index_files/figure-html/unnamed-chunk-86-1.png" width="50%" style="display: block; margin: auto;" />

3. **Use geometric shapes to represent other variables in your data.**

Map your variables onto the geometric features of the shapes. To define which shape should be used, use a `geom_*` command. Some options are, for example, `geom_point()`, `geom_boxplot()`, and `geom_violin()`. These functions should be added to your plot using the `+` sign. We can use a new line to keep the code from getting too wide, just make sure the `+` sign is at the end fo the top line. Let's try it:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_point()
```

<img src="index_files/figure-html/unnamed-chunk-87-1.png" width="50%" style="display: block; margin: auto;" />

In the same way that we mapped variables in our dataset to the plot axes, we can map variables in the dataset to the geometric features of the shapes we are using to represent our data. For this, again, use `aes()` to map your variables onto the geometric features of the shapes:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_point(aes(color = harvesting_regime))
```

<img src="index_files/figure-html/unnamed-chunk-88-1.png" width="50%" style="display: block; margin: auto;" />

In the plot above, the points are a bit small, how could we fix that? We can modify the features of the shapes by adding additional arguments to the `geom_*()` functions. To change the size of the points created by the `geom_point()` function, this means that we need to add the `size = ` argument. IMPORTANT! Please note that when we map a feature of a shape to a *variable* in our data(as we did with color/harvesting regime, above) then it goes *inside* aes(). In contrast, when we map a feature of a shape to a *constant*, it goes *outside* aes(). Here's an example:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_point(aes(color = harvesting_regime), size = 5)
```

<img src="index_files/figure-html/unnamed-chunk-89-1.png" width="50%" style="display: block; margin: auto;" />

One powerful aspect of `ggplot` is the ability to quickly change mappings to see if alternative plots are more effective at bringing out the trends in the data. For example, we could modify the plot above by switching how harvesting_regime is mapped:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_point(aes(size = harvesting_regime), color = "black")
```

<img src="index_files/figure-html/unnamed-chunk-90-1.png" width="50%" style="display: block; margin: auto;" />

** Important note: Inside the `aes()` function, map aesthetics (the features of the geom's shape) to a *variable*. Outside the `aes()` function, map aesthetics to *constants*. You can see this in the above two plots - in the first one, color is inside `aes()` and mapped to the variable called harvesting_regime, while size is outside the `aes()` call and is set to the constant 5. In the second plot, the situation is reversed, with size being inside the `aes()` function and mapped to the variable harvesting_regime, while color is outside the `aes()` call and is mapped to the constant "black".

We can also stack geoms on top of one another by using multiple `+` signs. We also don't have to assign the same mappings to each geom.


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_violin() +
  geom_point(aes(color = harvesting_regime), size = 5)
```

<img src="index_files/figure-html/unnamed-chunk-91-1.png" width="50%" style="display: block; margin: auto;" />

As you can probably guess right now, there are lots of mappings that can be done, and lots of different ways to look at the same data!


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_violin(aes(fill = algae_strain)) +
  geom_point(aes(color = harvesting_regime, size = replicate))
```

<img src="index_files/figure-html/unnamed-chunk-92-1.png" width="50%" style="display: block; margin: auto;" />


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_boxplot()
```

<img src="index_files/figure-html/unnamed-chunk-93-1.png" width="50%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-100-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-117-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-162-1.png" width="100%" style="display: block; margin: auto;" />

Also, please be aware of `geom_tile()`, which is nice for situations with two discrete variables and one continuous variable. `geom_tile()` makes what are often referred to as heat maps. Note that `geom_tile()` is somewhat similar to `geom_point(shape = 21)`, in that it has both `fill` and `color` aesthetics that control the fill color and the border color, respectively.


``` r
ggplot(
  data = filter(algae_data, harvesting_regime == "Heavy"),
  aes(x = algae_strain, y = chemical_species)
) + 
  geom_tile(aes(fill = abundance), color = "black", size = 1)
```

<img src="index_files/figure-html/unnamed-chunk-163-1.png" width="100%" style="display: block; margin: auto;" />

These examples should illustrate that there is, to some degree, correspondence between the type of data you are interested in plotting (number of discrete and continuous variables) and the types of geoms that can effectively be used to represent the data.

## facets {-}

As alluded to in Exercises 1, it is possible to map variables in your dataset to more than the geometric features of shapes (i.e. geoms). One very common way of doing this is with facets. Faceting creates small multiples of your plot, each of which shows a different subset of your data based on a categorical variable of your choice. Let's check it out.

Here, we can facet in the horizontal direction:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(.~replicate)
```

<img src="index_files/figure-html/unnamed-chunk-164-1.png" width="100%" style="display: block; margin: auto;" />

We can facet in the vertical direction:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(replicate~.)
```

<img src="index_files/figure-html/unnamed-chunk-165-1.png" width="100%" style="display: block; margin: auto;" />

And we can do both at the same time:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(harvesting_regime~replicate)
```

<img src="index_files/figure-html/unnamed-chunk-166-1.png" width="100%" style="display: block; margin: auto;" />

Faceting is a great way to describe more variation in your plot without having to make your geoms more complicated. For situations where you need to generate lots and lots of facets, consider `facet_wrap` instead of `facet_grid`:



``` r
ggplot(data = algae_data, aes(x = replicate, y = algae_strain)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_wrap(chemical_species~.)
```

<img src="index_files/figure-html/unnamed-chunk-167-1.png" width="100%" style="display: block; margin: auto;" />

## scales {-}

Every time you define an aesthetic mapping (e.g. aes(x = algae_strain)), you are defining a new scale that is added to your plot. You can control these scales using the `scale_*` family of commands. Consider our faceting example above. In it, we use `geom_tile(aes(fill = abundance))` to map the abundance variable to the fill aesthetic of the tiles. This creates a scale called `fill` that we can adjust using `scale_fill_*`. In this case, fill is mapped to a continuous variable and so the fill scale is a color gradient. Therefore, `scale_fill_gradient()` is the command we need to change it. Remember that you could always type `?scale_fill_` into the console and it will help you find relevant help topics that will provide more detail. Another option is to google: "How to modify color scale ggplot geom_tile", which will undoubtedly turn up a wealth of help.


``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(harvesting_regime~replicate) +
  scale_fill_gradient(low = "white", high = "black") +
  theme_classic()
```

<img src="index_files/figure-html/unnamed-chunk-168-1.png" width="100%" style="display: block; margin: auto;" />

One particularly useful type of scale are the color scales provided by RColorBrewer:


``` r
display.brewer.all()
```

<img src="index_files/figure-html/unnamed-chunk-169-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ggplot(mtcars) +
  geom_point(
    aes(x = mpg, y = factor(cyl), fill = factor(carb)), 
    shape = 21, size = 6
  ) +
  scale_fill_brewer(palette = "Set1")
```

<img src="index_files/figure-html/unnamed-chunk-170-1.png" width="100%" style="display: block; margin: auto;" />
  
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

<img src="index_files/figure-html/unnamed-chunk-171-1.png" width="100%" style="display: block; margin: auto;" />


``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point() +
  theme_dark()
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-172-1.png" width="100%" style="display: block; margin: auto;" />
  

``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point() +
  theme_void()
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-173-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-174-1.png" width="100%" style="display: block; margin: auto;" />

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
<img src="index_files/figure-html/unnamed-chunk-175-1.png" alt="Vapor pressure as a function of boiling point. A scatter plot with trendline showing the vapor pressure of thirty-two solvents (y-axis) a as a function of their boiling points (x-axis). Each point represents the boiling point and vapor pressure of one solvent. Data are from the 'solvents' dataset used in UMD CHEM5725." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-175)Vapor pressure as a function of boiling point. A scatter plot with trendline showing the vapor pressure of thirty-two solvents (y-axis) a as a function of their boiling points (x-axis). Each point represents the boiling point and vapor pressure of one solvent. Data are from the 'solvents' dataset used in UMD CHEM5725.</p>
</div>

In some cases, the following diagram illustrates a useful way to think about the `ggplot()` / `geom_*()` / `scale_*()` / `theme_*()` situation. It shows how we use these things together to achieve a sharp-looking plot:

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/what_is_ggplot.jpeg" width="100%" style="display: block; margin: auto;" />

## subplots {-}

We can make subplots using the `plot_grid()` function from the `cowplot` package, which comes with the `source()` command. Let's see:


``` r
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

<img src="index_files/figure-html/unnamed-chunk-177-1.png" width="100%" style="display: block; margin: auto;" />

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

There is a [handy cheat sheet](https://thebustalab.github.io/integrated_bioanalytics/images/ggplot2_geoms.pdf) that can help you identify the right geom for your situation. Please keep this cheat sheet in mind for your future plotting needs...

For additional explanations of ggplot2: [ggplot2-book](https://ggplot2-book.org/).

Check out some of the incredible geoms that are easy to access using R and ggplot2: [R Graph Gallery](https://r-graph-gallery.com/). Use these to make your figures attractive and easy to interpret!

For a challenge, try implementing these awesome color scales: [Famous R Color Palettes](https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/). Note that some of these are optimized for colorblind individuals and that other are optimized for continuous hue gradients, etc.

For a list of data visualization sins: [Friends Don't Let Friends](https://github.com/cxli233/FriendsDontLetFriends). Some interesting things in here!

For more information on data visualization and graphics theory, check out the works by Edward Tufte: [Edward Tufte](https://www.edwardtufte.com/tufte/). A digital text that covers similar topics is here: [Look At Data] (https://socviz.co/lookatdata.html).

Some examples of award winning data visualization: [Information Is Beautiful Awards](https://www.informationisbeautifulawards.com/showcase?award=2019&type=awards) and [Data Vis Inspiration](https://www.dataviz-inspiration.com/).

Additional color palettes: [MetBrewer](https://github.com/BlakeRMills/MetBrewer) and [Paletteer](https://github.com/EmilHvitfeldt/paletteer).

<!-- end -->

<!-- start data visualization III -->























































































