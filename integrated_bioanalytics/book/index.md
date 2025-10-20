--- 
title: "Integrated Bioanalytics"
author: "Lucas Busta and members of the Busta lab"
date: "2025-10-20"
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

<img src="index_files/figure-html/unnamed-chunk-64-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-87-1.png" width="50%" style="display: block; margin: auto;" />

2. **Define how your variables map onto the axes.**

This is called aesthetic mapping and is done with the `aes()` function. `aes()` should be placed inside the `ggplot` command. Now when we run it, we get our axes!


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance))
```

<img src="index_files/figure-html/unnamed-chunk-88-1.png" width="50%" style="display: block; margin: auto;" />

3. **Use geometric shapes to represent other variables in your data.**

Map your variables onto the geometric features of the shapes. To define which shape should be used, use a `geom_*` command. Some options are, for example, `geom_point()`, `geom_boxplot()`, and `geom_violin()`. These functions should be added to your plot using the `+` sign. We can use a new line to keep the code from getting too wide, just make sure the `+` sign is at the end fo the top line. Let's try it:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_point()
```

<img src="index_files/figure-html/unnamed-chunk-89-1.png" width="50%" style="display: block; margin: auto;" />

In the same way that we mapped variables in our dataset to the plot axes, we can map variables in the dataset to the geometric features of the shapes we are using to represent our data. For this, again, use `aes()` to map your variables onto the geometric features of the shapes:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_point(aes(color = harvesting_regime))
```

<img src="index_files/figure-html/unnamed-chunk-90-1.png" width="50%" style="display: block; margin: auto;" />

In the plot above, the points are a bit small, how could we fix that? We can modify the features of the shapes by adding additional arguments to the `geom_*()` functions. To change the size of the points created by the `geom_point()` function, this means that we need to add the `size = ` argument. IMPORTANT! Please note that when we map a feature of a shape to a *variable* in our data(as we did with color/harvesting regime, above) then it goes *inside* aes(). In contrast, when we map a feature of a shape to a *constant*, it goes *outside* aes(). Here's an example:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_point(aes(color = harvesting_regime), size = 5)
```

<img src="index_files/figure-html/unnamed-chunk-91-1.png" width="50%" style="display: block; margin: auto;" />

One powerful aspect of `ggplot` is the ability to quickly change mappings to see if alternative plots are more effective at bringing out the trends in the data. For example, we could modify the plot above by switching how harvesting_regime is mapped:


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_point(aes(size = harvesting_regime), color = "black")
```

<img src="index_files/figure-html/unnamed-chunk-92-1.png" width="50%" style="display: block; margin: auto;" />

** Important note: Inside the `aes()` function, map aesthetics (the features of the geom's shape) to a *variable*. Outside the `aes()` function, map aesthetics to *constants*. You can see this in the above two plots - in the first one, color is inside `aes()` and mapped to the variable called harvesting_regime, while size is outside the `aes()` call and is set to the constant 5. In the second plot, the situation is reversed, with size being inside the `aes()` function and mapped to the variable harvesting_regime, while color is outside the `aes()` call and is mapped to the constant "black".

We can also stack geoms on top of one another by using multiple `+` signs. We also don't have to assign the same mappings to each geom.


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) + 
  geom_violin() +
  geom_point(aes(color = harvesting_regime), size = 5)
```

<img src="index_files/figure-html/unnamed-chunk-93-1.png" width="50%" style="display: block; margin: auto;" />

As you can probably guess right now, there are lots of mappings that can be done, and lots of different ways to look at the same data!


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_violin(aes(fill = algae_strain)) +
  geom_point(aes(color = harvesting_regime, size = replicate))
```

<img src="index_files/figure-html/unnamed-chunk-94-1.png" width="50%" style="display: block; margin: auto;" />


``` r
ggplot(data = algae_data_small, aes(x = algae_strain, y = abundance)) +
  geom_boxplot()
```

<img src="index_files/figure-html/unnamed-chunk-95-1.png" width="50%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-102-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-119-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-164-1.png" width="100%" style="display: block; margin: auto;" />

Also, please be aware of `geom_tile()`, which is nice for situations with two discrete variables and one continuous variable. `geom_tile()` makes what are often referred to as heat maps. Note that `geom_tile()` is somewhat similar to `geom_point(shape = 21)`, in that it has both `fill` and `color` aesthetics that control the fill color and the border color, respectively.


``` r
ggplot(
  data = filter(algae_data, harvesting_regime == "Heavy"),
  aes(x = algae_strain, y = chemical_species)
) + 
  geom_tile(aes(fill = abundance), color = "black", size = 1)
```

<img src="index_files/figure-html/unnamed-chunk-165-1.png" width="100%" style="display: block; margin: auto;" />

These examples should illustrate that there is, to some degree, correspondence between the type of data you are interested in plotting (number of discrete and continuous variables) and the types of geoms that can effectively be used to represent the data.

## facets {-}

As alluded to in Exercises 1, it is possible to map variables in your dataset to more than the geometric features of shapes (i.e. geoms). One very common way of doing this is with facets. Faceting creates small multiples of your plot, each of which shows a different subset of your data based on a categorical variable of your choice. Let's check it out.

Here, we can facet in the horizontal direction:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(.~replicate)
```

<img src="index_files/figure-html/unnamed-chunk-166-1.png" width="100%" style="display: block; margin: auto;" />

We can facet in the vertical direction:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(replicate~.)
```

<img src="index_files/figure-html/unnamed-chunk-167-1.png" width="100%" style="display: block; margin: auto;" />

And we can do both at the same time:

``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(harvesting_regime~replicate)
```

<img src="index_files/figure-html/unnamed-chunk-168-1.png" width="100%" style="display: block; margin: auto;" />

Faceting is a great way to describe more variation in your plot without having to make your geoms more complicated. For situations where you need to generate lots and lots of facets, consider `facet_wrap` instead of `facet_grid`:



``` r
ggplot(data = algae_data, aes(x = replicate, y = algae_strain)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_wrap(chemical_species~.)
```

<img src="index_files/figure-html/unnamed-chunk-169-1.png" width="100%" style="display: block; margin: auto;" />

## scales {-}

Every time you define an aesthetic mapping (e.g. aes(x = algae_strain)), you are defining a new scale that is added to your plot. You can control these scales using the `scale_*` family of commands. Consider our faceting example above. In it, we use `geom_tile(aes(fill = abundance))` to map the abundance variable to the fill aesthetic of the tiles. This creates a scale called `fill` that we can adjust using `scale_fill_*`. In this case, fill is mapped to a continuous variable and so the fill scale is a color gradient. Therefore, `scale_fill_gradient()` is the command we need to change it. Remember that you could always type `?scale_fill_` into the console and it will help you find relevant help topics that will provide more detail. Another option is to google: "How to modify color scale ggplot geom_tile", which will undoubtedly turn up a wealth of help.


``` r
ggplot(data = algae_data, aes(x = algae_strain, y = chemical_species)) + 
  geom_tile(aes(fill = abundance), color = "black") + 
  facet_grid(harvesting_regime~replicate) +
  scale_fill_gradient(low = "white", high = "black") +
  theme_classic()
```

<img src="index_files/figure-html/unnamed-chunk-170-1.png" width="100%" style="display: block; margin: auto;" />

One particularly useful type of scale are the color scales provided by RColorBrewer:


``` r
display.brewer.all()
```

<img src="index_files/figure-html/unnamed-chunk-171-1.png" width="100%" style="display: block; margin: auto;" />

``` r
ggplot(mtcars) +
  geom_point(
    aes(x = mpg, y = factor(cyl), fill = factor(carb)), 
    shape = 21, size = 6
  ) +
  scale_fill_brewer(palette = "Set1")
```

<img src="index_files/figure-html/unnamed-chunk-172-1.png" width="100%" style="display: block; margin: auto;" />
  
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

<img src="index_files/figure-html/unnamed-chunk-173-1.png" width="100%" style="display: block; margin: auto;" />


``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point() +
  theme_dark()
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-174-1.png" width="100%" style="display: block; margin: auto;" />
  

``` r
ggplot(data = solvents, aes(x = boiling_point, y = vapor_pressure)) + 
  geom_smooth() +
  geom_point() +
  theme_void()
## `geom_smooth()` using method = 'loess' and formula = 'y ~
## x'
```

<img src="index_files/figure-html/unnamed-chunk-175-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-176-1.png" width="100%" style="display: block; margin: auto;" />

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
<img src="index_files/figure-html/unnamed-chunk-177-1.png" alt="Vapor pressure as a function of boiling point. A scatter plot with trendline showing the vapor pressure of thirty-two solvents (y-axis) a as a function of their boiling points (x-axis). Each point represents the boiling point and vapor pressure of one solvent. Data are from the 'solvents' dataset used in UMD CHEM5725." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-177)Vapor pressure as a function of boiling point. A scatter plot with trendline showing the vapor pressure of thirty-two solvents (y-axis) a as a function of their boiling points (x-axis). Each point represents the boiling point and vapor pressure of one solvent. Data are from the 'solvents' dataset used in UMD CHEM5725.</p>
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

<img src="index_files/figure-html/unnamed-chunk-179-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-200-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-201-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-202-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-203-1.png" width="100%" style="display: block; margin: auto;" />


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

<img src="index_files/figure-html/unnamed-chunk-204-1.png" width="100%" style="display: block; margin: auto;" />


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

<img src="index_files/figure-html/unnamed-chunk-209-1.png" width="100%" style="display: block; margin: auto;" />

Note that we can use `coord_map()` to do some pretty cool things!


``` r
ggplot(map_data("world")) +
  geom_point(aes(x = long, y = lat, color = group), size = 0.5) +
  theme_void() +
  coord_map(projection = "albers", lat0 = 39, lat1 = 45)
```

<img src="index_files/figure-html/unnamed-chunk-210-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-211-1.png" width="100%" style="display: block; margin: auto;" />

### maps with plots {-}

Please note that the Great Lakes are in map_data()!


``` r
filter(map_data("lakes"), region == "Great Lakes", subregion == "Superior") %>%
    ggplot() +
      geom_path(aes(x = long, y = lat)) +
      coord_map() +
      theme_minimal()
```

<img src="index_files/figure-html/unnamed-chunk-212-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-213-1.png" width="100%" style="display: block; margin: auto;" />

Now we could add some data. We could do something simple like plot total abundances as the size of a point:


``` r
lake_superior_PFAS <- readMonolist("/Users/bust0037/Documents/Websites/pfas_data_private.csv")
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

<img src="index_files/figure-html/unnamed-chunk-214-1.png" width="100%" style="display: block; margin: auto;" />

Or we could do something more sophisticated like add pie charts at each point:



``` r
lake_superior_PFAS <- readMonolist("/Users/bust0037/Documents/Websites/pfas_data_private.csv")

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

<img src="index_files/figure-html/unnamed-chunk-215-1.png" width="100%" style="display: block; margin: auto;" />

You can also access a high resolution shoreline dataset for Lake Superior directly from the source() command as `lake_superior_shoreline`:


``` r
shore <- readMonolist("/Users/bust0037/Documents/Websites/thebustalab.github.io/phylochemistry/sample_data/lake_superior_shoreline.csv")

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

<img src="index_files/figure-html/unnamed-chunk-216-1.png" width="100%" style="display: block; margin: auto;" />

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

________________________________________________________________________________________________
________________________________________________________________________________________________
________________________________________________________________________________________________

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

# To look at a single row (the second row)
head(alaska_lake_data[2,])
## # A tibble: 1 × 7
##   lake  park  water_temp    pH element mg_per_L element_type
##   <chr> <chr>      <dbl> <dbl> <chr>      <dbl> <chr>       
## 1 Devi… BELA        6.46  7.69 N          0.028 bound

# To look at select rows:
head(alaska_lake_data[2:5,])
## # A tibble: 4 × 7
##   lake  park  water_temp    pH element mg_per_L element_type
##   <chr> <chr>      <dbl> <dbl> <chr>      <dbl> <chr>       
## 1 Devi… BELA        6.46  7.69 N          0.028 bound       
## 2 Devi… BELA        6.46  7.69 P          0     bound       
## 3 Devi… BELA        6.46  7.69 Cl        10.4   free        
## 4 Devi… BELA        6.46  7.69 S          0.62  free

# To look at just a single column, by name
head(alaska_lake_data$pH)
## [1] 7.69 7.69 7.69 7.69 7.69 7.69

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

<img src="index_files/figure-html/unnamed-chunk-241-1.png" width="100%" style="display: block; margin: auto;" />

However, as our analyses get more complex, the code can get long and hard to read. We're going to use the pipe `%>%` to help us with this. Check it out:


``` r
alaska_lake_data %>%
  filter(park == "BELA") %>%
  ggplot(aes(x = pH, y = lake)) + geom_col()
```

<img src="index_files/figure-html/unnamed-chunk-242-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-249-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-251-1.png" width="100%" style="display: block; margin: auto;" />

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

<img src="index_files/figure-html/unnamed-chunk-252-1.png" width="100%" style="display: block; margin: auto;" />

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
##  6 White_Fish_Lake     BELA  S           0.18      9.37
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

- [Tidy Data Tutor](https://tidydatatutor.com/vis.html). This interactive site animates what happens to a data frame as you apply verbs like `select()`, `filter()`, or `pivot_longer()`, making it easier to see how each operation reshapes your tables.

- [tidyverse cheat sheet](https://posit.co/resources/cheatsheets/#r-programming). Posit maintains printable cheat sheets for dplyr, tidyr, and readr that summarize the verbs we used here; they’re handy for quick lookups when you forget an argument name or want to explore related functions.

- [tidyr and dplyr vignettes](https://tidyr.tidyverse.org/articles/tidy-data.html). The official package articles dig into the philosophy behind tidy data, the rationale for pivoting functions, and worked examples that go deeper than our walkthrough.

- [R for Data Science, Chapter 5–10](https://r4ds.hadley.nz/transform). Hadley Wickham and Mine Çetinkaya-Rundel’s open textbook expands on wrangling, reshaping, and summarizing with narratives, exercises, and datasets that reinforce the patterns we practiced.



# hierarchical clustering {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/clustering.png" width="100%" style="display: block; margin: auto;" />

"Which of my samples are most closely related?"

## {-}

## clustering

So far we have been looking at how to plot raw data, summarize data, and reduce a data set's dimensionality. It's time to look at how to identify relationships between the samples in our data sets. For example: in the Alaska lakes dataset, which lake is most similar, chemically speaking, to Lake Narvakrak? Answering this requires calculating numeric distances between samples based on their chemical properties. For this, the first thing we need is a distance matrix:

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/dist_matrix.jpg" width="100%" style="display: block; margin: auto;" />

Please note that we can get distance matrices directly from `runMatrixAnalyses` by specifying `analysis = "dist"`:


``` r
alaska_lake_data %>%
    select(-element_type) %>%
    pivot_wider(names_from = "element", values_from = "mg_per_L") -> alaska_lake_data_wide

dist <- runMatrixAnalyses(
    data = alaska_lake_data_wide,
    analysis = c("dist"),
    columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:15],
    columns_w_sample_ID_info = c("lake", "park")
)
## Replacing NAs in your data with mean

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

There is more that we can do with distance matrices though, lots more. Let's start by looking at an example of hierarchical clustering. For this, we just need to tell `runMatrixAnalyses()` to use `analysis = "hclust"`: 


``` r
AK_lakes_clustered <- runMatrixAnalyses(
    data = alaska_lake_data_wide,
    analysis = c("hclust"),
    columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:15],
    columns_w_sample_ID_info = c("lake", "park")
)
AK_lakes_clustered
## # A tibble: 38 × 26
##    sample_unique_ID   lake  park  parent  node branch.length
##    <chr>              <chr> <chr>  <int> <int>         <dbl>
##  1 Devil_Mountain_La… Devi… BELA      23     1         7.81 
##  2 Imuruk_Lake_BELA   Imur… BELA      25     2         6.01 
##  3 Kuzitrin_Lake_BELA Kuzi… BELA      24     3         3.27 
##  4 Lava_Lake_BELA     Lava… BELA      32     4         3.14 
##  5 North_Killeak_Lak… Nort… BELA      38     5       256.   
##  6 White_Fish_Lake_B… Whit… BELA      38     6         0.828
##  7 Iniakuk_Lake_GAAR  Inia… GAAR      36     7         2.59 
##  8 Kurupa_Lake_GAAR   Kuru… GAAR      30     8         6.00 
##  9 Lake_Matcharak_GA… Lake… GAAR      35     9         2.09 
## 10 Lake_Selby_GAAR    Lake… GAAR      29    10         3.49 
## # ℹ 28 more rows
## # ℹ 20 more variables: label <chr>, isTip <lgl>, x <dbl>,
## #   y <dbl>, branch <dbl>, angle <dbl>, bootstrap <lgl>,
## #   water_temp <dbl>, pH <dbl>, C <dbl>, N <dbl>, P <dbl>,
## #   Cl <dbl>, S <dbl>, F <dbl>, Br <dbl>, Na <dbl>,
## #   K <dbl>, Ca <dbl>, Mg <dbl>
```

It works! Now we can plot our cluster diagram with a ggplot add-on called ggtree. We've seen that ggplot takes a "data" argument (i.e. `ggplot(data = <some_data>) + geom_*()` etc.). In contrast, ggtree takes an argument called `tr`, though if you're using the `runMatrixAnalysis()` function, you can treat these two (`data` and `tr`) the same, so, use: `ggtree(tr = <output_from_runMatrixAnalyses>) + geom_*()` etc.

Note that `ggtree` also comes with several great new geoms: `geom_tiplab()` and `geom_tippoint()`. Let's try those out:


``` r
library(ggtree)
AK_lakes_clustered %>%
ggtree() +
  geom_tiplab() +
  geom_tippoint() +
  theme_classic() +
  scale_x_continuous(limits = c(0,700))
```

<img src="index_files/figure-html/unnamed-chunk-295-1.png" width="100%" style="display: block; margin: auto;" />

Cool! Though that plot could use some tweaking... let's try:


``` r
AK_lakes_clustered %>%
ggtree() +
    geom_tiplab(aes(label = lake), offset = 10, align = TRUE) +
    geom_tippoint(shape = 21, aes(fill = park), size = 4) +
    scale_x_continuous(limits = c(0,600)) +
    scale_fill_brewer(palette = "Set1") +
    # theme_classic() +
    theme(
      legend.position = c(0.4,0.2)
    )
```

<img src="index_files/figure-html/unnamed-chunk-296-1.png" width="100%" style="display: block; margin: auto;" />

Very nice! Since North Killeak and White Fish are so different from the others, we could re-analyze the data with those two removed:


``` r
alaska_lake_data_wide %>%
  filter(!lake %in% c("North_Killeak_Lake","White_Fish_Lake")) -> alaska_lake_data_wide_filtered

runMatrixAnalyses(
    data = alaska_lake_data_wide_filtered,
    analysis = c("hclust"),
    columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide_filtered)[3:15],
    columns_w_sample_ID_info = c("lake", "park")
) %>%
ggtree() +
    geom_tiplab(aes(label = lake), offset = 10, align = TRUE) +
    geom_tippoint(shape = 21, aes(fill = park), size = 4) +
    scale_x_continuous(limits = c(0,100)) +
    scale_fill_brewer(palette = "Set1") +
    # theme_classic() +
    theme(
      legend.position = c(0.05,0.9)
    )
## Replacing NAs in your data with mean
```

<img src="index_files/figure-html/unnamed-chunk-297-1.png" width="100%" style="display: block; margin: auto;" />

## Annotating trees {-}

Overlaying sample traits on a ggtree-based plot is straightforward when we combine `ggtree` with `ggplot2`. We begin by running the hierarchical clustering analysis and keeping its output for later plotting.


``` r
hclust_out <- runMatrixAnalyses(
  data = chemical_blooms,
  analysis = c("hclust"),
  columns_w_values_for_single_analyte = colnames(chemical_blooms)[2:10],
  columns_w_sample_ID_info = "label"
)
## ! The tree contained negative edge lengths. If you want to
## ignore the edges, you can set
## `options(ignore.negative.edge=TRUE)`, then re-run ggtree.
```

The object returned by `runMatrixAnalyses()` already contains the branch coordinates, so we can pass it directly to `ggtree()` and add tip labels while keeping the tree readable. Note that we should deliberately control the y-axis here, that will be key for aligning the tree with other plots later.


``` r
tree_plot <- ggtree(hclust_out) +
  geom_tiplab(size = 2, align = TRUE) +
  scale_x_continuous(limits = c(0, 300)) +
  scale_y_continuous(limits = c(0, 80)) +
  theme_classic()
## Scale for y is already present.
## Adding another scale for y, which will replace the existing
## scale.
tree_plot
```

<img src="index_files/figure-html/unnamed-chunk-299-1.png" width="100%" style="display: block; margin: auto;" />

Next, reshape the tip-level measurements to long form so each chemical becomes its own column of tiles. Because we reuse the `y` coordinate supplied by `ggtree`, the tiles inherit the same vertical order as the tips in the tree. Note that we remove the other columns in the hclust output for simplicity - they are only needed if we want to draw the full tree. Note that we also control the y-axis here to make sure it has the same bounds (limits) as the tree we made previously.


``` r
heat_plot <- hclust_out %>%
  filter(isTip) %>%
  select(-parent, -node, -branch.length, -label, -isTip, -x, -branch, -angle, -bootstrap) %>%
  pivot_longer(cols = 3:11, names_to = "chemical", values_to = "abundance") %>%
  ggplot(aes(x = chemical, y = y, fill = abundance)) +
  geom_tile() +
  scale_y_continuous(limits = c(0, 80)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
heat_plot
```

<img src="index_files/figure-html/unnamed-chunk-300-1.png" width="100%" style="display: block; margin: auto;" />

With matching y scales, `plot_grid()` can align the tree and the heat map so the tiles line up with the corresponding samples. Using `align = "h"` snaps them together horizontally, and `axis = "tb"` keeps the panel heights consistent.


``` r
plot_grid(tree_plot, heat_plot, axis = "tb", align = "h")
```

<img src="index_files/figure-html/unnamed-chunk-301-1.png" width="100%" style="display: block; margin: auto;" />

Note: if we were to instead build the heat map directly from the raw `chemical_blooms` table, the rows fall back to their alphabetical order and the heat map no longer matches the dendrogram ordering:


``` r
chemical_blooms %>%
  pivot_longer(cols = 2:10, names_to = "chemical", values_to = "abundance") %>%
  ggplot(aes(x = chemical, y = label, fill = abundance)) +
  geom_tile()
```

<img src="index_files/figure-html/unnamed-chunk-302-1.png" width="100%" style="display: block; margin: auto;" />


## further reading {-}

- [annotating phylogenetic trees with ggtree](https://yulab-smu.top/treedata-book/chapter10.html). This free, book-length reference from the ggtree authors walks through annotating dendrograms and phylogenetic trees in R, explaining how to layer metadata, add tip labels, and customize themes—perfect for polishing the cluster plots we generated in this chapter.

- [hierarchical clustering in R](https://uc-r.github.io/hc_clustering). UC Business Analytics’ tutorial provides a gentle, R-focused walkthrough of agglomerative clustering, covering distance matrices, linkage choices, dendrogram interpretation, and cluster cutting with reproducible code examples. Note that this text uses base R instead of our in-class functions.

<!-- ## exercises {-} -->


<!-- For this set of exercises, please use `runMatrixAnalyses()` to run and visualize a hierarchical cluster analysis with each of the main datasets that we have worked with so far, except for NY_trees. This means: `algae_data` (which algae strains are most similar to each other?), `alaska_lake_data` (which lakes are most similar to each other?). and `solvents` (which solvents are most similar to each other?). It also means you should use the periodic table (which elements are most similar to each other?), though please don't use the whole periodic table, rather, use `periodic_table_subset`. Please also conduct a heirarchical clustering analysis for a dataset of your own choice that is not provided by the `source()` code. For each of these, create (i) a tree diagram that shows how the "samples" in each data set are related to each other based on the numerical data associated with them, (ii) a caption for each diagram, and (iii) describe, in two or so sentences, an interesting trend you see in the diagram. You can ignore columns that contain categorical data, or you can list those columns as "additional_analyte_info". -->

<!-- For this assignment, you may again find the `colnames()` function and square bracket-subsetting useful. It will list all or a subset of the column names in a dataset for you. For example: -->

<!-- ```{r} -->
<!-- colnames(solvents) -->

<!-- colnames(solvents)[1:3] -->

<!-- colnames(solvents)[c(1,5,7)] -->
<!-- ``` -->


# dimensional reduction {-}

<div class="figure" style="text-align: center">
<img src="https://thebustalab.github.io/integrated_bioanalytics/images/dimensionality.png" alt="Overview of dimensional reduction. The schematic shows how high-dimensional measurements are projected into a lower-dimensional space so that dominant trends among samples can be visualized and interpreted." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-317)Overview of dimensional reduction. The schematic shows how high-dimensional measurements are projected into a lower-dimensional space so that dominant trends among samples can be visualized and interpreted.</p>
</div>

In the previous chapters, we looked at how to explore our data sets by visualizing many variables and manually identifying trends. Sometimes, we encounter data sets with so many variables, that it is not reasonable to manually select certain variables with which to create plots and manually search for trends. In these cases, we need dimensionality reduction - a set of techniques that helps us identify which variables are driving differences among our samples. In this course, we will conduct dimensionality reduction using `runMatrixAnalyses()`, a function that is loaded into your R Session when you run the source() command.

Matrix analyses can be a bit tricky to set up. There are two things that we can do to help us with this: (i) we will use a template for `runMatrixAnalyses()` (see below) and (ii) it is *critical* that we think about our data in terms of **samples** and **analytes**. Let's consider our Alaska lakes data set:


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

We can see that this dataset is comprised of measurements of various *analytes* (i.e. several chemical elements, as well as water_temp, and pH), in different *samples* (i.e. lakes). We need to tell the `runMatrixAnalyses()` function how each column relates to this samples and analytes structuree

## {-}

## pca {-}

"Which analytes are driving differences among my samples?"
"Which analytes in my data set are correlated?"

### theory {-}

PCA looks at all the variance in a high dimensional data set and chooses new axes within that data set that align with the directions containing highest variance. These new axes are called principal components. Let's look at an example:

<div class="figure" style="text-align: center">
<img src="https://thebustalab.github.io/integrated_bioanalytics/images/PCA.png" alt="Principal component rotation illustrated. The bold axes denote the new principal components that capture the largest variance directions, enabling us to describe complex data with fewer coordinates." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-319)Principal component rotation illustrated. The bold axes denote the new principal components that capture the largest variance directions, enabling us to describe complex data with fewer coordinates.</p>
</div>

In the example above, the three dimensional space can be reduced to a two dimensional space with the principal components analysis. New axes (principal components) are selected (bold arrows on left) that become the x and y axes in the principal components space (right).

We can run and visualize principal components analyses using the `runMatrixAnalyses()` function as in the example below. As you can see in the output, the command provides the sample_IDs, sample information, then the coordinates for each sample in the 2D projection (the "PCA plot") and the raw data, in case you wish to do further processing.


``` r

alaska_lake_data_wide <- pivot_wider(alaska_lake_data[,1:6], names_from = "element", values_from = "mg_per_L")
alaska_lake_data_wide
## # A tibble: 20 × 15
##    lake      park  water_temp    pH     C     N     P     Cl
##    <chr>     <chr>      <dbl> <dbl> <dbl> <dbl> <dbl>  <dbl>
##  1 Devil_Mo… BELA        6.46  7.69   3.4 0.028 0      10.4 
##  2 Imuruk_L… BELA       17.4   6.44   4.7 0.013 0       1.18
##  3 Kuzitrin… BELA        8.06  7.45   2   0     0       0.67
##  4 Lava_Lake BELA       20.2   7.42   8.3 0.017 0.001   2.53
##  5 North_Ki… BELA       11.3   8.04   4.3 0.037 0.001 337.  
##  6 White_Fi… BELA       12.0   7.82  12.3 0.034 0.006 105.  
##  7 Iniakuk_… GAAR        9.1   7.01   3.3 0.141 0       0.22
##  8 Kurupa_L… GAAR        9.3   7.03   2.1 0.043 0       0.13
##  9 Lake_Mat… GAAR       10.2   6.95   5.1 0     0       1.25
## 10 Lake_Sel… GAAR       15.1   7.15   4.2 0.107 0       0.11
## 11 Nutavukt… GAAR       17.6   6.88   4.5 0     0.001   0.18
## 12 Summit_L… GAAR       11.9   6.45   2.4 0     0.001   0.08
## 13 Takahula… GAAR        9.9   6.88   2.7 0.014 0       0.23
## 14 Walker_L… GAAR       15.3   7.22   1.3 0.19  0.001   0.19
## 15 Wild_Lake GAAR        5.5   6.98   6.5 0.13  0.001   0.31
## 16 Desperat… NOAT        2.95  6.34   2.1 0.005 0       0.2 
## 17 Feniak_L… NOAT        4.51  7.24   1.8 0     0       0.21
## 18 Lake_Kan… NOAT        5.36  6.56   8.5 0.005 0       0.55
## 19 Lake_Nar… NOAT       18.3   7.31   5.8 0     0       0.76
## 20 Okoklik_… NOAT        6.46  6.87   7.8 0     0       0.76
## # ℹ 7 more variables: S <dbl>, F <dbl>, Br <dbl>, Na <dbl>,
## #   K <dbl>, Ca <dbl>, Mg <dbl>

AK_lakes_pca <- runMatrixAnalyses(
  data = alaska_lake_data_wide,
  analysis = "pca",
  columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:15],
  columns_w_sample_ID_info = c("lake", "park")
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

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-321-1.png" alt="PCA scores for Alaskan lake chemistry. Points show each lake positioned by the first two principal components, with fill encoding the park and labels highlighting chemically distinct sites; distances capture multivariate differences across the analyte panel." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-321)PCA scores for Alaskan lake chemistry. Points show each lake positioned by the first two principal components, with fill encoding the park and labels highlighting chemically distinct sites; distances capture multivariate differences across the analyte panel.</p>
</div>

Great! In this plot we can see that White Fish Lake and North Killeak Lake, both in BELA park, are quite different from the other parks (they are separated from the others along dimension 1, i.e. the first principal component). At the same time, Wild Lake, Iniakuk Lake, Walker Lake, and several other lakes in GAAR park are different from all the others (they are separated from the others along dimension 2, i.e. the second principal component).

Important question: what makes the lakes listed above different from the others? Certainly some aspect of their chemistry, since that's the data that this analysis is built upon, but how do we determine which analyte(s) are driving the differences among the lakes that we see in the PCA plot?

### ordination plots {-}

Let's look at how to access the information about which analytes are major contributors to each principal component. This is important because it will tell you which analytes are associated with particular dimensions, and by extension, which analytes are associated with (and are markers for) particular groups in the PCA plot. This can be determined using an ordination plot. Let's look at an example. We can obtain the ordination plot information using `runMatrixAnalyses()` with `analysis = "pca_ord"`:


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
# AK_lakes_pca_ord <- runMatrixAnalysis(
#   data = alaska_lake_data,
#   analysis = c("pca_ord"),
#   column_w_names_of_multiple_analytes = "element",
#   column_w_values_for_multiple_analytes = "mg_per_L",
#   columns_w_values_for_single_analyte = c("water_temp", "pH"),
#   columns_w_additional_analyte_info = "element_type",
#   columns_w_sample_ID_info = c("lake", "park")
# )
# head(AK_lakes_pca_ord)

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

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-323-1.png" alt="Circular ordination plot for Alaskan lakes. Arrows mark analyte loadings scaled to the correlation circle, and labels flag the elements that dominate each principal axis so we can connect chemistry to lake groupings." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-323)Circular ordination plot for Alaskan lakes. Arrows mark analyte loadings scaled to the correlation circle, and labels flag the elements that dominate each principal axis so we can connect chemistry to lake groupings.</p>
</div>

Great! Here is how to read the ordination plot:

1. When considering one analyte's vector: the vector's projected value on an axis shows how much its variance is aligned with that principal component.

2. When considering two analyte vectors: the angle between two vectors indicates how correlated those two variables are. If they point in the same direction, they are highly correlated. If they meet each other at 90 degrees, they are not very correlated. If they meet at ~180 degrees, they are negatively correlated. If say that one analyte is "1.9" with respect to dimension 2 and another is "-1.9" with respect to dimension 2. Let's also say that these vectors are ~"0" with respect to dimension 1.

With the ordination plot above, we can now see that the abundances of K, Cl, Br, and Na are the major contributors of variance to the first principal component (or the first dimension). The abundances of these elements are what make White Fish Lake and North Killeak Lake different from the other lakes. We can also see that the abundances of N, S, and Ca are the major contributors to variance in the second dimension, which means that these elements ar what set Wild Lake, Iniakuk Lake, Walker Lake, and several other lakes in GAAR park apart from the rest of the lakes in the data set. It slightly easier to understand this if we look at an overlay of the two plots, which is often called a "biplot":


``` r
# AK_lakes_pca <- runMatrixAnalysis(
#   data = alaska_lake_data,
#   analysis = c("pca"),
#   column_w_names_of_multiple_analytes = "element",
#   column_w_values_for_multiple_analytes = "mg_per_L",
#   columns_w_values_for_single_analyte = c("water_temp", "pH"),
#   columns_w_additional_analyte_info = "element_type",
#   columns_w_sample_ID_info = c("lake", "park"),
#   scale_variance = TRUE
# )
# 
# AK_lakes_pca_ord <- runMatrixAnalysis(
#   data = alaska_lake_data,
#   analysis = c("pca_ord"),
#   column_w_names_of_multiple_analytes = "element",
#   column_w_values_for_multiple_analytes = "mg_per_L",
#   columns_w_values_for_single_analyte = c("water_temp", "pH"),
#   columns_w_additional_analyte_info = "element_type",
#   columns_w_sample_ID_info = c("lake", "park")
# )

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

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-324-1.png" alt="PCA biplot combining scores and loadings. Lakes are plotted as points coloured by park while analyte vectors overlay the same coordinate system, helping us link sample groupings to the drivers of chemical variance." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-324)PCA biplot combining scores and loadings. Lakes are plotted as points coloured by park while analyte vectors overlay the same coordinate system, helping us link sample groupings to the drivers of chemical variance.</p>
</div>

Note that you do not have to plot ordination data as a circular layout of segments. Sometimes it is much easier to plot (and interpret!) alternatives:


``` r
AK_lakes_pca_ord %>%
  ggplot(aes(x = Dim.1, y = analyte)) +
    geom_point(aes(fill = analyte), shape = 22, size = 3) +
    scale_fill_manual(values = discrete_palette) +
    theme_bw()
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-325-1.png" alt="Analyte loadings by principal component. The dot plot re-expresses the PCA loadings as coordinates along Dim.1, making it easy to compare how each element contributes relative to the others." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-325)Analyte loadings by principal component. The dot plot re-expresses the PCA loadings as coordinates along Dim.1, making it easy to compare how each element contributes relative to the others.</p>
</div>

### principal components {-}

We also can access information about the how much of the variance in the data set is explained by each principal component, and we can plot that using ggplot:


``` r
AK_lakes_pca_dim <- runMatrixAnalyses(
  data = alaska_lake_data_wide,
  analysis = c("pca_dim"),
  columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:15],
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

ggplot(
  data = AK_lakes_pca_dim, 
  aes(x = principal_component, y = percent_variance_explained)
) +
  geom_line() +
  geom_point() +
  theme_bw()
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-326-1.png" alt="Variance explained by principal components. The scree curve shows how much of the total chemical variability is captured by each component, informing how many dimensions to retain." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-326)Variance explained by principal components. The scree curve shows how much of the total chemical variability is captured by each component, informing how many dimensions to retain.</p>
</div>

Cool! We can see that the first principal component retains nearly 50% of the variance in the original dataset, while the second dimension contains only about 20%. We can derive an important notion about PCA visualization from this: the scales on the two axes need to be the same for distances between points in the x and y directions to be comparable. This can be accomplished using `coord_fixed()` as an addition to your ggplots.

### pcaVisualizer {-}

Static plots are great for reporting, but exploring PCA interactively can make it easier to understand the relationships between your samples and analytes. Our `source()` command provides an interactive app helper called `pcaVisualizer()` that wraps `runMatrixAnalyses()` and assembles a dashboard with coordinated plots.

<div class="figure" style="text-align: center">
<img src="https://thebustalab.github.io/integrated_bioanalytics/images/pca_visualizer.png" alt="Screenshot of the `pcaVisualizer()` dashboard showing the linked scores plot, loadings plot, and heatmap panels used to explore PCA interactively." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-327)Screenshot of the `pcaVisualizer()` dashboard showing the linked scores plot, loadings plot, and heatmap panels used to explore PCA interactively.</p>
</div>

The function takes three key arguments:

* `data`: a data frame that contains both the sample identifiers and the numeric analyte columns.
* `columns_w_sample_ID_info`: the column names that identify each sample (for example `c("lake", "park")`).
* `columns_w_values_for_single_analyte`: the numeric columns that should be analysed (for example `colnames(alaska_lake_data_wide)[3:15]`).

To launch the app for the Alaska lakes example:


``` r
pcaVisualizer(
  data = alaska_lake_data_wide,
  columns_w_sample_ID_info = c("lake", "park"),
  columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:15]
)
```

This opens a browser window (or RStudio viewer) with three linked panels:

* **PCA Plot** – the familiar scores plot. Use the sidebar dropdowns to map any sample ID column to colour or shape; the app reuses the `discrete_palette` defined earlier so colours stay consistent with the rest of this chapter.
* **Ordination Plot** – the loadings plot. Adjust the *Filter ordination plot* slider to hide vectors whose combined PC1/PC2 loading magnitude falls below your chosen threshold, making it easier to focus on the analytes that matter.
* **Analyte Abundance Heatmap** – a scaled (z-scored) heatmap for the analytes that pass the loading filter. You can order rows by Dim.1 or Dim.2 so the heatmap aligns with the direction of separation you care about in the scores plot.

Because `pcaVisualizer()` calls `runMatrixAnalyses()` under the hood, any transformations you perform on your data before launching the app (scaling, filtering, subsetting) carry through automatically. Use it as a quick sanity check while you are refining preprocessing steps, and once you are satisfied, you can reproduce the final view with scripted ggplot code for publication-quality figures.

## umap and tsne {-}

"How do non-linear dimensionality reduction techniques reveal hidden structure in my data?"
"Can I identify clusters or subtle gradients in my dataset that might be missed by PCA?"

UMAP (Uniform Manifold Approximation and Projection) is a non-linear dimensionality reduction technique that, unlike PCA, can capture both local and global data structure. It is especially useful when your data might have clusters or manifold structures that aren’t well-represented by linear combinations of features. You can think of umap as a technique that "unwinds" a dataset, rather than projecting it onto a plane. Compared with PCA, UMAP often reveals tighter clusters and curved trajectories because it optimizes neighborhood relationships rather than a straight-line variance objective. That flexibility is powerful, but it comes with trade-offs: UMAP axes have no direct numerical interpretation, results depend on hyperparameters such as `n_neighbors` and `min_dist`, and the algorithm is stochastic, so repeated runs can vary slightly. PCA, in contrast, is deterministic and easier to interpret because each principal component is a weighted combination of the original variables, but it may miss nonlinear structure that UMAP highlights.

We will work here with the `wine_quality` dataset. Our wine quality dataset doesn’t come with unique sample identifiers, which are essential for any matrix analysis. We can create these identifiers by adding a new column called ID to our data frame. The code below demonstrates how to do this:


``` r
wine_quality$ID <- seq(1, dim(wine_quality)[1], 1)
```

In this command, the dollar sign ($) is used to access (or create) the ID column within the wine_quality data frame. The seq() function generates a sequence of numbers starting at 1 and ending at the number of rows in the dataset (given by dim(wine_quality)[1]), with an increment of 1. This way, every sample is uniquely identified.

Before starting with non-linear techniques, it can be helpful to see how a linear method like PCA clusters our samples. We run the analysis using runMatrixAnalysis(), specifying the relevant columns that contain our continuous variables (columns 4 to 14), and passing along our sample identifier along with additional information such as wine type, quality_score, and quality_category.


``` r
wq <- runMatrixAnalyses(
  data = wine_quality,
  analysis = "pca",
  columns_w_values_for_single_analyte = colnames(wine_quality)[4:14],
  columns_w_sample_ID_info = c("ID", "type", "quality_score", "quality_category")
)
wq %>%
  arrange(quality_score) %>%
  ggplot(aes(x = Dim.1, y = Dim.2)) +
  geom_point(size = 3, aes(shape = type, fill = quality_score)) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_viridis() +
  theme_classic()
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-330-1.png" alt="PCA projection of wine chemistry. Samples are positioned by the first two components, with point shape distinguishing red and white wines and fill showing sensory quality scores; the layout highlights gradients that PCA captures." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-330)PCA projection of wine chemistry. Samples are positioned by the first two components, with point shape distinguishing red and white wines and fill showing sensory quality scores; the layout highlights gradients that PCA captures.</p>
</div>

In this PCA plot, each point represents a wine sample, with its position determined by the first two principal components. We’re using quality_score to fill the points with color, and different shapes to distinguish the wine type. This serves as a baseline for comparing how non-linear methods handle our data.

We can perform UMAP on the wine quality dataset just as easily as PCA. The code below shows how to run UMAP using runMatrixAnalysis():



``` r
runMatrixAnalyses(
  data = wine_quality,
  analysis = "umap",
  columns_w_values_for_single_analyte = colnames(wine_quality)[4:14],
  columns_w_sample_ID_info = c("ID", "type", "quality_score", "quality_category")
) %>%
  arrange(quality_score) %>%
  ggplot(aes(x = Dim_1, y = Dim_2)) +
  geom_point(size = 3, aes(shape = type, fill = quality_score)) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_viridis() +
  theme_classic()
```

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-331-1.png" alt="UMAP embedding of wine chemistry. The non-linear projection preserves neighbourhood relationships, revealing clusters driven by wine type and quality scores that complement the PCA view." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-331)UMAP embedding of wine chemistry. The non-linear projection preserves neighbourhood relationships, revealing clusters driven by wine type and quality scores that complement the PCA view.</p>
</div>

In the UMAP plot, each point’s coordinates (Dim_1 and Dim_2) are derived from UMAP’s algorithm, which strives to preserve the overall topology of the data. As a result, UMAP might reveal clusters or continuous gradients related to wine quality and type that aren’t as apparent with PCA.

t‑SNE (t‑Distributed Stochastic Neighbor Embedding) is another popular non-linear dimensionality reduction technique. It excels at revealing clusters but can sometimes distort the global structure in favor of preserving local relationships. Although we’re not showing a full t‑SNE example here, you can run t‑SNE using the runMatrixAnalysis() function by specifying analysis = "tsne" (see below). Note that tsne fails with samples that have duplicate analyte values, so we have to filter out any duplicates.


``` r
wine_quality_deduplicated <- wine_quality[!duplicated(wine_quality[4:14]),]
tsne_results <- runMatrixAnalyses(
  data = wine_quality_deduplicated,
  analysis = "tsne",
  columns_w_values_for_single_analyte = colnames(wine_quality)[4:14],
  columns_w_sample_ID_info = c("ID", "type", "quality_score", "quality_category")
)
```

From there, you could plot the t‑SNE dimensions in a similar fashion to the PCA and UMAP examples.

Both UMAP and t‑SNE provide powerful alternatives to PCA when your data’s structure is non-linear. They can help uncover hidden patterns and clusters by focusing on preserving local relationships—UMAP while maintaining a sense of global structure, and t‑SNE by emphasizing the neighborhood structure of the data.

## {-}

## further reading {-}

- PCA Explanation Video: This YouTube video provides a detailed and visually intuitive explanation of Principal Component Analysis (PCA), breaking down complex concepts with clear examples and graphics. It is part of a curated playlist that covers a variety of topics related to data visualization and statistical analysis. [Watch the PCA video](https://www.youtube.com/watch?v=FgakZw6K1QQ&list=PLblh5JKOoLUIcdlgu78MnlATeyx4cEVeR).

- Understanding UMAP: This blog post from Pair Code delves into the fundamentals of UMAP, explaining both the intuition behind the algorithm and its practical applications in data analysis. It provides an accessible overview that bridges the gap between theoretical concepts and real-world use cases, making it a valuable read for both beginners and advanced users. [Explore the Understanding UMAP Article](https://pair-code.github.io/understanding-umap/).

- UMAP: Mathematical Details (clearly explained!!!) This YouTube video offers a detailed explanation of the mathematical underpinnings of UMAP, breaking down the algorithm in a clear and approachable manner. It is an excellent resource for viewers who want to deepen their understanding of how UMAP works behind the scenes.
[Watch the UMAP Mathematical Details Video](https://www.youtube.com/watch?v=jth4kEvJ3P8).


# flat clustering {-}

"Do my samples fall into definable clusters?"

## {-}

## kmeans {-}

One of the questions we've been asking is "which of my samples are most closely related?". We've been answering that question using clustering. However, now that we know how to run principal components analyses, we can use another approach. This alternative approach is called k-means, and can help us decide how to assign our data into clusters. It is generally desirable to have a small number of clusters, however, this must be balanced by not having the variance within each cluster be too big. To strike this balance point, the elbow method is used. For it, we must first determine the maximum within-group variance at each possible number of clusters. An illustration of this is shown in **A** below:

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/kmeans.png" width="100%" style="display: block; margin: auto;" />

One we know within-group variances, we find the "elbow" point - the point with minimum angle theta - thus picking the outcome with a good balance of cluster number and within-cluster variance (illustrated above in **B** and **C**.)

Let's try k-means using `runMatrixAnalysis`. For this example, let's run it on the PCA projection of the alaska lakes data set. We can set `analysis = "kmeans"`. When we do this, an application will load that will show us the threshold value for the number of clusters we want. We set the number of clusters and then close the app. In the context of markdown document, simply provide the number of clusters to the `parameters` argument:



``` r
alaska_lake_data %>%
  select(-element_type) %>%
  pivot_wider(names_from = "element", values_from = "mg_per_L") -> alaska_lake_data_wide

alaska_lake_data_pca <- runMatrixAnalyses(
    data = alaska_lake_data_wide,
    analysis = c("pca"),
    columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:dim(alaska_lake_data_wide)[2]],
    columns_w_sample_ID_info = c("lake", "park")
)

alaska_lake_data_pca_clusters <- runMatrixAnalyses(
    data = alaska_lake_data_pca,
    analysis = c("kmeans"),
    parameters = c(5),
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

<img src="index_files/figure-html/unnamed-chunk-353-1.png" width="100%" style="display: block; margin: auto;" />

## dbscan {-}

There is another method to define clusters that we call dbscan. In this method, not all points are necessarily assigned to a cluster, and we define clusters according to a set of parameters, instead of simply defining the number of clusteres, as in kmeans. In interactive mode, `runMatrixAnalysis()` will again load an interactive means of selecting parameters for defining dbscan clusters ("k", and "threshold"). In the context of markdown document, simply provide "k" and "threshold" to the `parameters` argument:


``` r
alaska_lake_data_pca <- runMatrixAnalyses(
    data = alaska_lake_data_wide,
    analysis = c("pca"),
    columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:dim(alaska_lake_data_wide)[2]],
    columns_w_sample_ID_info = c("lake", "park")
) 

alaska_lake_data_pca_clusters <- runMatrixAnalyses(
    data = alaska_lake_data_pca,
    analysis = c("dbscan"),
    parameters = c(4, 0.45),
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

<img src="index_files/figure-html/unnamed-chunk-355-1.png" width="100%" style="display: block; margin: auto;" />

## summarize by cluster {-}

One more important point: when using kmeans or dbscan, we can use the clusters as groupings for summary statistics. For example, suppose we want to see the differences in abundances of certain chemicals among the clusters:


``` r
alaska_lake_data_pca <- runMatrixAnalysis(
  data = alaska_lake_data_wide,
  analysis = c("pca"),
  columns_w_values_for_single_analyte = colnames(alaska_lake_data_wide)[3:dim(alaska_lake_data_wide)[2]],
  columns_w_sample_ID_info = c("lake", "park")
)

alaska_lake_data_pca_clusters <- runMatrixAnalyses(
  data = alaska_lake_data_pca,
  analysis = c("dbscan"),
  parameters = c(4, 0.45),
  columns_w_values_for_single_analyte = c("Dim.1", "Dim.2"),
  columns_w_sample_ID_info = "sample_unique_ID"
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

<img src="index_files/figure-html/unnamed-chunk-356-1.png" width="100%" style="display: block; margin: auto;" />
 
## {-}

## further reading {-}

- STHDA DBSCAN tutorial: STHDA walks through the DBSCAN algorithm, explains how `eps` and `MinPts` control density-based clusters, and demonstrates complete R examples using the `dbscan` and `fpc` packages. [Read the STHDA DBSCAN guide](http://www.sthda.com/english/wiki/wiki.php?id_contents=7940).

- Ryan Wingate on hierarchical vs. density-based clustering: Ryan Wingate compares agglomerative methods with DBSCAN, showing how linkage choices and parameter tuning change the resulting partitions and offering intuition for when to use each approach. [Review the walkthrough](https://ryanwingate.com/intro-to-machine-learning/unsupervised/hierarchical-and-density-based-clustering/).

- Ryan Wingate dendrogram reference: This figure distills the merge sequence from the accompanying tutorial into a color-coded dendrogram, making it easy to point out where to cut the tree when discussing flat clusters. [Open the dendrogram figure](https://ryanwingate.com/intro-to-machine-learning/unsupervised/hierarchical-and-density-based-clustering/hierarchical-4.png).

- GeeksforGeeks DBSCAN in R: GeeksforGeeks uses a step-by-step R workflow to run DBSCAN, visualize clusters, and explain how adjusting `eps` and `MinPts` affects noise handling. [Follow the tutorial](https://www.geeksforgeeks.org/dbscan-clustering-in-r-programming/).


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

ggplot(aquifers_summarized) + geom_col(aes(x = n_wells, y = aquifer_code))
```

<img src="index_files/figure-html/unnamed-chunk-383-1.png" width="100%" style="display: block; margin: auto;" />

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

ggplot(K_data_1_6, aes(x = aquifer_code, y = abundance)) +
    geom_boxplot() +
    geom_point()
```

<img src="index_files/figure-html/unnamed-chunk-386-1.png" width="100%" style="display: block; margin: auto;" />

Are these data normally distributed? Do they have similar variance? Let's get a first approximation by looking at a plot:


``` r
K_data_1_6 %>%
  ggplot(aes(x = abundance)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~aquifer_code) +
    geom_density(aes(y = ..density..*10), color = "blue")
```

<img src="index_files/figure-html/unnamed-chunk-387-1.png" width="100%" style="display: block; margin: auto;" />

Based on this graphic, it's hard to say! Let's use a statistical test to help. When we want to run the Shaprio test, we are looking to see if each group has normally distributed here (here group is "aquifer_code", i.e. aquifer_1 and aquifer_6). This means we need to `group_by(aquifer_code)` before we run the test:


``` r
K_data_1_6 %>%
  group_by(aquifer_code) %>% 
  shapiroTest(abundance)
## # A tibble: 2 × 4
##   aquifer_code variable statistic     p
##   <chr>        <chr>        <dbl> <dbl>
## 1 aquifer_1    values       0.885 0.102
## 2 aquifer_6    values       0.914 0.239
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

A p-value of 0.012! This is below 0.05, meaning that if the two aquifers truly had the same mean, we would only observe a difference at least this extreme about 1.2% of the time. That gives us fairly strong evidence that the means differ. Suppose that our data had not passed the Shapiro and/or Levene tests. We would then need to use a Wilcox test. The Wilcox test is a non-parametric test, which means that it does not use a normal distribution to model the data in order to make comparisons. This means that is a less powerful test than the t-test, which means that it is less likely to detect a difference in the means, assuming there is one. For fun, let's try that one out and compare the p-values from the two methods:


``` r
K_data_1_6 %>%
  wilcoxTest(abundance ~ aquifer_code)
## # A tibble: 1 × 7
##   .y.       group1    group2       n1    n2 statistic      p
## * <chr>     <chr>     <chr>     <int> <int>     <dbl>  <dbl>
## 1 abundance aquifer_1 aquifer_6    12    12      33.5 0.0282
```

A p-value of 0.028! This is higher than the value given by the t-test (0.012). That is because the Wilcox test is a less powerful test: it is less likely to detect differences in means, assuming they exist. Interpreting the number in the same way, if the two groups really had identical distributions, we would get results this extreme only about 2.8% of the time.

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

ggplot(data = K_data, aes(y = aquifer_code, x = abundance)) +
  geom_boxplot() +
  geom_point(color = "maroon", alpha = 0.6, size = 3)
```

<img src="index_files/figure-html/unnamed-chunk-392-1.png" width="100%" style="display: block; margin: auto;" />

Let's check visually to see if each group is normally distributed and to see if they have roughly equal variance:


``` r
K_data %>%
  group_by(aquifer_code) %>%
  ggplot(aes(x = abundance)) + 
    geom_histogram(bins = 30) +
    facet_wrap(~aquifer_code) +
    geom_density(aes(y = ..density..*10), colour = "blue")
```

<img src="index_files/figure-html/unnamed-chunk-393-1.png" width="100%" style="display: block; margin: auto;" />

Again, it is somewhat hard to tell visually if these data are normally distributed. It seems pretty likely that they have different variances about the means, but let's check using the Shapiro and Levene tests. Don't forget: with the Shaprio test, we are looking within each group and so need to `group_by()`, with the Levene test, we are looking across groups, and so need to provide a `y~x` formula:


``` r
K_data %>%
  group_by(aquifer_code) %>% 
  shapiroTest(abundance)
## # A tibble: 10 × 4
##    aquifer_code variable statistic         p
##    <chr>        <chr>        <dbl>     <dbl>
##  1 aquifer_1    values       0.885 0.102    
##  2 aquifer_10   values       0.864 0.163    
##  3 aquifer_2    values       0.913 0.459    
##  4 aquifer_3    values       0.893 0.363    
##  5 aquifer_4    values       0.948 0.421    
##  6 aquifer_5    values       0.993 0.972    
##  7 aquifer_6    values       0.914 0.239    
##  8 aquifer_7    values       0.915 0.355    
##  9 aquifer_8    values       0.842 0.220    
## 10 aquifer_9    values       0.790 0.0000214
```


``` r
K_data %>%
  leveneTest(abundance ~ aquifer_code)
## # A tibble: 1 × 4
##     df1   df2 statistic       p
##   <int> <int>     <dbl>   <dbl>
## 1     9    96      2.95 0.00387
```

Based on these tests, it looks like the data for aquifer 9 is significantly different from a normal distribution (Shaprio test p < 0.05), and the variances are certainly different from one another (Levene test p < 0.05).

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

A pretty small p-value! Under the null hypothesis that all aquifers have the same mean, observing differences at least this extreme would be very unlikely, so we have evidence that at least one mean differs. But, WHICH are different from one another though? For this, we need to run Tukey's Honest Significant Difference test (implemented using `tukey_hsd`). This will essentially run t-test on all the pairs of study subjects that we can derive from our data set (in this example, aquifer_1 vs. aquifer_2, aquifer_1 vs. aquifer_3, etc.). After that, it will correct the p-values according to the number of comparisons that it performed. This controls the rate of type I error that we can expect from the test. These corrected values are provided to us in the `p.adj` column.


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

<img src="index_files/figure-html/unnamed-chunk-399-1.png" width="100%" style="display: block; margin: auto;" />

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

A pretty small p-value! This is higher than the p-value from running ANOVA on the same data (remember, the Kruskal test is less powerful). Never the less, the value is still well below 0.05, meaning that data this extreme would be rare if all group medians were identical. So, how do we determine WHICH are different from one another? When we ran ANOVA the follow-up test (the post hoc test) was Tukey's HSD. After the Kruskal test, the post hoc test we use is the Dunn test. Let's try:


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

ggplot(data = K_data, aes(y = aquifer_code, x = abundance)) +
  geom_boxplot() +
  geom_point(color = "black", alpha = 0.4, size = 2) +
  scale_x_continuous(name = "Potassium abundance", breaks = seq(0,10,1)) +
  scale_y_discrete(name = "Aquifer code") +
  geom_text(data = groups_based_on_dunn, aes(y = treatment, x = 9, label = group)) +
  theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-402-1.png" width="100%" style="display: block; margin: auto;" />

Note that these groupings are different from those generated by ANOVA/Tukey.

## pairs of means {-}

Oftentimes we have more than two means to compare, but rather than wanting to compare all means at once, we want to compare them in a pairwise fashion. For example, suppose we want to know if any of the aquifers contain different amounts of Na and Cl. We are not interested in testing for differences among *all* values of Na and Cl, rather, we want to test all *pairs* of Na and Cl values arising from each aquifer. That is to say, we want to compare the means in each facet of the plot below:


``` r
hawaii_aquifers %>%
  filter(analyte %in% c("Na", "Cl")) %>%
  ggplot(aes(x = analyte, y = abundance)) + geom_violin() + geom_point() + facet_grid(.~aquifer_code)
```

<img src="index_files/figure-html/unnamed-chunk-403-1.png" width="100%" style="display: block; margin: auto;" />

Fortunately, we can use an approach that is very similar to the what we've learned in the earlier portions of this chapter, just with minor modifications. Let's have a look! We start with the Shapiro and Levene tests, as usual (note that we group using two variables when using the Shapiro test so that each analyte within each aquifer is considered as an individual distribution):


``` r
hawaii_aquifers %>%
  filter(analyte %in% c("Na", "Cl")) %>%
  group_by(analyte, aquifer_code) %>%
  shapiroTest(abundance)
## Note: Found a group with exactly 3 data points where at least two were identical.
## A tiny bit of noise was added to these data points, because the Shapiro–Wilk test
## p-value can be artificially low in such edge cases due to test mechanics.
## Consider carefully whether you want to interpret the results of the test, even with noise added.
## # A tibble: 20 × 5
##    analyte aquifer_code variable statistic        p
##    <chr>   <chr>        <chr>        <dbl>    <dbl>
##  1 Cl      aquifer_1    values       0.900 1.59e- 1
##  2 Cl      aquifer_10   values       0.486 1.09e- 5
##  3 Cl      aquifer_2    values       0.869 2.24e- 1
##  4 Cl      aquifer_3    values       0.750 5.19e- 6
##  5 Cl      aquifer_4    values       0.903 7.49e- 2
##  6 Cl      aquifer_5    values       0.849 2.24e- 1
##  7 Cl      aquifer_6    values       0.741 2.15e- 3
##  8 Cl      aquifer_7    values       0.893 2.12e- 1
##  9 Cl      aquifer_8    values       0.878 3.17e- 1
## 10 Cl      aquifer_9    values       0.420 2.68e-10
## 11 Na      aquifer_1    values       0.886 1.06e- 1
## 12 Na      aquifer_10   values       0.593 2.26e- 4
## 13 Na      aquifer_2    values       0.884 2.88e- 1
## 14 Na      aquifer_3    values       0.822 1.69e- 1
## 15 Na      aquifer_4    values       0.933 2.41e- 1
## 16 Na      aquifer_5    values       0.827 1.61e- 1
## 17 Na      aquifer_6    values       0.764 3.80e- 3
## 18 Na      aquifer_7    values       0.915 3.51e- 1
## 19 Na      aquifer_8    values       0.855 2.53e- 1
## 20 Na      aquifer_9    values       0.531 3.97e- 9
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

<img src="index_files/figure-html/unnamed-chunk-408-1.png" width="100%" style="display: block; margin: auto;" />

## {-}

## further reading {-}

- Datanovia on comparing multiple means: This step-by-step tutorial shows how to run one-way and repeated-measures ANOVA in R, diagnose assumptions, and follow up with Tukey or pairwise comparisons using tidyverse-friendly code. [Open the Datanovia guide](https://www.datanovia.com/en/courses/comparing-multiple-means-in-r/).

- Statistics by Jim on parametric vs nonparametric tests: Jim Frost breaks down when parametric tests are appropriate, the trade-offs of switching to rank-based alternatives, and how violation of assumptions affects power and interpretation. [Read the comparison](https://statisticsbyjim.com/hypothesis-testing/nonparametric-parametric-tests/).

- Ulrich Dirnagl on the p value wars: This short commentary captures the ongoing debates about redefining statistical significance, highlights the historical context of *p* < 0.05, and weighs the risks of rigid thresholds versus broader inferential thinking. [Consider the perspective](https://link.springer.com/article/10.1007/s00259-019-04467-5).

- eLife forum on common statistical mistakes: The authors catalog frequent errors such as double dipping, underpowered designs, and misreported effect sizes, and pair each with pragmatic fixes for both authors and reviewers. [Review the checklist](https://elifesciences.org/articles/48175).

- Wolfgang Huber on tiny *p*-values: Huber explains why extremely small *p*-values emerge in high-throughput experiments, how to think about multiplicity and effect sizes, and what to report alongside adjusted significance. [Study the reporting guidance](https://www.cell.com/cell-systems/fulltext/S2405-4712(19)30071-7?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2405471219300717%3Fshowall%3Dtrue).


<!-- ## exercises {-}

Using the `hawaii_aquifers` data set or the `tequila_chemistry` data set, please complete the following:

1. Choose one analyte and filter the data so only the rows for that analyte are shown.

2. Choose two of the aquifers (or bottles of tequila). Are the mean abundances for your chosen analyte different in these two aquifers? Don't forget to test your data for normality and homogeneity of variance before selecting a statistical test. Use a plot to illustrate whether the means are similar or different. If you wish, you can use the figure caption to help illustrate whether means are significantly different or not. Either way, be sure to describe the tests you ran in your figure caption, and how you interpreted the results of the tests. For example, you should include a statment like "Groups sharing a letter label (a, b) are not significantly different (ANOVA/Tukey HSD tests, p > 0.05). Groups with different letters indicate a significant difference (p < 0.05).", or, if you are using asterisks then perhaps something like "Significant differences between groups are indicated by asterisks: \* = p < 0.05, \*\* = p < 0.01 (Kruskal / Dunn tests)."

3. Choose a second analyte, different from the first one you chose. Considering all the aquifers (or all bottles of tequila) in the dataset, do any of them have the same abundance of this analyte? Again, don't forget about normality and homogeneity of variance tests. Use a plot to illustrate your answer. Be sure to describe the tests you ran in your figure caption, and how you interpreted the results of the tests.

4. Repeat #3 above, but switch the type of test used (i.e. use non-parametric if you used parametric for #3 and vice-versa). Compare the *p* values and *p* groups obtained by the two methods. Use a graphic to illustrate this. Why are they different? -->

<!-- end -->

________________________________________________________________________________________________
________________________________________________________________________________________________
________________________________________________________________________________________________

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
metabolomics_data[] <- lapply(metabolomics_data, function(x) Hmisc::impute(x, median(x, na.rm = TRUE)))
any(is.na(metabolomics_data))
## [1] FALSE
```

## single linear regression {-}

We will start with some of the simplest models - linear models. There are a variety of ways to build linear models in R, but we will use `buildModel`, as mentioned above. First, we will use least squares regression to model the relationship between input and output variables. Suppose we want to know if the abundances of iso-Leucine and Valine are related in our metabolomics dataset:


``` r
ggplot(metabolomics_data) +
  geom_point(aes(x = `iso-Leucine`, y = Valine))
```

<img src="index_files/figure-html/unnamed-chunk-441-1.png" width="100%" style="display: block; margin: auto;" />

It looks like there might be a relationship! Let's build an linear regression model and use it inferentially to examine the details of that that relationship:


``` r
basic_regression_model <- buildModel2(
  data = metabolomics_data,
  model_type = "linear_regression",
  input_variables = "iso-Leucine",
  output_variable = "Valine"
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

It shows us the r-squared, the total and residual sum of squares, the intercept (b in y = mx + b), and the coefficient for iso-Leucine (i.e. the slope, m), as well some other things (we will talk about them in a second).

We can also use a function called `predictWithModel` to make some predictions using the model. Let's try that for iso-Leucine and Valine. What we do is give it the model, and then tell it what values we want to predict for. In this case, we want to predict the abundance of Valine for each value of iso-Leucine in our data set. We can do that like this:


``` r
predicted_Valine_values <- predictWithModel(
  data = metabolomics_data,
  model_type = "linear_regression",
  model = basic_regression_model$model
)
head(predicted_Valine_values)
## # A tibble: 6 × 1
##   value
##   <dbl>
## 1  4.61
## 2  4.80
## 3  4.96
## 4  4.78
## 5  4.66
## 6  4.69
```

So, `predictWithModel` is using the model to predict Valine values from iso-Leucine. However, note that we have the measured Valine values in our data set. We can compare the predicted values to the measured values to see how well our model is doing. We can do that like this:


``` r
predictions_from_basic_linear_model <- data.frame(
    iso_Leucine_values = metabolomics_data$`iso-Leucine`,
    predicted_Valine_values = predicted_Valine_values,
    measured_Valine_values = metabolomics_data$Valine
)
names(predictions_from_basic_linear_model)[2] <- "predicted_Valine_values"

plot1 <- ggplot() +
    geom_line(
        data = predictions_from_basic_linear_model,
        aes(x = iso_Leucine_values, y = predicted_Valine_values), color = "red"
    ) +
    geom_point(
        data = predictions_from_basic_linear_model,
        aes(x = iso_Leucine_values, y = predicted_Valine_values), color = "red"
    ) +
    geom_point(
        data = metabolomics_data,
        aes(x = `iso-Leucine`, y = Valine), color = "blue"
    )
plot1
```

<img src="index_files/figure-html/unnamed-chunk-446-1.png" width="100%" style="display: block; margin: auto;" />

Very good. Now let's talk about evaluating the quality of our model. For this we need some means of assessing how well our line fits our data. We will use residuals - the distance between each of our points and our line.


``` r
ggplot(predictions_from_basic_linear_model) +
  geom_point(aes(x = iso_Leucine_values, y = measured_Valine_values)) +
  geom_line(aes(x = iso_Leucine_values, y = predicted_Valine_values)) +
  geom_segment(aes(x = iso_Leucine_values, y = measured_Valine_values, xend = iso_Leucine_values, yend = predicted_Valine_values))
```

<img src="index_files/figure-html/unnamed-chunk-447-1.png" width="100%" style="display: block; margin: auto;" />

We can calculate the sum of the squared residuals:


``` r
sum(
  (predictions_from_basic_linear_model$measured_Valine_values - predictions_from_basic_linear_model$predicted_Valine_values)^2
, na.rm = TRUE)
## [1] 1.285459
```

Cool! Let's call that the "residual sum of the squares". So... does that mean our model is good? I don't know. We have to compare that number to something. Let's compare it to a super simple model that is just defined by the mean y value of the input data.


``` r
ggplot(metabolomics_data) +
  geom_point(aes(x = `iso-Leucine`, y = Valine)) +
  geom_hline(aes(yintercept = mean(Valine, na.rm = TRUE)))
```

<img src="index_files/figure-html/unnamed-chunk-449-1.png" width="100%" style="display: block; margin: auto;" />

A pretty bad model, I agree. How much better is our linear model that the flat line model? Let's create a measure of the distance between each point and the point predicted for that same x value on the model:


``` r
ggplot(metabolomics_data) +
  geom_point(aes(x = `iso-Leucine`, y = Valine)) +
  geom_hline(aes(yintercept = mean(Valine, na.rm = TRUE))) +
  geom_segment(aes(x = `iso-Leucine`, y = Valine, xend = `iso-Leucine`, yend = mean(Valine, na.rm = TRUE)))
```

<img src="index_files/figure-html/unnamed-chunk-450-1.png" width="100%" style="display: block; margin: auto;" />

``` r

sum(
  (metabolomics_data$Valine - mean(metabolomics_data$Valine, na.rm = TRUE))^2
, na.rm = TRUE)
## [1] 2.154464
```

Cool. Let's call that the "total sum of the squares", and now we can compare that to our "residual sum of the squares": 


``` r
residual_sum_of_squares <- sum(
  (predictions_from_basic_linear_model$measured_Valine_values - predictions_from_basic_linear_model$predicted_Valine_values)^2,
  na.rm = TRUE
)
total_sum_of_squares <- sum(
  (metabolomics_data$Valine - mean(metabolomics_data$Valine, na.rm = TRUE))^2,
  na.rm = TRUE
)
1 - (residual_sum_of_squares / total_sum_of_squares)
## [1] 0.403351
```

Alright. That is our R squared value. It is equal to 1 minus the ratio of the "residual sum of the squares" to the "total sum of the squares". You can think of the R squared value as:
- The amount of variance in the response explained by the dependent variable.
- How much better the line of best fit describes the data than the flat line.
Now, let's put it all together and make it pretty:


``` r
top <- ggplot() +
    geom_line(
        data = predictions_from_basic_linear_model,
        aes(x = iso_Leucine_values, y = predicted_Valine_values), color = "red"
    ) +
    geom_point(
        data = predictions_from_basic_linear_model,
        aes(x = iso_Leucine_values, y = predicted_Valine_values), color = "red"
    ) +
    geom_point(
        data = metabolomics_data,
        aes(x = `iso-Leucine`, y = Valine), color = "blue"
    ) +
    annotate(geom = "table",
      x = 3.25,
      y = 5,
      label = list(select(basic_regression_model$metrics, variable, value))
    ) +
    coord_cartesian(ylim = c(4.25,5.25)) +
    theme_bw()

bottom <- ggplot(predictions_from_basic_linear_model) +
  geom_col(
    aes(x = iso_Leucine_values, y = measured_Valine_values - predicted_Valine_values),
    width = 0.03, color = "black", position = "dodge", alpha = 0.5
  ) +
  theme_bw()

cowplot::plot_grid(top, bottom, ncol = 1, labels = "AUTO", rel_heights = c(2,1))
```

<img src="index_files/figure-html/unnamed-chunk-452-1.png" width="100%" style="display: block; margin: auto;" />

## multiple linear regression {-}

Cool! Now let's try a multiple linear regression model. This is the same as a simple linear regression model, but with more than one predictor variable. Simple and multiple linear regression are both statistical methods used to explore the relationship between one or more independent variables (predictor variables) and a dependent variable (outcome variable). Simple linear regression involves one independent variable to predict the value of one dependent variable, utilizing a linear equation of the form y = mx + b. Multiple linear regression extends this concept to include two or more independent variables, with a typical form of  y = m1x1 + m2x2 + ... + b, allowing for a more complex representation of relationships among variables. While simple linear regression provides a straight-line relationship between the independent and dependent variables, multiple linear regression can model a multi-dimensional plane in the variable space, providing a more nuanced understanding of how the independent variables collectively influence the dependent variable. The complexity of multiple linear regression can offer more accurate predictions and insights, especially in scenarios where variables interact or are interdependent, although it also requires a more careful consideration of assumptions and potential multicollinearity among the independent variables. Let's try it with the first 30 metabolites in our data set:


``` r

single_input_regression_model <- buildModel2(
  data = metabolomics_data,
  model_type = "linear_regression",
  input_variables = "iso-Leucine",
  output_variable = "Valine"
)

single_input_regression_model$metrics %>%
  filter(type == "coefficient") %>%
  arrange(desc(abs(value)))
##        variable  value std_err        type p_value
## 1   (Intercept) 2.7843  0.2459 coefficient       0
## 2 `iso-Leucine` 0.5080  0.0648 coefficient       0
##   p_value_adj
## 1           0
## 2           0

multiple_input_regression_model <- buildModel2(
  data = metabolomics_data,
  model_type = "linear_regression",
  input_variables = colnames(metabolomics_data)[3:32],
  output_variable = "Valine"
)

multiple_input_regression_model$metrics %>% filter(type == "statistic")
##               variable  value std_err      type p_value
## 1            r_squared 0.8942      NA statistic      NA
## 2    total_sum_squares 2.1545      NA statistic      NA
## 3 residual_sum_squares 0.2280      NA statistic      NA
##   p_value_adj
## 1          NA
## 2          NA
## 3          NA

multiple_input_regression_model$metrics %>%
  filter(type == "coefficient") %>%
  arrange(desc(abs(value)))
##                        variable   value std_err        type
## 1                   (Intercept) -8.4392  2.1225 coefficient
## 2         `5-Aminovaleric Acid`  0.8629  0.0584 coefficient
## 3                  Homocysteine -0.1139  0.0558 coefficient
## 4     `Alpha-Ketoglutaric Acid` -0.0882  0.0780 coefficient
## 5           `Pyroglutamic Acid`  0.0750  0.0507 coefficient
## 6                    Cadaverine -0.0636  0.0531 coefficient
## 7                     Glycerate -0.0594  0.0544 coefficient
## 8                     Sarcosine  0.0547  0.0276 coefficient
## 9                     Carnitine  0.0521  0.0433 coefficient
## 10        `Glucose 1-phosphate`  0.0490  0.0361 coefficient
## 11       `Phosphoglyceric Acid` -0.0465  0.0227 coefficient
## 12              MethylSuccinate  0.0394  0.0508 coefficient
## 13                     Tyramine  0.0264  0.0556 coefficient
## 14                     Pyruvate  0.0183  0.0302 coefficient
## 15                   Pipecolate  0.0175  0.0507 coefficient
## 16   `2-Hydroxyisovaleric Acid` -0.0145  0.0198 coefficient
## 17                      Betaine -0.0144  0.0244 coefficient
## 18                      Taurine -0.0144  0.0394 coefficient
## 19            `isoValeric Acid` -0.0141  0.0789 coefficient
## 20          `1-Methylhistidine`  0.0137  0.0091 coefficient
## 21                   Cysteamine  0.0116  0.0123 coefficient
## 22     `2-Aminoisobutyric acid` -0.0106  0.0400 coefficient
## 23             Guanidinoacetate  0.0105  0.0285 coefficient
## 24               `Malonic Acid`  0.0093  0.0137 coefficient
## 25                     Creatine  0.0087  0.0173 coefficient
## 26            `N-AcetylGlycine` -0.0087  0.0161 coefficient
## 27          `1-Methylhistamine` -0.0034  0.0140 coefficient
## 28 `3-Methyl-2-Oxovaleric Acid`  0.0030  0.0152 coefficient
## 29               `Fumaric Acid` -0.0026  0.0129 coefficient
## 30                   Creatinine -0.0024  0.0381 coefficient
## 31           `4-Hydroxyproline`  0.0001  0.0173 coefficient
##       p_value p_value_adj
## 1  0.00018570    0.005571
## 2  0.00000000    0.000000
## 3  0.04539385    1.000000
## 4  0.26252957    1.000000
## 5  0.14424229    1.000000
## 6  0.23545690    1.000000
## 7  0.27908387    1.000000
## 8  0.05203798    1.000000
## 9  0.23333675    1.000000
## 10 0.17966959    1.000000
## 11 0.04448459    1.000000
## 12 0.44116148    1.000000
## 13 0.63648128    1.000000
## 14 0.54731576    1.000000
## 15 0.73164052    1.000000
## 16 0.46726689    1.000000
## 17 0.55755863    1.000000
## 18 0.71615475    1.000000
## 19 0.85827398    1.000000
## 20 0.13825735    1.000000
## 21 0.34897953    1.000000
## 22 0.79070612    1.000000
## 23 0.71355085    1.000000
## 24 0.50039850    1.000000
## 25 0.61601809    1.000000
## 26 0.59213968    1.000000
## 27 0.81043509    1.000000
## 28 0.84546515    1.000000
## 29 0.84339908    1.000000
## 30 0.95064942    1.000000
## 31 0.99539331    1.000000
```

Even though the multiple-input model predicts Valine more accurately overall (based on r-squared values), it is helpful to quantify which metabolites drive that improvement. The regression summary (`multiple_input_regression_model$metrics`) reports the coefficient estimate, standard error, and p-value for each predictor. To decide which variables are the most important contributors:

- Coefficient value: magnitude and sign describe the expected change in the output when that predictor moves by one unit while all other variables are held constant. However, please note that because our predictors are on different scales, centering and scaling them (e.g., with `metabolomics_data %>% mutate(across(3:32, scale))`) lets us compare the absolute size of coefficients directly.
- `p_value` (or an adjusted p-value when many metabolites are examined at once) is the p-value obtained by testing whether the predictors coefficient is significantly different from zero. Small p-values highlight predictors that are likely to matter after accounting for the rest of the variables.

In this Valine model, `5-Aminovaleric Acid` has the largest absolute coefficient, so the third panel below plots its raw relationship with Valine alongside the two sets of predictions to highlight why it stands out in the summary table.



``` r
model_comparison_data <- data.frame(
  measured_Valine_values = metabolomics_data$Valine,
  single_input_predicted_Valine_values = predictWithModel(
    data = metabolomics_data,
    model_type = "linear_regression",
    model = basic_regression_model$model
  )$value,
  multiple_input_predicted_Valine_values = predictWithModel(
    data = metabolomics_data,
    model_type = "linear_regression",
    model = multiple_input_regression_model$model
  )$value,
  measured_5AA_values = metabolomics_data$`5-Aminovaleric Acid`
)
  
plot1 <- ggplot(model_comparison_data) + geom_point(aes(
  x = measured_Valine_values, y = single_input_predicted_Valine_values
))

plot2 <- ggplot(model_comparison_data) + geom_point(aes(
  x = measured_Valine_values, y = multiple_input_predicted_Valine_values
))

plot3 <- ggplot(model_comparison_data) + geom_point(aes(
  x = measured_Valine_values, y = measured_5AA_values
))

plot_grid(plot1, plot2, plot3, nrow = 1)
```

<img src="index_files/figure-html/unnamed-chunk-454-1.png" width="100%" style="display: block; margin: auto;" />





<!-- # ggplot() + -->
<!-- #   geom_point( -->
<!-- #     data = metabolomics_data, -->
<!-- #     aes(x = `iso-Leucine`, y = Valine), fill = "gold", shape = 21, color = "black" -->
<!-- #   ) + -->
<!-- #   geom_line(aes( -->
<!-- #     x = metabolomics_data$`iso-Leucine`, -->
<!-- #     y = mean(metabolomics_data$Valine) -->
<!-- #   ), color = "grey") + -->
<!-- #   geom_line(aes( -->
<!-- #     x = metabolomics_data$`iso-Leucine`, -->
<!-- #     y = unlist(predictWithModel( -->
<!-- #       data = metabolomics_data, -->
<!-- #       model_type = "linear_regression", -->
<!-- #       model = basic_regression_model$model -->
<!-- #     ))), -->
<!-- #     color = "maroon", size = 1 -->
<!-- #   ) + -->
<!-- #   geom_line(aes( -->
<!-- #     x = metabolomics_data$`iso-Leucine`, -->
<!-- #     y = unlist(predictWithModel( -->
<!-- #       data = metabolomics_data, -->
<!-- #       model_type = "linear_regression", -->
<!-- #       model = multiple_regression_model$model -->
<!-- #     ))), -->
<!-- #     color = "black", size = 1, alpha = 0.6 -->
<!-- #   ) + -->
<!-- #   theme_bw() -->

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
  output_variable = "Valine"
)

check_model(multiple_regression_model$model)
```

<img src="index_files/figure-html/unnamed-chunk-455-1.png" width="100%" style="display: block; margin: auto;" />

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
- `input_variables` selects the features or predictors for the model, here using columns 3 to 32 from metabolomics_data as predictors.
- `output_variable` is the target variable for prediction, in this case, "Valine".

The optimization_parameters argument takes a list to define the grid of parameters for optimization, including n_vars_tried_at_split, n_trees, and min_leaf_size. The seq() function generates sequences of numbers and is used here to create ranges for each parameter:

- `n_vars_tried_at_split` = seq(1,24,3) generates a sequence for the number of variables tried at each split, starting at 1, ending at 24, in steps of 3 (e.g., 1, 4, 7, ..., 24).
- `n_trees` = seq(1,40,2) creates a sequence for the number of trees in the forest, from 1 to 40 in steps of 2.
- `min_leaf_size` = seq(1,3,1) defines the minimal size of leaf nodes, ranging from 1 to 3 in steps of 1.

This setup creates a grid of parameter combinations where each combination of n_vars_tried_at_split, n_trees, and min_leaf_size defines a unique random forest model. The function will test each combination within this grid to identify the model that performs best according to a given evaluation criterion, effectively searching through a defined parameter space to optimize the random forest's performance. This approach allows for a systematic exploration of how different configurations affect the model's ability to predict the output variable, enabling the selection of the most effective model configuration based on the dataset and task at hand.


``` r
random_forest_model <- buildModel2(
    data = metabolomics_data,
    model_type = "random_forest_regression",
    input_variables = colnames(metabolomics_data)[3:32],
    output_variable = "Valine",
    optimization_parameters = list(
      n_vars_tried_at_split = seq(1,20,5),
      n_trees = seq(1,40,10),
      min_leaf_size = seq(1,3,1)
    )
)

names(random_forest_model)
## [1] "model_type" "metrics"    "model"
```

The above code builds our random forest model. It's output provides both the model itself and key components indicating the performance and configuration of the model. Here's a breakdown of each part of the output:

`$model_type` tells us what type of model this is.

`$model` shows the configuration of the best random forest model that was created.

`$metrics` provides detailed results of model performance across different combinations of the random forest parameters n_vars_tried_at_split (the number of variables randomly sampled as candidates at each split) and n_trees (the number of trees in the forest). For each combination, it shows:
- n_vars_tried_at_split and n_trees: The specific values used in that model configuration.
- .metric: The performance metric used, here it's accuracy, which measures how often the model correctly predicts the patient status.
- .estimator: Indicates the type of averaging used for the metric, here it's binary for binary classification tasks.
- mean: The average accuracy across the cross-validation folds.
- fold_cross_validation: Indicates the number of folds used in cross-validation, here it's 3 for all models.
- std_err: The standard error of the mean accuracy, providing an idea of the variability in model performance.
- .config: A unique identifier for each model configuration tested.

We can thus inspect the performance of the model based on the specific parameters used during configuration. This can help us understand if we are exploring the right parameter space - do we have good values for n_vars_tried_at_split and n_trees? In this case we are doing regression, and the performance metric reported is RMSE: root mean squared error. We want that value to be small! So smaller values for that metric indicate a better model.


``` r
random_forest_model$metrics %>%
    ggplot(aes(x = n_vars_tried_at_split, y = n_trees, fill = mean)) +
    facet_grid(.~min_leaf_size) +
    scale_fill_viridis(direction = -1) +
    geom_tile() +
    theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-458-1.png" width="100%" style="display: block; margin: auto;" />

We can easily use the model to make predictions by using the `predictWithModel()` function:


``` r
ggplot() +
  geom_point(
    data = metabolomics_data,
    aes(x = `iso-Leucine`, y = Valine), fill = "gold", shape = 21, color = "black"
  ) +
  geom_line(aes(
    x = metabolomics_data[["iso-Leucine"]],
    y = mean(metabolomics_data$Valine)
  ), color = "grey") +
  geom_line(aes(
    x = metabolomics_data[["iso-Leucine"]],
    y = unlist(predictWithModel(
      data = metabolomics_data,
      model_type = "random_forest_regression",
      model = random_forest_model$model
    ))),
    color = "maroon", size = 1
  ) +
  theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-459-1.png" width="100%" style="display: block; margin: auto;" />

In addition to regression modeling, random forests can also be used to do classification modeling. In classification modeling, we are trying to predict a categorical outcome variable from a set of predictor variables. For example, we might want to predict whether a patient has a disease or not based on their metabolomics data. All we have to do is set the model_type to "random_forest_classification" instead of "random_forest_regression". Let's try that now:


``` r
set.seed(123)
unknown <- metabolomics_data[35:40,]

rfc <- buildModel2(
  data = metabolomics_data[c(1:34, 41:93),],
  model_type = "random_forest_classification",
  input_variables = colnames(metabolomics_data)[3:75],
  output_variable = "patient_status",
  optimization_parameters = list(
    n_vars_tried_at_split = seq(5,50,5),
    n_trees = seq(10,100,10),
    min_leaf_size = seq(1,3,1)
  )
)
rfc$metrics %>% arrange(desc(mean))
## # A tibble: 300 × 9
##    n_vars_tried_at_split n_trees min_leaf_size .metric 
##                    <dbl>   <dbl>         <dbl> <chr>   
##  1                    25      20             2 accuracy
##  2                     5      50             1 accuracy
##  3                    10      60             2 accuracy
##  4                     5      80             2 accuracy
##  5                    10     100             2 accuracy
##  6                     5      60             3 accuracy
##  7                    10      70             3 accuracy
##  8                    30      20             1 accuracy
##  9                    25      50             2 accuracy
## 10                    20      90             2 accuracy
## # ℹ 290 more rows
## # ℹ 5 more variables: .estimator <chr>, mean <dbl>,
## #   std_err <dbl>, .config <chr>,
## #   fold_cross_validation <int>
```

Cool! Our best settings lead to a model with >90% accuracy! We can also make predictions on unknown data with this model:


``` r
rfc$metrics %>%
    ggplot(aes(x = n_vars_tried_at_split, y = n_trees, fill = mean)) +
    facet_grid(.~min_leaf_size) +
    scale_fill_viridis(direction = -1) +
    geom_tile() +
    theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-461-1.png" width="100%" style="display: block; margin: auto;" />


``` r
predictions <- predictWithModel(
  data = unknown,
  model_type = "random_forest_classification",
  model = rfc$model
)

data.frame(
  real_status = metabolomics_data[35:40,]$patient_status,
  predicted_status = unlist(predictions)
)
##                          real_status predicted_status
## .pred_healthy1               healthy             0.90
## .pred_healthy2               healthy             0.80
## .pred_healthy3               healthy             0.90
## .pred_healthy4        kidney_disease             0.10
## .pred_healthy5        kidney_disease             0.00
## .pred_healthy6        kidney_disease             0.05
## .pred_kidney_disease1        healthy             0.10
## .pred_kidney_disease2        healthy             0.20
## .pred_kidney_disease3        healthy             0.10
## .pred_kidney_disease4 kidney_disease             0.90
## .pred_kidney_disease5 kidney_disease             1.00
## .pred_kidney_disease6 kidney_disease             0.95
```

## {-}

## further reading {-}

- More on assessing regression models: [performance R package](https://github.com/easystats/performance). The performance package helps check how well your statistical models work by providing simple tools to evaluate things like fit quality and performance. It offers easy ways to spot problems in models, like if they're too complex or not fitting the data properly, and works with different types of models, including mixed-effects and Bayesian ones.

- [common machine learning tasks](https://pythonprogramminglanguage.com/machine-learning-tasks/). Machine learning involves using algorithms to learn from data, with key tasks including classification, regression, and clustering. Classification categorizes data, such as recognizing images of animals, regression predicts continuous values like sales forecasts, and clustering groups data based on similarities without prior labels.

Classification with random forests:
- http://www.rebeccabarter.com/blog/2020-03-25_machine_learning/
- https://hansjoerg.me/2020/02/09/tidymodels-for-machine-learning/
- https://towardsdatascience.com/dials-tune-and-parsnip-tidymodels-way-to-create-and-tune-model-parameters-c97ba31d6173


<!-- ## exercises {-}

To practice creating models, try the following:

1. Choose one of the datasets we have used so far, and run a principal components analysis on it. Note that the output of the analysis when you run "pca" contains the Dimension 1 coordinate "Dim.1" for each sample, as well as the abundance of each analyte in that sample.

2. Using the information from the ordination plot, identify two analytes: one that has a variance that is strongly and positively correlated with the first principal component (i.e. dimension 1), and one that has a variance that is slightly less strongly, but still positively correlated with the first principal component. Using `buildModel`, create and plot two linear regression models, one that regresses each of those analytes against dimension 1 (in other words, the x-axis should be the Dim.1 coordinate for each sample, and the y-axis should be the values for one of the two selected analytes). Which has the greater r-squared value? Based on what you know about PCA, does that make sense?

3. Choose two analytes: one should be one of the analytes from question 2 above, the other should be an analyte that, according to your PCA ordination analysis, is negatively correlated with the first principal component. Using `buildModel` and `predictWithModel` create plots showing how those two analytes are correlated with dimension 1. One should be positively correlated, and the other negatively correlated.

4. Have a look at the dataset `metabolomics_unknown`. It is metabolomics data from patients with an unknown healthy/kidney disease status. Build a classification model using the `metabolomics_data` dataset and diagnose each patient in `metabolomics_unknown`.

5. (optional) Explore the wine_quality dataset. What are the most important factors in determining what is a good wine? Use a model in an inferential way to provide evidence for your answer.


{r, echo = FALSE, eval = FALSE}
runMatrixAnalysis(
  data = wine_grape_data,
  analysis = "pca_ord",
  column_w_names_of_multiple_analytes = "metabolite",
  column_w_values_for_multiple_analytes = "log_abundance",
  columns_w_values_for_single_analyte = NULL,
  columns_w_additional_analyte_info = NULL,
  columns_w_sample_ID_info = c("cultivar", "treatment"),
  na_replacement = "drop"
) %>%
filter(Dim.1 > 0.9) %>%
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2)) +
  geom_label_repel(aes(x = Dim.1, y = Dim.2, label = analyte))

runMatrixAnalysis(
  data = wine_grape_data,
  analysis = "pca_ord",
  column_w_names_of_multiple_analytes = "metabolite",
  column_w_values_for_multiple_analytes = "log_abundance",
  columns_w_values_for_single_analyte = NULL,
  columns_w_additional_analyte_info = NULL,
  columns_w_sample_ID_info = c("cultivar", "treatment"),
  na_replacement = "drop"
) %>%
filter(Dim.1 < -0.9) %>%
ggplot() +
  geom_segment(aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2)) +
  geom_label_repel(aes(x = Dim.1, y = Dim.2, label = analyte))

out <- runMatrixAnalysis(
  data = wine_grape_data,
  analysis = "pca",
  column_w_names_of_multiple_analytes = "metabolite",
  column_w_values_for_multiple_analytes = "log_abundance",
  columns_w_values_for_single_analyte = NULL,
  columns_w_additional_analyte_info = NULL,
  columns_w_sample_ID_info = c("cultivar", "treatment"),
  na_replacement = "drop"
)

model <- buildLinearModel(
  data = solvents,
  # formula = "Dim.1 = Aspartic_acid_3TMS"
  # formula = "Dim.1 = Tyrosine_2TMS"
  formula = "refractive_index = formula_weight + melting_point"
)


top <- ggplot(model$data) +
  geom_point(aes(x = input_x, y = input_y)) +
  geom_line(aes(x = model_x, y = model_y)) +
  # annotate(geom = "table",
  #   x = 13,
  #   y = 16,
  #   label = list(model$metrics)
  # ) +
  # coord_cartesian(ylim = c(10,16)) +
  theme_bw()

bottom <- ggplot(model$data) +
  geom_col(
    aes(x = input_x, y = residuals),
    width = 0.03, color = "black", position = "dodge", alpha = 0.5
  ) +
  theme_bw()

cowplot::plot_grid(top, bottom, ncol = 1, labels = "AUTO", rel_heights = c(2,1))

model$metrics

``` -->


<!-- end -->

<!-- start embedding models -->


# embedding models {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/embedding.jpeg" width="100%" style="display: block; margin: auto;" />

To run the analyses in this chapter, you will need four things. 

1. Please ensure that your computer can run the following R script. It may prompt you to install additional R packages.

``` r
source("https://thebustalab.github.io/phylochemistry/modules/language_model_analysis.R")
## Loading language model module...
## Done with language model loading!
```
2. Please create an account at and obtain an API key from https://pubmed.ncbi.nlm.nih.gov/ (Login > Account Settings > API Key Management)
3. Please create an account at and obtain an API key from https://huggingface.co (Login > Settings > Access Tokens, then configure your access token/key to "Make calls to the serverless Inference API" and "Make calls to Inference Endpoints")
4. Please create an account at and obtain an API key from https://biolm.ai/ (Login > Account > API Tokens)
5. Please create an account (you may also need to create an NVIDIA cloud account if prompted) at and obtain an API key from https://build.nvidia.com/. (To get API key, go to: https://build.nvidia.com/meta/esm2-650m, switch "input" to python and click "Get API Key" > Generate Key)

Keep your API keys (long sequences of numbers and letters, like a password) handy for use in these analyses.

In the last chapter, we looked at models that use numerical data to understand the relationships between different aspects of a data set (inferential model use) and models that make predictions based on numerical data (predictive model use). In this chapter, we will explore a set of models called language models that transform non-numerical data (such as written text or protein sequences) into the numerical domain, enabling the non-numerical data to be analyzed using the  techniques we have already covered. Language models are algorithms that are trained on large amounts of text (or, in the case of protein language models, many sequences) and can perform a variety of tasks related to their training data. In particular, we will focus on embedding models, which convert language data into numerical data. An embedding is a numerical representation of data that captures its essential features in a lower-dimensional space or in a different domain. In the context of language models, embeddings transform text, such as words or sentences, into vectors of numbers, enabling machine learning models and other statistical methods to process and analyze the data more effectively. 

A basic form of an embedding model is a neural network called an autoencoder. Autoencoders consist of two main parts: an encoder and a decoder. The encoder takes the input data and compresses it into a lower-dimensional representation, called an embedding. The decoder then reconstructs the original input from this embedding, and the output from the decoder is compared against the original input. The model (the encoder and the decoder) are then iteratively optimized with the objective of minimizing a loss function that measures the difference between the original input and its reconstruction, resulting in an embedding model that creates meaningful embeddings that capture the important aspects of the original input.

## pre-reading {-}

Please read over the following:

- [Text Embeddings: Comprehensive Guide](https://towardsdatascience.com/text-embeddings-comprehensive-guide-afd97fce8fb5). In her article, "Text Embeddings: Comprehensive Guide", Mariya Mansurova explores the evolution, applications, and visualization of text embeddings. Beginning with early methods like Bag of Words and TF-IDF, she traces how embeddings have advanced to capture semantic meaning, highlighting significant milestones such as word2vec and transformer-based models like BERT and Sentence-BERT. Mansurova explains how these embeddings transform text into vectors that computers can analyze for tasks like clustering, classification, and anomaly detection. She provides practical examples using tools like OpenAI’s embedding models and dimensionality reduction techniques, making this article an in-depth resource for both theoretical and hands-on understanding of text embeddings.

- [ESM3: Simulating 500 million years of evolution with a language model](https://www.evolutionaryscale.ai/blog/esm3-release#simulating-500-million-years-of-evolution). The 2024 blog article "ESM3: Simulating 500 million years of evolution with a language model" by EvolutionaryScale introduces ESM3, a revolutionary language model trained on billions of protein sequences. This article explores how ESM3 marks a major advancement in computational biology by enabling researchers to reason over protein sequences, structures, and functions. With massive datasets and powerful computational resources, ESM3 can generate entirely new proteins, including esmGFP, a green fluorescent protein that differs significantly from known natural variants. The article highlights the model's potential to transform fields like medicine, synthetic biology, and environmental sustainability by making protein design programmable. Please note the "Open Model" section of the blog, which highlights applications of ESM models in the natural sciences.

## text embeddings {-}

Here, we will create text embeddings using publication data from PubMed. Text embeddings are numerical representations of text that preserve important information and allow us to apply mathematical and statistical analyses to textual data. Below, we use a series of functions to obtain titles and abstracts from PubMed, create embeddings for their titles, and analyze them using principal component analysis.

First, we use the searchPubMed function to extract relevant publications from PubMed based on specific search terms. This function interacts with the PubMed website via a tool called an API. An API, or Application Programming Interface, is a set of rules that allows different software programs to communicate with each other. In this case, the API allows our code to access data from the PubMed database directly, without needing to manually search through the website.  An API key is a unique identifier that allows you to authenticate yourself when using an API. It acts like a password, giving you permission to access the API services. Here, I am reading my API key from a local file. You can obtain by signing up for an NCBI account at https://pubmed.ncbi.nlm.nih.gov/. Once you have an API key, pass it to the searchPubMed function along with your search terms. Here I am using "beta-amyrin synthase," "friedelin synthase," "Sorghum bicolor," and "cuticular wax biosynthesis." I also specify that I want the results to be sorted according to relevance (as opposed to sorting by date) and I only want three results per term (the top three most relevant hits) to be returned:


``` r
search_results <- searchPubMed(
  search_terms = c("beta-amyrin synthase", "friedelin synthase", "sorghum bicolor", "cuticular wax biosynthesis"),
  pubmed_api_key = readLines("/Users/bust0037/Documents/Websites/pubmed_api_key.txt"),
  retmax_per_term = 3,
  sort = "relevance"
)
colnames(search_results)
## [1] "entry_number" "term"         "date"        
## [4] "journal"      "title"        "doi"         
## [7] "abstract"
select(search_results, term, title)
## # A tibble: 12 × 2
##    term                       title                         
##    <chr>                      <chr>                         
##  1 beta-amyrin synthase       β-Amyrin synthase from Conyza…
##  2 beta-amyrin synthase       Ginsenosides in Panax genus a…
##  3 beta-amyrin synthase       β-Amyrin synthase (EsBAS) and…
##  4 friedelin synthase         Friedelin in Maytenus ilicifo…
##  5 friedelin synthase         Friedelin Synthase from Mayte…
##  6 friedelin synthase         Functional characterization o…
##  7 sorghum bicolor            Sorghum (Sorghum bicolor).    
##  8 sorghum bicolor            Progress in Optimization of A…
##  9 sorghum bicolor            Grain and sweet sorghum (Sorg…
## 10 cuticular wax biosynthesis Regulatory mechanisms underly…
## 11 cuticular wax biosynthesis Cuticular wax in wheat: biosy…
## 12 cuticular wax biosynthesis Update on Cuticular Wax Biosy…
```

From the output here, you can see that we've retrieved records for various publications, each containing information such as the title, journal, and search term used. This gives us a dataset that we can further analyze to gain insights into the relationships between different research topics.

Next, we use the embedText function to create embeddings for the titles of the extracted publications. Just like PubMed, the Hugging Face API requires an API key, which acts as a unique identifier and grants you access to their services. You can obtain an API key by signing up at https://huggingface.co and following the instructions to generate your own key. Once you have your API key, you will need to specify it when using the embedText function. In the example below, I am reading the key from a local file for convenience.

To set up the embedText function, provide the dataset containing the text you want to embed (in this case, search_results, the output from the PubMed search above), the column with the text (title), and your Hugging Face API key. This function will then generate numerical embeddings for each of the publication titles. By default, the embeddings are generated using a pre-trained embedding language model called 'BAAI/bge-small-en-v1.5', available through the Hugging Face API at https://api-inference.huggingface.co/models/BAAI/bge-small-en-v1.5. This model is designed to create compact, informative numerical representations of text, making it suitable for a wide range of downstream tasks, such as clustering or similarity analysis. If you would like to know more about the model and its capabilities, you can visit the Hugging Face website at https://huggingface.co, where you will find detailed documentation and additional resources.


``` r
search_results_embedded <- embedText(
  df = search_results,
  column_name = "title",
  hf_api_key = readLines("/Users/bust0037/Documents/Websites/hf_api_key.txt")
)
##   |                                                          |                                                  |   0%  |                                                          |==================================================| 100%
search_results_embedded[1:3,1:10]
## # A tibble: 3 × 10
##   entry_number term  date       journal title doi   abstract
##          <dbl> <chr> <date>     <chr>   <chr> <chr> <chr>   
## 1            1 beta… 2019-11-20 FEBS o… β-Am… 10.1… Conyza …
## 2            2 beta… 2024-04-03 Acta p… Gins… 10.1… Ginseno…
## 3            3 beta… 2017-03-16 Phytoc… β-Am… 10.1… Siberia…
## # ℹ 3 more variables: embedding_1 <dbl>, embedding_2 <dbl>,
## #   embedding_3 <dbl>
```

The output of the embedText function is a data frame where the 384 appended columns represent the embedding variables. These embeddings capture the features of each publication title. These embeddings are like a bar codes:


``` r
search_results_embedded %>%
  pivot_longer(
    cols = grep("embed",colnames(search_results_embedded)),
    names_to = "embedding_variable",
    values_to = "value"
  ) %>%
  ggplot() +
    geom_tile(aes(x = embedding_variable, y = factor(entry_number), fill = value)) +
    scale_y_discrete(name = "article") +
    scale_fill_gradient(low = "white", high = "black") +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
```

<img src="index_files/figure-html/unnamed-chunk-493-1.png" width="100%" style="display: block; margin: auto;" />

To examine the relationships between the publication titles, we perform PCA on the text embeddings. We use the runMatrixAnalysis function, specifying PCA as the analysis type and indicating which columns contain the embedding values. We visualize the results using a scatter plot, with each point representing a publication title, colored by the search term it corresponds to. The `grep` function is used here to search for all column names in the `search_results` data frame that contain the word 'embed'. This identifies and selects the columns that hold the embedding values, which will be used as the columns with values for single analytes for the PCA and enable the visualization below. While we've seen lots of PCA plots over the course of our explorations, note that this one is different in that it represents the relationships between the meaning of text passages (!) as opposed to relationships between samples for which we have made many measurements of numerical attributes.


``` r
runMatrixAnalysis(
  data = search_results_embedded,
  analysis = "pca",
  columns_w_values_for_single_analyte = colnames(search_results_embedded)[grep("embed", colnames(search_results_embedded))],
  columns_w_sample_ID_info = c("title", "journal", "term")
) %>%
  ggplot() +
    geom_label_repel(
      aes(x = Dim.1, y = Dim.2, label = str_wrap(title, width = 35)),
      size = 2, min.segment.length = 0.5, force = 50
    ) +  
    geom_point(aes(x = Dim.1, y = Dim.2, fill = term), shape = 21, size = 5, alpha = 0.7) +
    scale_fill_brewer(palette = "Set1") +
    scale_x_continuous(expand = c(0,1)) +
    scale_y_continuous(expand = c(0,5)) +
    theme_minimal()
```

<img src="index_files/figure-html/unnamed-chunk-494-1.png" width="100%" style="display: block; margin: auto;" />

We can also use embeddings to examine data that are not full sentences but rather just lists of terms, such as the descriptions of odors in the `beer_components` dataset:


``` r
n <- 31

odor <- data.frame(
  sample = seq(1,n,1),
  odor = dropNA(unique(beer_components$analyte_odor))[sample(1:96, n)]
)

out <- embedText(
  odor, column_name = "odor",
  hf_api_key = readLines("/Users/bust0037/Documents/Websites/hf_api_key.txt")
)
##   |                                                          |                                                  |   0%  |                                                          |=========================                         |  50%  |                                                          |==================================================| 100%

runMatrixAnalysis(
  data = out,
  analysis = "pca",
  columns_w_values_for_single_analyte = colnames(out)[grep("embed", colnames(out))],
  columns_w_sample_ID_info = c("sample", "odor")
) -> pca_out

pca_out$color <- rgb(
  scales::rescale(pca_out$Dim.1, to = c(0, 1)),
  0,
  scales::rescale(pca_out$Dim.2, to = c(0, 1))
)

ggplot(pca_out) +
  geom_label_repel(
    aes(x = Dim.1, y = Dim.2, label = str_wrap(odor, width = 35)),
    size = 2, min.segment.length = 0.5, force = 25
  ) +  
  geom_point(aes(x = Dim.1, y = Dim.2), fill = pca_out$color, shape = 21, size = 3, alpha = 0.7) +
  # scale_x_continuous(expand = c(1,0)) +
  # scale_y_continuous(expand = c(1,0)) +
  theme_minimal()
```

<img src="index_files/figure-html/unnamed-chunk-495-1.png" width="100%" style="display: block; margin: auto;" />

## protein embeddings {-}

<!-- Protein Language Models: -->
<!--     The goal in these models is to train them so that the embeddings they create capture important biological features of proteins. -->
<!--     The attention mechanism in transformer models allows capturing both local and global information in a protein sequence: -->
<!--         Local Information: Might include interactions between neighboring amino acids. -->
<!--         Global Information: Could encompass long-range relationships between distant parts of the sequence. -->
<!--     While embedding models can be simple autoencoders, many embedding models, especially in protein language modeling, use transformers with attention mechanisms to capture complex patterns in the data. -->

<!-- Attention Mechanism: -->
<!--     The attention mechanism works within the encoder and decoder, allowing each element of the input (e.g., an amino acid) to compare itself to every other element. -->
<!--     It generates attention scores to weigh how much attention one amino acid should give to another. -->
<!--     The attention mechanism helps capture both local and long-range dependencies in protein sequences, enabling the model to focus on important areas regardless of their position in the sequence. -->

<!-- Why Attention is Beneficial: -->
<!--     Long-Range Dependencies: Captures interactions between distant amino acids. -->
<!--     Structural Complexity: Weighs relationships between amino acids to account for protein folding and interactions. -->
<!--     Handling Variable Sequence Lengths: Adjusts focus across sequences of varying lengths. -->
<!--     Multi-Dimensional Relationships: Multi-head attention allows capturing different kinds of relationships, like hydrophobic interactions or secondary structures. -->
<!--     Contextualized Embeddings: Embeddings reflect the broader sequence environment, not just local motifs. -->

<!-- Additional Mechanisms in Protein Language Models: -->
<!--     Positional Encoding: Adds position information to the sequence so that the model can differentiate between identical amino acids at different positions. -->
<!--     Masked Language Modeling (MLM): Trains the model to predict masked amino acids, learning patterns in the sequence. -->
<!--     Multiscale Representations: Allows capturing both fine-grained and coarse-grained structural information. -->
<!--     Evolutionary Information: Incorporates multiple sequence alignments (MSAs) to learn from conserved regions. -->
<!--     Residual Connections: Helps information flow through the network and stabilizes training by allowing the model to retain original input data as it processes through layers. -->
<!--     Normalization and Regularization: Techniques like layer normalization and dropout are used to stabilize training and prevent overfitting. -->

Autoencoders can be trained to accept various types of inputs, such as text (as shown above), images, audio, videos, sensor data, and sequence-based information like peptides and DNA. Protein language models convert protein sequences into numerical representations that can be used for a variety of downstream tasks, such as structure prediction or function annotation. Protein language models, like their text counterparts, are trained on large datasets of protein sequences to learn meaningful patterns and relationships within the sequence data.

Protein language models offer several advantages over traditional approaches, such as multiple sequence alignments (MSAs). One major disadvantage of MSAs is that they are computationally expensive and become increasingly slow as the number of sequences grows. While language models are also computationally demanding, they are primarily resource-intensive during the training phase, whereas applying a trained language model is much faster. Additionally, protein language models can capture both local and global sequence features, allowing them to identify complex relationships that span across different parts of a sequence. Furthermore, unlike MSAs, which rely on evolutionary information, protein language models can be applied to proteins without homologous sequences, making them suitable for analyzing sequences where little evolutionary data is available. This flexibility broadens the scope of proteins that can be effectively studied using these models.

Beyond the benefits described above, protein language models have an additional, highly important capability: the ability to capture information about connections between elements in their input, even if those elements are very distant from each other in the sequence. This capability is achieved through the use of a model architecture called a transformer, which is a more sophisticated version of an autoencoder. For example, amino acids that are far apart in the primary sequence may be very close in the 3D, folded protein structure. Proximate amino acids in 3D space can play crucial roles in protein stability, enzyme catalysis, or binding interactions, depending on their spatial arrangement and interactions with other residues. Embedding models with transformer architecture can effectively capture these functionally important relationships.

By adding a mechanism called an "attention mechanism" to an autoencoder, we can create a simple form of a transformer. The attention mechanism works within the encoder and decoder, allowing each element of the input (e.g., an amino acid) to compare itself to every other element, generating attention scores that weigh how much attention one amino acid should give to another. This mechanism helps capture both local and long-range dependencies in protein sequences, enabling the model to focus on important areas regardless of their position in the sequence. Attention is beneficial because it captures interactions between distant amino acids, weighs relationships to account for protein folding and interactions, adjusts focus across sequences of varying lengths, captures different types of relationships like hydrophobic interactions or secondary structures, and provides contextualized embeddings that reflect the broader sequence environment rather than just local motifs. For more on attention mechanisms, check out the further reading section of this chapter.

In this section, we will explore how to generate embeddings for protein sequences using a pre-trained protein language model and demonstrate how these embeddings can be used to analyze and visualize protein data effectively. First, we need some data. You can use the `OSC_sequences` object provided by the `source()` code, though you can also use the `searchNCBI()` function to retrieve your own sequences. For example:


``` r
ncbi_results <- searchNCBI(search_term = "oxidosqualene cyclase", retmax = 100)
ncbi_results
## AAStringSet object of length 100:
##       width seq                         names               
##   [1]   513 MVPYALPIHPGR...LGEFRRRLLANK CAN6220183.1 unna...
##   [2]   586 MFGSVLTYVSLR...EYRCRVLAAGKQ CAN6181702.1 unna...
##   [3]   329 MTDTDGTLAWLL...WDAAAAQPQPPR WP_289330907.1 MU...
##   [4]   575 MNKRNEIEEMIR...LNALKKVQTVID WP_446787674.1 pr...
##   [5]   290 MKDTLSIKLFKT...YTFYGLLALGTI WP_446787673.1 pr...
##   ...   ... ...
##  [96]   735 MSGTHIAPWRTP...LYARKFGNDSLL XP_020051874.1 te...
##  [97]   717 MADPIAPWRTAA...LYSRKFGNEELL XP_015409667.1 te...
##  [98]   700 MPPTLPEKTDYT...KFARTYPDYRLT XP_015402884.1 te...
##  [99]   749 MPDHIGRWKNGG...LYSRKFGNEELQ XP_681529.2 terpe...
## [100]   735 MAADHIGPWRTD...LYSRKFGNEELI XP_043139417.1 te...
```

Once you have some sequences, we can embed them with the function `embedAminoAcids()`. An example is below. Note that we need to provide either a biolm API key or an NVIDIA api key, and specify which platform we wish to use. We also need to provide the amino acid sequences as an AAStringSet object. If you use the NVIDIA platform, the model esm2-650m will be used (note: esm2 truncates sequences longer than 1022 AA in length). If you use bioLM, you can pick between a number of models.


``` r
embedded_OSCs <- embedAminoAcids(
  amino_acid_stringset = OSC_sequences,
  biolm_api_key = readLines("/Users/bust0037/Documents/Websites/biolm_api_key.txt"),
  nvidia_api_key = readLines("/Users/bust0037/Documents/Websites/nvidia_api_key.txt"),
  platform = "biolm"
)
embedded_OSCs$product <- tolower(gsub(".*_", "", embedded_OSCs$name))
embedded_OSCs <- select(embedded_OSCs, name, product, everything())
embedded_OSCs[1:3,1:4]
```

Nice! Once we've bot the embeddings, we can run a PCA analysis to visualize them in 2D space:


``` r
runMatrixAnalysis(
  data = embedded_OSCs,
  analysis = "pca",
  columns_w_values_for_single_analyte = colnames(embedded_OSCs)[3:dim(embedded_OSCs)[2]],
  columns_w_sample_ID_info = c("name", "product")
) %>%
  ggplot() +
    geom_jitter(
      aes(x = Dim.1, y = Dim.2, fill = product),
      shape = 21, size = 5, height = 2, width = 2, alpha = 0.6
    ) +
    theme_minimal()

```

## {-}

## further reading {-}

- [creating knowledge graphs with LLMs](https://bratanic-tomaz.medium.com/constructing-knowledge-graphs-from-text-using-openai-functions-096a6d010c17). This blog post explains how to create knowledge graphs from text using OpenAI functions combined with LangChain and Neo4j. It highlights how large language models (LLMs) have made information extraction more accessible, providing step-by-step instructions for setting up a pipeline to extract structured information and construct a graph from unstructured data.

- [creating RAG systems with LLMs](https://medium.com/enterprise-rag/a-first-intro-to-complex-rag-retrieval-augmented-generation-a8624d70090f). This article provides a technical overview of implementing complex Retrieval Augmented Generation (RAG) systems, focusing on key concepts like chunking, query augmentation, document hierarchies, and knowledge graphs. It highlights the challenges in data retrieval, multi-hop reasoning, and query planning, while also discussing opportunities to improve RAG infrastructure for more accurate and efficient information extraction.

- [using protein embeddings in biochemical research](https://www.biorxiv.org/content/10.1101/2024.01.29.577750v3). This study presents a machine learning pipeline that successfully identifies and characterizes terpene synthases (TPSs), a challenging task due to the limited availability of labeled protein sequences. By combining a curated TPS dataset, advanced structural domain segmentation, and language model techniques, the authors discovered novel TPSs, including the first active enzymes in Archaea, significantly improving the accuracy of substrate prediction across TPS classes.

- [attention mechanims and transformers explained](https://ig.ft.com/generative-ai/). This Financial Times article explains the development and workings of large language models (LLMs), emphasizing their foundation on the transformer model created by Google researchers in 2017. These models use self-attention mechanisms to understand context, allowing them to respond to subtle relationships between elements in their input, even if those elements are far from one another in the linear input sequence.

- [other types of protein language models](https://build.nvidia.com/nim?q=protein). *3D Protein Structure Prediction* deepmind / alphafold2-multimer: Predicts the 3D structure of protein complexes from amino acid sequences. deepmind / alphafold2: Predicts the 3D structure of single proteins from amino acid sequences. meta / esmfold: Predicts the 3D structure of proteins based on amino acid sequences. *Protein Embedding Generation* meta / esm2-650m: Generates protein embeddings from amino acid sequences. *Protein Sequence Design* ipd / proteinmpnn: Predicts amino acid sequences for given protein backbone structures. *Generative Protein Design* ipd / rfdiffusion: A generative model for designing protein backbones, particularly for protein binder design. *Molecule-Protein Interaction Prediction* mit / diffdock: Predicts the 3D interactions between molecules and proteins (docking simulations).

<!-- ## exercises {-}

1. Recreate the PubMed search and subsequent analysis described in this chapter using search terms that relate to research you are involved in or are interested in. Use multiple search terms and retrieve publications over a period of several years (you may need to set `sort` = "date"). Embed the titles and visualize the changes in clustering over time using PCA or an x-axis that is the date. Discuss how research trends might evolve and reflect broader changes in the scientific community or societal challenges. Below is an example to help you:


``` r
search_results_ex <- searchPubMed(
  search_terms = c("oxidosqualene cyclase", "chemotaxonomy", "protein engineering"),
  pubmed_api_key = readLines("/Users/bust0037/Documents/Science/Websites/pubmed_api_key.txt"),
  retmax_per_term = 50,
  sort = "date"
)

search_results_ex_embed <- embedText(
  search_results_ex, column_name = "abstract",
  hf_api_key = readLines("/Users/bust0037/Documents/Science/Websites/hf_api_key.txt")
)

runMatrixAnalysis(
  data = search_results_ex_embed,
  analysis = "pca",
  columns_w_values_for_single_analyte = colnames(search_results_ex_embed)[grep("embed", colnames(search_results_ex_embed))],
  columns_w_sample_ID_info = c("title", "journal", "term", "date")
) -> search_results_ex_embed_pca

search_results_ex_embed_pca %>%
    ggplot() +
      geom_point(aes(x = Dim.1, y = date, fill = date, shape = term), size = 5, alpha = 0.7) +
    scale_shape_manual(values = c(21, 22, 23)) +
    scale_fill_viridis() +
    scale_x_continuous(expand = c(0,1)) +
    scale_y_continuous(expand = c(0.1,0)) +
    theme_minimal()
```


2. Using the hops_components dataset, determine whether there are any major clusters of hops that are grouped by aroma. To do this, compute embeddings for the hop_aroma column of the dataset, then use a dimensional reduction (pca, if you like) to determine if any clear clusters are present.



3. Generate and visualize a set of protein embeddings. You can use `OSC_sequences` dataset provided by the source() command, or you can create your own protein sequence dataset using the `searchNCBI()` function.


asdf -->


________________________________________________________________________________________________
________________________________________________________________________________________________
________________________________________________________________________________________________

# (PART) SEQUENCE ANALYSIS {-}


# homology {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/homology2.png" width="100%" style="display: block; margin: auto;" />

## blastNCBI {-}

The `blastNCBI()` function provides a quick, programmatic way to search the NCBI database with custom sequences. By inputting a DNA sequence, the desired BLAST program (such as blastn or megablast), and selecting a database like "core_nt," you can initiate a remote BLAST search directly through NCBI’s servers. The function then handles the submission, waits for the search to complete, and retrieves the results for further analysis. The output is a tibble with details about each significant hit, allowing you to analyze identity percentages, gaps, and scores for matches in an accessible format. This function is useful for researchers needing automated access to NCBI BLAST results directly in R without needing to manually manage BLAST requests or wait on a browser. These examples assume the `Biostrings` package is available because it supplies `DNAStringSet`.


``` r
blastNCBI(
    query_sequence = "AAGCTTGGAAATATTAAGTGAACAGGGAATAGAAAGGATACAACAAAAGGGAAGAACTTAGAGCA",
    program = "blastn",
    database = "core_nt"
)

seqs <- DNAStringSet("AAGCTTGGAAATATTAAGTGAACAGGGAATAGAAAGGATACAACAAAAGGGAAGAACTTAGAGCA")
blastNCBI(
    query_sequence = as.character(seqs[1]),
    program = "blastn",
    database = "core_nt"
)

```

## local blast {-}

On NCBI, you can search various sequence collections with one or more queries. However, often we want to search a custom library, or multiple libraries. For example, maybe we have downloaded some genomes of interest and want to run blast searches on them. That is what polyBlast() is designed to do. polyBlast() relies on the BLAST+ program available from [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download). Download the program and then point this function to the executable via the `blast_module_directory_path` argument. You can search multiple sequence libraries at once using multiple queries, and all the usual blast configurations (blastp, blastn, tblastn, etc.) are available. Please note that searches with protein sequences or translated DNA sequences are 5–10-fold more sensitive than DNA:DNA sequence comparison.

Let's check out `polyBlast()` by looking at an example. For this example, we need to set up a few things:

* "named_subjects_list": A named list of sequence collections (often transcriptomes) to search (one fasta for each collection, often one collection for each species or accession).
* "query_in_path": One or more queries, all listed in a single fasta file.
* "sequences_of_interest_directory_path": The path to a directory where the BLAST hits will be written as individual files (this will be useful later on).
* "blast_module_directory_path": The path to the folder of BLAST+ executable.
* "blast_mode": The default shorthand is XYblastZ where X describes the subject type (`n` for nucleotide, `p` for amino acid), Y indicates whether subjects are translated (`t`) or not (`n`), and Z describes the query type (`n` or `p`). For example, `ntblastp` translates nucleotide subjects and compares them to protein queries. You can also supply specific BLAST+ presets such as `megablast`, `dc-megablast`, `blastx`, `tblastn`, or `tblastx` when you want BLAST to manage the translation behaviour for you.
* "e_value_cutoff": Hits with an e-value greater than this cutoff are discarded. Default = 1.
* "queries_in_output": TRUE/FALSE, should the queries be included in the output? Set to TRUE if you need the query sequences for downstream analyses like tree building.
* "monolist_out_path": The path to where we want a summary file of the BLAST hits to be written.

Different BLAST modes tailor the search to the molecule types you are comparing. In the XYblastZ naming scheme, `nnblastn` reproduces the classic nucleotide-vs-nucleotide search (`blastn`), `ntblastp` lines up with `tblastn`, and `pnblastp` matches `blastp`. Preset names such as `blastx`, `megablast`, or `dc-megablast` extend these basics to translated queries or alternate nucleotide heuristics, while `tblastx` translates both query and subject to compare inferred proteins. Pick the mode that matches the biology of your query and target sequences so the scoring and heuristics work in your favor.

Once we have those things, we can set up the search (see below). There are two main outputs from the search: a list of the hits ("monolist_out", which is written to "monolist_out_path"), and the hits themselves, written as individual files to "sequences_of_interest_directory_path". These two things can be used in downstream analyses, such as alignments. The function does not return an object.


``` r
the_transcriptomes <- c(
  "/path_to/the_transcriptomes_or_proteomes/Nicotiana_glauca.fa",
  "/path_to/the_transcriptomes_or_proteomes/Nicotiana_tabacum.fa",
  "/path_to/the_transcriptomes_or_proteomes/Nicotiana_benthamiana.fa"
)

names(the_transcriptomes) <- c(
  "Nicotiana_glauca",
  "Nicotiana_tabacum",
  "Nicotiana_benthamiana"
)

polyBlast(
  named_subjects_list = the_transcriptomes,
  query_in_path = "/path_to/sequences_you_want_to_find_in_the_transcriptomes.fa",
  sequences_of_interest_directory_path = "/path_to/a_folder_for_hit_sequences/",
  blast_module_directory_path = "/path_to/the_blast_module/", # on bustalab server this is /project_data/shared/general_lab_resources/blast/
  blast_mode = c("nnblastn", "ntblastp", "pnblastp", "dc-megablast"), 
  e_value_cutoff = 1,
  queries_in_output = TRUE,
  monolist_out_path = "/path_to/a_csv_file_that_will_list_all_blast_hits.csv"
)
```

## hmmer {-}

HMM, which stands for Hidden Markov Model, is a statistical model often used in various applications involving sequences, including speech recognition, natural language processing, and bioinformatics. In the context of bioinformatics, HMMs are frequently applied for sequence similarity searching, notably in the analysis of protein or DNA sequences. When we talk about using HMMs for sequence similarity searching, we're often referring to identifying conserved patterns or domains within biological sequences. These conserved regions can be indicative of functional or structural properties of the molecule. One of the advantages of using HMMs over traditional sequence similarity searching tools like BLAST is that HMMs can be more sensitive in detecting distant homologues. They take into account the position-specific variability within a protein family, as opposed to just looking for stretches of similar sequence.

Here's a general idea of how HMMs are used for sequence similarity searching:

1. Build a library of HMM domains: In bioinformatics, a typical application is the construction of library of HMM domains. These are HMMs built from multiple sequence alignments of a family of related proteins or genes. The alignments help highlight the conserved and variable positions in the sequence family. Once you have an alignment, the HMM can be 'trained' on this data. The training process estimates the probabilities of different events, like a particular amino acid (in the case of proteins) appearing at a specific position.

2. Predict domains in unknown sequences: After training, you can then use the HMMs to score other sequences. If a sequence scores above a certain threshold, it suggests that the sequence may be a member of the protein or gene family represented by the HMM. You can search databases of uncharacterized sequences using the HMM. Sequences in the database that get a high score against the HMM are potential new members of the family, and thus might share similar functional or structural properties.

We can implement these two steps using the `buildDomainLibrary()` function and the `predictDomains()` function. See below:


``` r
buildDomainLibrary(
    alignment_in_paths = c(
        "/project_data/shared/kalanchoe_transporters/alignments/subset_cluster_1_amin_seqs_aligned.fa",
        "/project_data/shared/kalanchoe_transporters/alignments/subset_cluster_2_amin_seqs_aligned.fa"
    ),
    domain_library_out_path = "/project_data/shared/kalanchoe_transporters/test.hmm"
)

predictDomains(
    fasta_in_path = "/project_data/shared/kalanchoe_transporters/alignments/subset_cluster_3_amin_seqs.fa",
    domain_library_in_path = "/project_data/shared/kalanchoe_transporters/test.hmm"
)
```

<!-- ## 3D similarity {-} -->

<!-- `foldseek easy-search /project_data/shared/kalanchoe_phylogeny/protein_structures/structures/AHF22083.1.pdb /project_data/shared/kalanchoe_phylogeny/protein_structures/structures/ result.m8 /project_data/shared/kalanchoe_phylogeny/protein_structures/tmp --exhaustive-search 1` -->

## interpreting homology data {-}

Bitscore versus e-value: when the expected value becomes extremely small (for example < 1e-250), R prints it as zero and you lose the ability to rank hits. In those cases rely on the bit score, which remains informative across the full range of alignments.

The "30% identity rule-of-thumb" is too conservative. Statistically significant (E < 10−6 – 10−3) protein homologs can share less than 20% identity. E-values and bit scores (bits > 50) are far more sensitive and reliable than percent identity for inferring homology.

The expect value (E-value) can be changed in order to limit the number of hits to the most significant ones. The lower the E-value, the better the hit. The E-value is dependent on the length of the query sequence and the size of the database. For example, an alignment obtaining an E-value of 0.05 means that there is a 5 in 100 chance of occurring by chance alone. E-values are very dependent on the query sequence length and the database size. Short identical sequence may have a high E-value and may be regarded as "false positive" hits. This is often seen if one searches for short primer regions, small domain regions etc. The default threshold for the E-value on the BLAST web page is 10, the default for polyBlast is 1. Increasing this value will most likely generate more hits. Below are some rules of thumb from the CLC Genomics Workbench documentation—treat them as approximate, database-dependent guidelines:

* E-value < 1e-100 Identical sequences. You will get long alignments across the entire query and hit sequence.
* 1e-100 < E-value < 1e-50 Almost identical sequences. A long stretch of the query protein is matched to the database.
* 1e-50 < E-value < 1e-10 Closely related sequences, could be a domain match or similar.
* 1e-10 < E-value < 1 Could be a true homologue but it is a gray area.
* E-value > 1 Proteins are most likely not related
* E-value > 10 Hits are most likely junk unless the query sequence is very short.

reference: https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/_E_value.html

reference: Pearson W. R. (2013). An introduction to sequence similarity ("homology") searching. Current protocols in bioinformatics, Chapter 3, Unit3.1. https://doi.org/10.1002/0471250953.bi0301s42.

* E-values and Bit-scores: Pfam-A is based around hidden Markov model (HMM) searches, as provided by the HMMER3 package. In HMMER3, like BLAST, E-values (expectation values) are calculated. The E-value is the number of hits that would be expected to have a score equal to or better than this value by chance alone. A good E-value is much less than 1. A value of 1 is what would be expected just by chance. In principle, all you need to decide on the significance of a match is the E-value.

E-values are dependent on the size of the database searched, so we use a second system in-house for maintaining Pfam models, based on a bit score (see below), which is independent of the size of the database searched. For each Pfam family, we set a bit score gathering (GA) threshold by hand, such that all sequences scoring at or above this threshold appear in the full alignment. It works out that a bit score of 20 equates to an E-value of approximately 0.1, and a score 25 of to approximately 0.01. From the gathering threshold both a “trusted cutoff” (TC) and a “noise cutoff” (NC) are recorded automatically. The TC is the score for the next highest scoring match above the GA, and the NC is the score for the sequence next below the GA, i.e. the highest scoring sequence not included in the full alignment.

* Sequence versus domain scores: There’s an additional wrinkle in the scoring system. HMMER3 calculates two kinds of scores, the first for the sequence as a whole and the second for the domain(s) on that sequence. The “sequence score” is the total score of a sequence aligned to the model (the HMM); the “domain score” is the score for a single domain — these two scores are virtually identical where only one domain is present on a sequence. Where there are multiple occurrences of the domain on a sequence any individual match may be quite weak, but the sequence score is the sum of all the individual domain scores, since finding multiple instances of a domain increases our confidence that that sequence belongs to that protein family, i.e. truly matches the model.

* Meaning of bit-score for non-mathematicians: A bit score of 0 means that the likelihood of the match having been emitted by the model is equal to that of it having been emitted by the Null model (by chance). A bit score of 1 means that the match is twice as likely to have been emitted by the model than by the Null. A bit score of 2 means that the match is 4 times as likely to have been emitted by the model than by the Null. So, a bit score of 20 means that the match is 2 to the power 20 times as likely to have been emitted by the model than by the Null.


# alignments {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/alignments.jpeg" width="100%" style="display: block; margin: auto;" />

Multiple sequence alignments form the bridge between raw sequences and downstream phylogenetic analysis. In this chapter we first use `alignSequences()` to assemble consistent nucleotide, amino-acid, or codon alignments. This is easily accomplished from the monolist generated by `polyBlast()`, if your sequences were generated that way. Then we turn to `analyzeAlignment()` to diagnose gap-heavy or poorly conserved sites, interactively trim them away, and export a cleaned alignment for tree building.

## alignSequences {-}

There are, of course, many tools for aligning sequences. `alignSequences()`, from the phylochemistry toolkit, is designed to be flexible: it will align nucleotide, amino-acid, or codon sequences, and it can restrict the alignment to any subset of records in your BLAST monolist. If you already ran `polyBlast()`, most of the inputs outlined below should be ready to go. The function writes the aligned FASTA to `alignment_out_path` and returns nothing (invisibly).

* `"monolist"`: a data frame where each row represents a hit of interest. It must contain an `accession` column whose values match individual FASTA files in `sequences_of_interest_directory_path`. The helper `readMonolist()` returns this object from the CSV generated by `polyBlast()`.
* `"subset"`: the name of a logical column in the monolist (for example `subset_all`). Only rows where that column is `TRUE` are aligned. Create additional `subset_*` columns if you want multiple alignment configurations.
* `"alignment_out_path"`: full file path (including filename) for the alignment that will be written, e.g. `/tmp/my_alignment.fa`. Existing files at this path are overwritten.
* `"sequences_of_interest_directory_path"`: directory that stores one FASTA per accession, as produced by `polyBlast()`. The function normalises the path internally, so either absolute or relative locations are fine.
* `"input_sequence_type"`: set to `"nucl"` when the subject FASTA files contain nucleotide sequences and `"amin"` for amino-acid sequences. Codon alignments start from nucleotides, while amino-acid alignments can start from either input type (translation happens automatically when needed).
* `"mode"`: chooses the alignment strategy. `"nucl_align"` performs a straight nucleotide multiple sequence alignment. `"amin_align"` aligns amino-acid sequences (translated from nucleotides if required). `"codon_align"` translates to protein, aligns there, and then projects the alignment back to nucleotides. `"fragment_align"` is reserved for fragment-to-fragment workflows and currently expects a base fragment defined via `base_fragment`.
* `"base_fragment"`: optional path to a FASTA file containing the fragment that other sequences should be aligned against. Only used when `mode = "fragment_align"`.


``` r
alignSequences(
  monolist = readMonolist("/path_to/a_csv_file_that_will_list_all_blast_hits.csv"),
  subset = "subset_all",
  alignment_out_path = "/path_to/a_folder_for_alignments/subset_all_amin_seqs_aligned.fa",
  sequences_of_interest_directory_path = "/path_to/a_folder_for_hit_sequences/",
  input_sequence_type = "amin",
  mode = "amin_align",
  base_fragment = NULL
)
```

## analyzeAlignment {-}

Once you have an alignment on disk, `analyzeAlignment()` helps you decide which positions to keep before phylogeny building. It reads a FASTA alignment, profiles each column for gap content and conservation, and launches a Shiny interface that lets you tune thresholds while inspecting the original and trimmed alignments alongside neighbour-joining trees. The function returns the trimmed alignment (as a `DNAStringSet` or `AAStringSet`) invisibly and also writes it to `<alignment_in_path>_trimmed`.

* `"alignment_in_path"`: path to the FASTA alignment to review. This is typically the file produced by `alignSequences()`, e.g. `/tmp/subset_all_amin_seqs_aligned.fa`.
* `"type"`: `"DNA"` for nucleotide alignments or `"AA"` for amino acids. This determines how the alignment is read, which tree scaffold is built, and which Biostrings container is used when writing the trimmed alignment.
* `"jupyter"`: set to `TRUE` when launching from a Jupyter environment so the app claims an available internal port and prints the connection URL. Leave at the default `FALSE` for regular RStudio or terminal sessions.

The interface plots gap percentages and conservation across the alignment, highlights the positions that survive the current filters, and previews the trimmed tree. When you click **Return Trimmed Alignment & Close App**, the filtered alignment is written to `<alignment_in_path>_trimmed` (with `_trimmed` appended) and returned to the caller.


``` r
analyzeAlignment(
  alignment_in_path = "/path_to/a_folder_for_alignments/subset_all_amin_seqs_aligned.fa",
  type = "AA",
  jupyter = FALSE
)
```


# phylogenies {-}

## buildTree {-}

This function is a swiss army knife for tree building. It takes as input alignments or existing phylogenies from which to derive a phylogeny of interest, it can use neighbor-joining or maximum liklihood methods (with model optimization), it can run bootstrap replicates, and it can calculate ancestral sequence states. To illustrate, let's look at some examples:

### newick input {-}

Let's use the Busta lab's plant phylogeny [derived from Qian et al., 2016] to build a phylogeny with five species in it.


``` r
tree <- buildTree(
  scaffold_type = "newick",
  scaffold = "https://thebustalab.github.io/data/plant_phylogeny.newick",
  members = c("Sorghum_bicolor", "Zea_mays", "Setaria_viridis", "Arabidopsis_thaliana", "Amborella_trichopoda")
)
## Pro tip: most tree read/write functions reset node numbers. Fortify your tree and save it as a csv file to preserve node numbering.

tree
## 
## Phylogenetic tree with 5 tips and 4 internal nodes.
## 
## Tip labels:
##   Amborella_trichopoda, Zea_mays, Sorghum_bicolor, Setaria_viridis, Arabidopsis_thaliana
## Node labels:
##   , , , 
## 
## Rooted; includes branch length(s).

plot(tree)
```

<img src="index_files/figure-html/unnamed-chunk-534-1.png" width="100%" style="display: block; margin: auto;" />

Cool! We got our phylogeny. What happens if we want to build a phylogeny that has a species on it that isn't in our scaffold? For example, what if we want to build a phylogeny that includes *Arabidopsis neglecta*? We can include that name in our list of members:


``` r
tree <- buildTree(
  scaffold_type = "newick",
  scaffold_in_path = "https://thebustalab.github.io/data/plant_phylogeny.newick",
  members = c("Sorghum_bicolor", "Zea_mays", "Setaria_viridis", "Arabidopsis_neglecta", "Amborella_trichopoda")
)
## IMPORTANT: Some species substitutions or removals were made as part of buildTree. Run build_tree_substitutions() to see them all.
## Pro tip: most tree read/write functions reset node numbers. Fortify your tree and save it as a csv file to preserve node numbering.

tree
## 
## Phylogenetic tree with 5 tips and 4 internal nodes.
## 
## Tip labels:
##   Amborella_trichopoda, Zea_mays, Sorghum_bicolor, Setaria_viridis, Arabidopsis_neglecta
## Node labels:
##   , , , 
## 
## Rooted; includes branch length(s).

plot(tree)
```

<img src="index_files/figure-html/unnamed-chunk-535-1.png" width="100%" style="display: block; margin: auto;" />

Note that `buildTree` informs us: "Scaffold newick tip Arabidopsis_thaliana substituted with Arabidopsis_neglecta". This means that *Arabidopsis neglecta* was grafted onto the tip originally occupied by *Arabidopsis thaliana*. This behaviour is useful when operating on a large phylogenetic scale (i.e. where *exact* phylogeny topology is not critical below the family level). However, if a person is interested in using an existing newick tree as a scaffold for a phylogeny where genus-level topology *is* critical, then beware! Your scaffold may not be appropriate if you see that message. When operating at the genus level, you probably want to use sequence data to build your phylogeny anyway. So let's look at how to do that:

### alignment input {-}

Arguments in this case are:

* "scaffold_type": "amin_alignment" or "nucl_alignment" for amino acids or nucleotides.
* "scaffold_in_path": path to the fasta file that contains the alignment from which you want to build a tree.
* "ml": Logical, TRUE if you want to use maximum liklihood, FALSE if not, in which case neighbor joining will ne used.
* "model_test": if you say TRUE to "ml", should buildTree test different maximum liklihood models and then use the "best" one? 
* "bootstrap": TRUE or FALSE, whether you want bootstrap values on the nodes.
* "ancestral_states": TRUE or FALSE, should buildTree() compute the ancestral sequence at each node?
* "root": NULL, or the name of an accession that should form the root of the tree.


``` r
buildTree(
  scaffold_type = "amin_alignment",
  scaffold_in_path = "/path_to/a_folder_for_alignments/all_amin_seqs.fa",
  ml = FALSE, 
  model_test = FALSE,
  bootstrap = FALSE,
  ancestral_states = FALSE,
  root = NULL
)
```

## plotting trees {-}

There are several approaches to plotting trees. A simple one is using the base `plot` function:


``` r
test_tree_small <- buildTree(
  scaffold_type = "newick",
  scaffold_in_path = "https://thebustalab.github.io/data/plant_phylogeny.newick",
  members = c("Sorghum_bicolor", "Zea_mays", "Setaria_viridis")
)
## Pro tip: most tree read/write functions reset node numbers. Fortify your tree and save it as a csv file to preserve node numbering.

plot(test_tree_small)
```

<img src="index_files/figure-html/unnamed-chunk-537-1.png" width="100%" style="display: block; margin: auto;" />

Though this can get messy when there are lots of tip labels:


``` r
set.seed(122)
test_tree_big <- buildTree(
  scaffold_type = "newick",
  scaffold_in_path = "https://thebustalab.github.io/data/plant_phylogeny.newick",
  members = plant_species$Genus_species[abs(floor(rnorm(60)*100000))]
)

plot(test_tree_big)
```

<img src="index_files/figure-html/unnamed-chunk-538-1.png" width="100%" style="display: block; margin: auto;" />

One solution is to use `ggtree`, which by default doesn't show tip labels. `plot` can do that too, but `ggtree` does a bunch of other useful things, so I recommend that:


``` r
ggtree(test_tree_big)
```

<img src="index_files/figure-html/unnamed-chunk-539-1.png" width="100%" style="display: block; margin: auto;" />

Another convenient fucntion is ggplot's `fortify`. This will convert your `phylo` object into a data frame:


``` r
test_tree_big_fortified <- fortify(test_tree_big)
test_tree_big_fortified
## # A tbl_tree abstraction: 101 × 9
## # which can be converted to treedata or phylo 
## # via as.treedata or as.phylo
##    parent  node branch.length label isTip     x     y branch
##     <int> <int>         <dbl> <chr> <lgl> <dbl> <dbl>  <dbl>
##  1     54     1          83.0 Wolf… TRUE   188.     1   147.
##  2     54     2          83.0 Spat… TRUE   188.     2   147.
##  3     55     3         138.  Dios… TRUE   188.     3   120.
##  4     58     4          42.7 Bulb… TRUE   188.     5   167.
##  5     58     5          42.7 Ober… TRUE   188.     6   167.
##  6     59     6          32.0 Poma… TRUE   188.     7   172.
##  7     59     7          32.0 Teli… TRUE   188.     8   172.
##  8     56     8         135.  Cala… TRUE   188.     4   121.
##  9     61     9         147.  Pepe… TRUE   188.     9   115.
## 10     62    10         121.  Endl… TRUE   188.    10   128.
## # ℹ 91 more rows
## # ℹ 1 more variable: angle <dbl>
```
`ggtree` can still plot this dataframe, and it allows metadata to be stored in a human readable format by using mutating joins (explained below). This metadata can be plotted with standard ggplot geoms, and these dataframes can also conveniently be saved as .csv files:


``` r

## Note that "plant_species" comes with the phylochemistry source.

test_tree_big_fortified_w_data <- left_join(test_tree_big_fortified, plant_species, by = c("label" = "Genus_species"))

test_tree_big_fortified_w_data
## # A tbl_tree abstraction: 101 × 14
## # which can be converted to treedata or phylo 
## # via as.treedata or as.phylo
##    parent  node branch.length label isTip     x     y branch
##     <int> <int>         <dbl> <chr> <lgl> <dbl> <dbl>  <dbl>
##  1     54     1          83.0 Wolf… TRUE   188.     1   147.
##  2     54     2          83.0 Spat… TRUE   188.     2   147.
##  3     55     3         138.  Dios… TRUE   188.     3   120.
##  4     58     4          42.7 Bulb… TRUE   188.     5   167.
##  5     58     5          42.7 Ober… TRUE   188.     6   167.
##  6     59     6          32.0 Poma… TRUE   188.     7   172.
##  7     59     7          32.0 Teli… TRUE   188.     8   172.
##  8     56     8         135.  Cala… TRUE   188.     4   121.
##  9     61     9         147.  Pepe… TRUE   188.     9   115.
## 10     62    10         121.  Endl… TRUE   188.    10   128.
## # ℹ 91 more rows
## # ℹ 6 more variables: angle <dbl>, Phylum <chr>,
## #   Order <chr>, Family <chr>, Genus <chr>, species <chr>

ggtree(test_tree_big_fortified_w_data) + 
  geom_point(
    data = filter(test_tree_big_fortified_w_data, isTip == TRUE),
    aes(x = x, y = y, fill = Order), size = 3, shape = 21, color = "black") +
  geom_text(
    data = filter(test_tree_big_fortified_w_data, isTip == TRUE),
    aes(x = x, y = y, label = y), size = 2, color = "white") +
  geom_tiplab(aes(label = label), offset = 10, size = 2) +
  theme_void() +
  scale_fill_manual(values = discrete_palette) +
  coord_cartesian(xlim = c(0,280)) +
  theme(
    legend.position = c(0.15, 0.75)
  )
```

<img src="index_files/figure-html/unnamed-chunk-541-1.png" width="100%" style="display: block; margin: auto;" />

## collapseTree {-}

Sometimes we want to view a tree at a higher level of taxonomical organization, or some other higher level. This can be done easily using the `collapseTree` function. It takes two arguments: an un-fortified tree (`tree`), and a two-column data frame (`associations`). In the first column of the data frame are all the tip labels of the tree, and in the second column are the higher level of organization to which each tip belongs. The function will prune the tree so that only one member of the higher level of organization is included in the output. For example, let's look at the tree from the previous section at the family level:


``` r
collapseTree(
  tree = test_tree_big,
  associations = data.frame(
    tip.label = test_tree_big$tip.label,
    family = plant_species$Family[match(test_tree_big$tip.label, plant_species$Genus_species)]
  )
) -> test_tree_big_families
## Branch lengths have been set to one.

ggtree(test_tree_big_families) + geom_tiplab() + coord_cartesian(xlim = c(0,300))
```

<img src="index_files/figure-html/unnamed-chunk-542-1.png" width="100%" style="display: block; margin: auto;" />

## trees and traits {-}

To plot traits alongside a tree, we can use ggtree in combination with ggplot. Here is an example. First, we make the tree:


``` r
chemical_bloom_tree <- buildTree(
  scaffold_type = "newick",
  scaffold_in_path = "http://thebustalab.github.io/data/angiosperms.newick",
  members = unique(chemical_blooms$label)
)
## IMPORTANT: Some species substitutions or removals were made as part of buildTree. Run build_tree_substitutions() to see them all.
## Pro tip: most tree read/write functions reset node numbers. Fortify your tree and save it as a csv file to preserve node numbering.
```

Next we join the tree with the data:

``` r
data <- left_join(fortify(chemical_bloom_tree), chemical_blooms)
## Joining with `by = join_by(label)`
head(data)
## # A tibble: 6 × 18
##   parent  node branch.length label  isTip     x     y branch
##    <int> <int>         <dbl> <chr>  <lgl> <dbl> <dbl>  <dbl>
## 1     80     1         290.  Ginkg… TRUE   352.     1   207.
## 2     81     2         267.  Picea… TRUE   352.     2   219.
## 3     81     3         267.  Cupre… TRUE   352.     3   219.
## 4     84     4         135.  Eryth… TRUE   352.     5   285.
## 5     86     5          16.2 Iris_… TRUE   352.     6   344.
## 6     86     6          16.2 Iris_… TRUE   352.     7   344.
## # ℹ 10 more variables: angle <dbl>, Alkanes <dbl>,
## #   Sec_Alcohols <dbl>, Others <dbl>, Fatty_acids <dbl>,
## #   Alcohols <dbl>, Triterpenoids <dbl>, Ketones <dbl>,
## #   Other_compounds <dbl>, Aldehydes <dbl>
```

Now we can plot the tree:

``` r
tree_plot <- ggtree(data) +
  geom_tiplab(
    align = TRUE, hjust = 1, offset = 350,
    geom = "label", label.size = 0, size = 3
  ) +
  scale_x_continuous(limits = c(0,750))
```

IMPORTANT! When we plot the traits, we need to reorder whatever is on the shared axis (in this case, the y axis) so that it matches the order of the tree. In this case, we need to reorder the species names so that they match the order of the tree. We can do this by using the `reorder` function, which takes two arguments: the thing to be reordered, and the thing to be reordered by. In this case, we want to reorder the species names by their y coordinate on the tree. We can do this by using the `y` column of the data frame that we created when we fortified the tree. We can then plot the traits:


``` r
trait_plot <- ggplot(
    data = pivot_longer(
      dplyr::filter(data, isTip == TRUE),
      cols = 10:18, names_to = "compound", values_to = "abundance"
    ),
    aes(x = compound, y = reorder(label, y), size = abundance)
  ) +
  geom_point() +
  scale_y_discrete(name = "") +
  theme(
    plot.margin = unit(c(1,1,1,1), "cm")
  )
```

Finally, we can plot the two plots together using `plot_grid`. It is important to manually inspect the tree tips and the y axis text to make sure that everything lines up. We don't want to be plotting the abundance of one species on the y axis of another species. In this case, everything looks good:


``` r
plot_grid(
  tree_plot,
  trait_plot,
  nrow = 1, align = "h", axis = "tb"
)
```

<img src="index_files/figure-html/unnamed-chunk-547-1.png" width="100%" style="display: block; margin: auto;" />


Once our manual inspection is complete, we can make a new version of the plot in which the y axis text is removed from the trait plot and we can reduce the margin on the left side of the trait plot to make it look nicer:


``` r
tree_plot <- ggtree(data) +
  geom_tiplab(
    align = TRUE, hjust = 1, offset = 350,
    geom = "label", label.size = 0, size = 3
  ) +
  scale_x_continuous(limits = c(0,750))

trait_plot <- ggplot(
    data = pivot_longer(
      filter(data, isTip == TRUE),
      cols = 10:18, names_to = "compound", values_to = "abundance"
    ),
    aes(x = compound, y = reorder(label, y), size = abundance)
  ) +
  geom_point() +
  scale_y_discrete(name = "") +
  theme(
    axis.text.y = element_blank(),
    plot.margin = unit(c(1,1,1,-1.5), "cm")
  )

plot_grid(
  tree_plot,
  trait_plot,
  nrow = 1, align = "h", axis = "tb"
)
```

<img src="index_files/figure-html/unnamed-chunk-548-1.png" width="100%" style="display: block; margin: auto;" />


# phylogenetic analyses {-}

We use a wrapper function to run phylogenetic comparative analyses. It 


``` r
chemical_bloom_tree <- buildTree(
  scaffold_type = "newick",
  scaffold_in_path = "http://thebustalab.github.io/data/angiosperms.newick",
  members = unique(chemical_blooms$label)
)
## IMPORTANT: Some species substitutions or removals were made as part of buildTree. Run build_tree_substitutions() to see them all.
## Pro tip: most tree read/write functions reset node numbers. Fortify your tree and save it as a csv file to preserve node numbering.

runPhylogeneticAnalyses(
    traits = pivot_longer(chemical_blooms[,1:4], cols = c(3:4), names_to = "trait", values_to = "value"),
    column_w_names_of_tiplabels = "label",
    column_w_names_of_traits = "trait",
    column_w_values_for_traits = "value",
    tree = chemical_bloom_tree
)
## Joining with `by = join_by(trait)`
## # A tibble: 310 × 17
##    parent  node branch.length label isTip     x     y branch
##     <int> <dbl>         <dbl> <chr> <lgl> <dbl> <dbl>  <dbl>
##  1     80     1         290.  Gink… TRUE   352.     1   207.
##  2     80     1         290.  Gink… TRUE   352.     1   207.
##  3     81     2         267.  Pice… TRUE   352.     2   219.
##  4     81     2         267.  Pice… TRUE   352.     2   219.
##  5     81     3         267.  Cupr… TRUE   352.     3   219.
##  6     81     3         267.  Cupr… TRUE   352.     3   219.
##  7     84     4         135.  Eryt… TRUE   352.     5   285.
##  8     84     4         135.  Eryt… TRUE   352.     5   285.
##  9     86     5          16.2 Iris… TRUE   352.     6   344.
## 10     86     5          16.2 Iris… TRUE   352.     6   344.
## # ℹ 300 more rows
## # ℹ 9 more variables: angle <dbl>, trait <chr>,
## #   value <dbl>, trait_type <chr>,
## #   phylogenetic_signal_k_value <dbl>,
## #   phylogenetic_signal_k_p_value <dbl>,
## #   phylogenetic_signal_lambda_value <dbl>,
## #   phylogenetic_signal_lambda_p_value <dbl>, pic <dbl>
```


For all the below, there are some structural requirements: (i) the tree needs to be a phylo object (ii) the traits need to be a data.frame in which each row is a species and each column is a variable, and (iii) the first column in the data.frame needs to be the names of the species and they must exactly match the tip labels of the tree (though they don't have to be in the same order), for example:


``` r
chemical_bloom_tree <- buildTree(
  scaffold_type = "newick",
  scaffold_in_path = "http://thebustalab.github.io/data/angiosperms.newick",
  members = unique(chemical_blooms$label)
)
## IMPORTANT: Some species substitutions or removals were made as part of buildTree. Run build_tree_substitutions() to see them all.
## Pro tip: most tree read/write functions reset node numbers. Fortify your tree and save it as a csv file to preserve node numbering.
```

## phylogeneticSignal {-}

Phylogenetic signal is a measure of the degree to which related species share similar trait values. It is used to determine whether a trait has evolved in a manner that is consistent with the species' evolutionary history. `phylochemistry` provides the `phylogeneticSignal` function, which can be used to calculate phylogenetic signal for a given set of traits and a phylogenetic tree. Here is an example:


``` r
phylogeneticSignal(
  traits = pivot_longer(chemical_blooms, cols = c(2:10), names_to = "compound", values_to = "value"),
  column_w_names_of_tiplabels = "label",
  column_w_names_of_traits = "compound",
  column_w_values_for_traits = "value",
  tree = chemical_bloom_tree
)
##             trait trait_type n_species number_of_levels
## 1         Alkanes continuous        78               NA
## 2    Sec_Alcohols continuous        78               NA
## 3          Others continuous        78               NA
## 4     Fatty_acids continuous        78               NA
## 5        Alcohols continuous        78               NA
## 6   Triterpenoids continuous        78               NA
## 7         Ketones continuous        78               NA
## 8 Other_compounds continuous        78               NA
## 9       Aldehydes continuous        78               NA
##   evolutionary_transitions_observed
## 1                                NA
## 2                                NA
## 3                                NA
## 4                                NA
## 5                                NA
## 6                                NA
## 7                                NA
## 8                                NA
## 9                                NA
##   median_evolutionary_transitions_in_randomization
## 1                                               NA
## 2                                               NA
## 3                                               NA
## 4                                               NA
## 5                                               NA
## 6                                               NA
## 7                                               NA
## 8                                               NA
## 9                                               NA
##   minimum_evolutionary_transitions_in_randomization
## 1                                                NA
## 2                                                NA
## 3                                                NA
## 4                                                NA
## 5                                                NA
## 6                                                NA
## 7                                                NA
## 8                                                NA
## 9                                                NA
##   evolutionary_transitions_in_randomization
## 1                                        NA
## 2                                        NA
## 3                                        NA
## 4                                        NA
## 5                                        NA
## 6                                        NA
## 7                                        NA
## 8                                        NA
## 9                                        NA
##   phylogenetic_signal_k_value phylogenetic_signal_k_p_value
## 1                 0.066016688                         0.167
## 2                 2.045611108                         0.001
## 3                 0.029719595                         0.711
## 4                 0.069056092                         0.292
## 5                 0.053761730                         0.354
## 6                 0.566806846                         0.001
## 7                 0.239488396                         0.030
## 8                 0.018420367                         0.868
## 9                 0.008613879                         0.928
##   phylogenetic_signal_lambda_value
## 1                           0.0001
## 2                           0.9999
## 3                           0.0001
## 4                           0.0001
## 5                           0.0001
## 6                           0.9999
## 7                           0.7874
## 8                           0.0001
## 9                           0.0001
##   phylogenetic_signal_lambda_p_value
## 1                              1.000
## 2                              0.000
## 3                              1.000
## 4                              1.000
## 5                              1.000
## 6                              0.000
## 7                              0.045
## 8                              1.000
## 9                              1.000
```

## independentContrasts {-}

Phylogenetic independent contrasts are a method for analyzing the relationship between two or more traits while taking into account the evolutionary history of the species being studied. This method involves transforming the data in to "independent contrasts" to remove the effects of shared ancestry, allowing for more accurate analysis of the relationship between traits. `phylochemistry` provides the `independentContrasts` function to calculate phylogenetic independent contrasts for a given set of traits and a phylogenetic tree. Here is an example of calculating independent contrasts for an example dataset, followed by generating a linear model based on the contrasts.


``` r
contrasts <- independentContrasts(
  traits = pivot_longer(chemical_blooms, cols = c(2:10), names_to = "compound", values_to = "value"),
  column_w_names_of_tiplabels = "label",
  column_w_names_of_traits = "compound",
  column_w_values_for_traits = "value",
  tree = chemical_bloom_tree
)

# buildLinearModel(
#   data = contrasts,
#   formula = "Fatty_acids = Alkanes + 0"
# ) -> model
# 
# ggplot(model$data) +
#   geom_point(aes(x = input_x, y = input_y)) +
#   geom_line(aes(x = model_x, model_y))
```

## ancestralTraits {-}

Ancestral trait reconstruction is a method to infer the characteristics (or "traits") of ancestral organisms based on the traits of their modern descendants. By examining the traits of present-day species and using phylogenetic trees, we can estimate or "reconstruct" the traits of common ancestors. This method can be applied to various types of traits, including continuously varying and discrete traits. Ancestral trait reconstruction helps us gain insights into the evolutionary processes and the historical transitions that led to current biodiversity. `phylochemistry` provides the function `ancestralTraits` to perform these operations. Note that `ancestralTraits` is different from `buildTree`s "ancestral_states". "ancestral_states"  estimates ancestral sequence states at phylogeny nodes, while `ancestralTraits` will estimate the traits of an ancestor, given the traits of extant species that are present on the leaves of a phylogeny. Here is an example. Please note that ancestralTraits accepts data in a long-style data frame.


``` r
anc_traits_tree <- ancestralTraits(
  traits = pivot_longer(chemical_blooms, cols = -1),
  column_w_names_of_tiplabels = "label",
  column_w_names_of_traits = "name",
  column_w_values_for_traits = "value",
  tree = chemical_bloom_tree
)
head(anc_traits_tree)
## # A tibble: 6 × 11
##   parent  node branch.length label  isTip     x     y branch
##    <int> <dbl>         <dbl> <chr>  <lgl> <dbl> <dbl>  <dbl>
## 1     80     1          290. Ginkg… TRUE   352.     1   207.
## 2     80     1          290. Ginkg… TRUE   352.     1   207.
## 3     80     1          290. Ginkg… TRUE   352.     1   207.
## 4     80     1          290. Ginkg… TRUE   352.     1   207.
## 5     80     1          290. Ginkg… TRUE   352.     1   207.
## 6     80     1          290. Ginkg… TRUE   352.     1   207.
## # ℹ 3 more variables: angle <dbl>, trait <chr>, value <dbl>
```

In addition to providing ancestral state estimations, there is also a function for plotting those estimations on a phylogeny: `geom_ancestral_pie`. Here is an example. Note that `cols` is a vector of column numbers that correspond to the traits of interest. `pie_size` is the size of the pie chart that will be plotted at each node. `geom_ancestral_pie` relies on having columns in its input called `trait` and `value`, such as those output by   `ancestralTraits`. Note that if you are passing an object to ggtree() that has duplicate node names, you will need to use the `distinct` function to remove the duplicates, otherwise geom_ancestral_pie will get confused about where to place the pies.


``` r
ggtree(
  distinct(anc_traits_tree, node, .keep_all = TRUE)
) +
  geom_ancestral_pie(
    data = filter(anc_traits_tree, isTip == FALSE),
    pie_size = 0.1, pie_alpha = 1
  ) +
  geom_tiplab(offset = 20, align = TRUE) +
  scale_x_continuous(limits = c(0,650)) +
  theme_void()
```

<img src="index_files/figure-html/unnamed-chunk-570-1.png" width="100%" style="display: block; margin: auto;" />

________________________________________________________________________________________________
________________________________________________________________________________________________
________________________________________________________________________________________________


# (PART) ANALYTICAL REPORTS {-}

<!-- start WRITING -->

<!-- # overview {-} -->

<!-- For your final project in this course you will use the techniques we have learned in class to analyze a large dataset, prepare high quality figures, and write a miniature manuscript describing the results: -->

<!-- ```{r fig.align='center', echo=FALSE, include=identical(knitr:::pandoc_to(), 'html'), results="markup"} -->
<!-- knitr:::include_graphics('https://thebustalab.github.io/integrated_bioanalytics/images/project_overview.png', dpi = NA) -->
<!-- ``` -->

<!-- 1. **Find a data set: Large!** >10ish variables, >5ish categories -->
<!-- + Sources: your research supervisor, CHEM5725 database spreadsheet, Google searches, Kaggle.com! -->
<!-- + Relevant to your research or interests (ideally). -->
<!-- + Requires approval from Dr. Busta. -->

<!-- 2. **Ask at least three scientific questions.** -->
<!-- + These should drive your data analyses. -->
<!-- + Requires approval from Dr. Busta. -->

<!-- 3. **Analyze your data using what you learned in class.** -->
<!-- + Refer to our book. -->
<!-- + Ask Dr. Busta for assistance. -->

<!-- 4. **Create a written overview of your analysis**. A mini-manuscript in R Markdown: -->
<!-- + *Content* similar to the articles we looked at in class, though shorter. -->
<!-- + *Layout* similar to this example ([pdf](https://github.com/thebustalab/thebustalab.github.io/blob/master/integrated_bioanalytics/final_project/final_project_example.pdf), [rmd](https://github.com/thebustalab/thebustalab.github.io/blob/master/integrated_bioanalytics/final_project/final_project_example.Rmd)). -->

<!-- ## scope {-} -->

<!-- When conducting a project of this type, it is very common for there to be mismatches in the scope of how the project was conducted and how the written report is presented (see image below). We often spend LOTS of time exploring our data and running into dead ends, conclusions that are mundane, or questions we can't answer. When we write a report on the project, we often focus the report of a specific discovery we made during our vast avenues of exploration, rather than boring the reader with all the mundane details. -->

<!-- ```{r fig.align='center', echo=FALSE, include=identical(knitr:::pandoc_to(), 'html'), results="markup"} -->
<!-- knitr:::include_graphics('https://thebustalab.github.io/integrated_bioanalytics/images/scope.jpeg', dpi = NA) -->
<!-- ``` -->

<!-- ## order and content {-} -->

<!-- The manuscript will be comprised of a title, abstract, introduction, results and discussion section, figures and captions, conclusions section, and at least five references. Please note the following when preparing your manuscript: the orders of presentation and preparation do not have to be the same (see the images below)! While in some instances a scientist may choose to write the components of a manuscript in the same order in which they appear on the page, this is not always the case. The order of preparation suggsted above is designed to minimize the amount of revision / re-writing that needs to be performed during the manuscript preparation process. Note that the suggested order of composition is in line with the class schedule for the rest of the semester. -->

<!-- ```{r fig.align='center', echo=FALSE, include=identical(knitr:::pandoc_to(), 'html'), results="markup"} -->
<!-- knitr:::include_graphics('https://thebustalab.github.io/integrated_bioanalytics/images/writing_order2.png', dpi = NA) -->
<!-- ``` -->

<!-- Here is a guide for documenting analysis in a format that is polished and comprehensible. Note that each section is written for a specific and slightly unique audience. Utilizing R Markdown, the document created will convey findings and narrate the research process. In crafting your R Markdown document, ensure that code remains concealed in the final presentation by including echo = FALSE within chunk headers. This action hides the R code blocks in the output, yet permits their execution for generating figures and results. -->

<!-- TITLE: Craft a title that’s both clear and descriptive. It should be accessible to specialists in the field as well as the wider scientific audience, avoiding overly technical language that might limit its broader appeal. -->

<!-- ABSTRACT: Summarize the study in the abstract, providing detail that will be informative to experts and also comprehensible to those outside the field. It should briefly outline the study’s aims, methods, key results, and conclusions, offering a snapshot of the entire project. -->

<!-- INTRODUCTION: In the introduction, present the context and significance of the research in a way that's understandable to researchers and those not as versed in the subject. Clearly enumerate the multiple research questions that the you aim to answer, highlighting the research's relevance and framing the inquiry. -->

<!-- METHODS: Describe where the dataset was sourced, any data processing steps taken, and provide a detailed description of the analysis procedures utilized in RStudio. This section should be detailed enough to allow for replication and validation, yet clearly written to be accessible to non-experts. The explanation of the methods should help readers understand how the research questions were addressed. -->

<!-- RESULTS AND DISCUSSION WITH FIGURES AND CAPTIONS: Integrate your findings with clear, illustrative figures and captions within this section, ensuring they can be understood independently of the text for visual learners. The written portion should concisely interpret the results in light of the research questions, discussing the significance of the findings in a way that appeals to both those interested in the analytical nuances and those seeking to understand the overall implications. -->

<!-- CONCLUSION: Summarize the primary insights and their relevance, keeping this section succinct and to the point. It should crystallize the main findings and their contribution to the field, suited for readers who seek a quick synopsis without delving into the full text. -->

<!-- REFERENCES: Include at least five references to substantiate your research, ensuring they are pertinent and formatted to facilitate easy follow-up for interested readers. The references should be organized to serve both as a trail for fellow researchers and a resource for those who are less experienced in the academic discourse. -->

<!-- ## scientific questions {-} -->

<!-- Scientific questions are pivotal in plant chemistry research, especially when examining the quantification of compounds in various plants, tissues, or environments. They orient the scope of inquiry and dictate the choice of analytical methods to draw pertinent conclusions from data. -->

<!-- DESCRIPTIVE questions might catalog the variety or concentration of phytochemicals present in a given species, often employing summary statistics to encapsulate the data: -->

<!-- - What is the average concentration of alkaloids found in the leaves of nightshade plants in temperate zones? -->

<!-- - How do the levels of flavonoids vary among different tissues of the grapevine? -->

<!-- - What is the frequency distribution of terpene profiles in pine populations across different altitudes? -->

<!-- CORRELATIVE questions investigate the relationships between environmental factors and chemical expression in plants, typically utilizing regression modeling: -->

<!-- - Does the level of UV radiation correlate with the production of protective anthocyanins in vineyard grape varieties? -->

<!-- - How is the accumulation of heavy metals in fern tissues related to soil pollution levels? -->

<!-- COMPARATIVE questions explore variations or consistencies across groups or conditions, often answered through statistical comparisons or pattern recognition methods like clustering: -->

<!-- - Which are more similar in their secondary metabolite profiles, medicinal herbs grown in greenhouse conditions or in the wild? -->

<!-- - What distinguishes the phenolic compound content in shade-grown coffee plants versus those grown in direct sunlight? -->

<!-- - Is there a significant difference in essential oil compositions between lavender plants cultivated in different soil types? -->

<!-- Unclear: How should social networking sites address the harm they cause?
Clear: What action should social networking sites like MySpace and Facebook take to protect users’ personal information and privacy?

The unclear version of this question doesn’t specify which social networking sites or suggest what kind of harm the sites might be causing. It also assumes that this “harm” is proven and/or accepted. The clearer version specifies sites (MySpace and Facebook), the type of potential harm (privacy issues), and who may be experiencing that harm (users). A strong research question should never leave room for ambiguity or interpretation.

Unfocused: What is the effect on the environment from global warming?
Focused: What is the most significant effect of glacial melting on the lives of penguins in Antarctica?

The unfocused research question is so broad that it couldn’t be adequately answered in a book-length piece, let alone a standard college-level paper. The focused version narrows down to a specific effect of global warming (glacial melting), a specific place (Antarctica), and a specific animal that is affected (penguins). It also requires the writer to take a stance on which effect has the greatest impact on the affected animal. When in doubt, make a research question as narrow and focused as possible.

Too simple: How are doctors addressing diabetes in the U.S.?
Appropriately Complex:  What main environmental, behavioral, and genetic factors predict whether Americans will develop diabetes, and how can these commonalities be used to aid the medical community in prevention of the disease? -->


# figures & captions {-}

## {-}

## figures {-}

One of the first components in preparing a scientific manuscript is creating high quality figures. Considering the following for your figures:

- General Appearance:

Create plots that are clean, professional, and easy to view from a distance. Ensure axes tick labels are clear, non-overlapping, and utilize the available space efficiently for enhanced readability and precision. Use an appealing (and color blind-friendly) color palette to differentiate data points or categories. Tailor axes labels to be descriptive, and select an appropriate theme that complements the data and maintains professionalism.

- Representing Data:

Appropriate Geoms and Annotations: Choose geoms that best represent the data and help the viewer evaluate the hypothesis or make the desired comparison. Include raw data points where possible for detailed data distribution understanding. Consider apply statistical transformations like smoothing lines or histograms where appropriate to provide deeper insights into the data. Consider using facets for visualizing multiple categories or groups, allowing for easier comparison while maintaining a consistent scale and layout. Adhere to specific standards or conventions relevant to your field, including the representation of data, error bars, or statistical significance markers.

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/plot_quality.jpg" width="100%" style="display: block; margin: auto;" />

### advanced figure elements {-}

### insets {-}

- zoomed insets

Zoom in on certain plot regions


``` r
p <- ggplot(mpg, aes(displ, hwy, colour = factor(cyl))) +
  geom_point() 

data.tb <- 
  tibble(x = 7, y = 44, 
         plot = list(p + 
                       coord_cartesian(xlim = c(4.9, 6.2), 
                                       ylim = c(13, 21)) +
                       labs(x = NULL, y = NULL) +
                       theme_bw(8) +
                       scale_colour_discrete(guide = "none")))

ggplot(mpg, aes(displ, hwy, colour = factor(cyl))) +
  geom_plot(data = data.tb, aes(x, y, label = plot)) +
  annotate(geom = "rect", 
           xmin = 4.9, xmax = 6.2, ymin = 13, ymax = 21,
           linetype = "dotted", fill = NA, colour = "black") +
  geom_point() 
```

<img src="index_files/figure-html/unnamed-chunk-21-1.png" width="100%" style="display: block; margin: auto;" />

- plot insets


``` r
p <- ggplot(mpg, aes(factor(cyl), hwy, fill = factor(cyl))) +
  stat_summary(geom = "col", fun = mean, width = 2/3) +
  labs(x = "Number of cylinders", y = NULL, title = "Means") +
  scale_fill_discrete(guide = "none")

data.tb <- tibble(x = 7, y = 44, 
                  plot = list(p +
                                theme_bw(8)))

ggplot(mpg, aes(displ, hwy, colour = factor(cyl))) +
  geom_plot(data = data.tb, aes(x, y, label = plot)) +
  geom_point() +
  labs(x = "Engine displacement (l)", y = "Fuel use efficiency (MPG)",
       colour = "Engine cylinders\n(number)") +
  theme_bw()
```

<img src="index_files/figure-html/unnamed-chunk-22-1.png" width="100%" style="display: block; margin: auto;" />

- image insets


``` r
Isoquercitin_synthase <- magick::image_read("https://thebustalab.github.io/integrated_bioanalytics/images/homology2.png")
grobs.tb <- tibble(x = c(0, 10, 20, 40), y = c(4, 5, 6, 9),
                   width = c(0.05, 0.05, 0.01, 1),
                   height =  c(0.05, 0.05, 0.01, 0.3),
                   grob = list(grid::circleGrob(), 
                               grid::rectGrob(), 
                               grid::textGrob("I am a Grob"),
                               grid::rasterGrob(image = Isoquercitin_synthase)))

ggplot() +
  geom_grob(data = grobs.tb, 
            aes(x, y, label = grob, vp.width = width, vp.height = height),
            hjust = 0.7, vjust = 0.55) +
  scale_y_continuous(expand = expansion(mult = 0.3, add = 0)) +
  scale_x_continuous(expand = expansion(mult = 0.2, add = 0)) +
  theme_bw(12)
```

<img src="index_files/figure-html/unnamed-chunk-23-1.png" width="100%" style="display: block; margin: auto;" />


``` r
# ggplot() +
#   annotate("grob", x = 1, y = 3, vp.width = 0.5,
#            label = grid::rasterGrob(image = Isoquercitin_synthase, width = 1)) +
#   theme_bw(12)
```



``` r
# bloom_example_pics <- ggplot(data = data.frame(x = c(0,1), y = c(0.5,0.5))) +
#   geom_point(aes(x = x, y = y), color = "white") +
#   theme_void() +
#   annotation_custom(
#       rasterGrob(
#           png::readPNG(
#               "https://thebustalab.github.io/integrated_bioanalytics/images/homology2.png"
#           ), interpolate=TRUE
#       ), xmin=0, xmax=1, ymin=0, ymax=1
#   )

```

### composite figures {-}

Many high quality figures are composite figures in which there is more than one panel. Here is a simple way to make such figures in R. First, make each component of the composite figure and send the plot to a new object:


``` r
color_palette <- RColorBrewer::brewer.pal(11, "Paired")
names(color_palette) <- unique(alaska_lake_data$element)

plot1 <- ggplot(
  data = filter(alaska_lake_data, element_type == "bound"),
  aes(y = lake, x = mg_per_L)
) +
  geom_col(
    aes(fill = element), size = 0.5, position = "dodge",
    color = "black"
  ) +
  facet_grid(park~., scales = "free", space = "free") +
  theme_bw() + 
  scale_fill_manual(values = color_palette) +
  scale_y_discrete(name = "Lake Name") +
  scale_x_continuous(name = "Abundance mg/L)") +
  theme(
    text = element_text(size = 14)
  )

plot2 <- ggplot(
  data = filter(alaska_lake_data, element_type == "free"),
  aes(y = lake, x = mg_per_L)
) +
  geom_col(
    aes(fill = element), size = 0.5, position = "dodge",
    color = "black"
  ) +
  facet_grid(park~., scales = "free", space = "free") +
  theme_bw() + 
  scale_fill_manual(values = color_palette) +
  scale_y_discrete(name = "Lake Name") +
  scale_x_continuous(name = "Abundance mg/L)") +
  theme(
    text = element_text(size = 14)
  )
```

Now, add them together to lay them out. Let's look at various ways to lay this out:


``` r
plot_grid(plot1, plot2)
```

<img src="index_files/figure-html/unnamed-chunk-27-1.png" width="100%" style="display: block; margin: auto;" />


``` r
plot_grid(plot1, plot2, ncol = 1)
```

<img src="index_files/figure-html/unnamed-chunk-28-1.png" width="100%" style="display: block; margin: auto;" />


``` r
plot_grid(plot_grid(plot1,plot2), plot1, ncol = 1)
```

<img src="index_files/figure-html/unnamed-chunk-29-1.png" width="100%" style="display: block; margin: auto;" />

### exporting graphics {-}

To export graphics from R, consider the code below. The <path_to_file_you_want_to_create> should be something like: "C:\\Desktop\\the_file.png" (i.e. a path to a specific file with a .png suffix. It should be a file that does not yet exist - if it does already exist, it will be overwritten. You should adjust with height and width to get the image to look how you want, then once you have that dialed in, crank the resolution to 1200 or 2400 and export a final version.


``` r
plot <- ggplot(data, aes(x = x, y = y)) + geom_point()

png(filename = <path_to_file_you_want_to_create>, width = 8, height = 8, res = 600, units = "in")

plot

dev.off()
```


``` r
plot <- ggplot(data, aes(x = x, y = y)) + geom_point()

pdf(filename = <path_to_file_you_want_to_create>, width = 8, height = 8)

plot

dev.off()
```

## captions {-}

Figures are critical tools for clearly and effectively communicating scientific results. However, as Reviewer 2 will tell you, a figure is only as good as its caption. Captions provide essential context, guiding the reader through the significance, structure, and details of the visual information presented. Below are some guidelines to help you craft informative captions. The recommendations are organized into categories, covering essential components like figure titles, panel descriptions, variable definitions, data representation details, statistical analyses, and data sources. Some example captions and a helpful interactive tool (`buildCaption()`) are also included to streamline caption construction and ensure consistency in your scientific communication.

### title and text {-}

- **Figure Title:**  
- Provide a concise, descriptive title that summarizes the overall message or purpose of the figure.
- Ensure the title quickly informs the reader about the main topic, experimental system, or hypothesis addressed by the figure.

- **In All Caption Text:**
- Avoid using unexplained abbreviations or jargon. If abbreviations are necessary, provide definitions at first use.

### panel-by-panel descriptions {-}

- **Panel Identification:**  
  - Label each panel (e.g., A, B, C, etc.) and refer to these labels consistently in the caption.

- **Graph Type and Layout:**  
  - State the type of plot (line plot, bar chart, scatter plot, histogram, etc.). When in doubt, "plot" is okay.
  - Describe any special features such as insets, overlays, or embedded plots (e.g., zoomed regions, additional mini-panels).

- **Axes and Variables:**  
  - Clearly define what is on each axis (x vs. y) in descriptive terms, including units of measurement (e.g., time in seconds, concentration in µM).
  - Explain if additional dimensions (such as color coding, marker sizes, or symbols) are used to represent extra variables.

- **Data Representation Details:**
  - Describe what the individual data points, bars, or error bars represent. For example:
    - **Data Points/Bars:** Explain whether they indicate individual measurements, means, medians, or other summary statistics.
    - **Error Bars:** Specify whether these indicate standard error, standard deviation, 95% confidence intervals, or another metric.
  - Note any graphical elements like trend lines or regression lines and what model or fit has been applied.

- **Sample Size and Replicates:**  
  - Indicate the number of independent samples or experimental replicates underlying each element of the graph.

- **Statistical Analysis and Comparisons:**
  - Describe any control experiments or baseline data presented in the figure, including how they were used to validate or compare with experimental results.
  - State any statistical tests used (e.g., t-test, ANOVA, regression analysis) and the significance level(s).
  - Describe how statistical significance is indicated in the figure (e.g., asterisks, brackets, p-value annotations).
  - Provide any necessary details about data normalization, transformation, or curve fitting that influence data interpretation.
  - For figures with scale bars or reference markers (e.g., microscopy images), specify the scale explicitly.
  
### data source and methodology {-}

- **Data Origins and Methodology:**  
  - Clearly state where the data come from (e.g., experimental assays, clinical samples, simulations, or databases).
  - If the data are derived from previously published work or a public repository, include proper references or accession numbers, if reasonable / possible.
  - Consider including a summary of the methods used to obtain or generate the data. This is standard practice in some fields - have a look at what articles from your field typically include.
  - Consider including information on experimental conditions (e.g., treatment concentrations, temperature, environmental conditions) or computational parameters (e.g., algorithm settings), assuming these weren't already mentioned when describing the axes.
  - Mention any image processing steps (e.g., brightness/contrast adjustments, background subtraction) if these steps are critical for understanding the visual data.

### example captions {-}

<div class="figure" style="text-align: center">
<img src="index_files/figure-html/unnamed-chunk-32-1.png" alt="Figure 1: Carbon, nitrogen, and phosphorous in Alaskan lakes. A) A bar chart showing the abundance (in mg per L, x-axis) of the bound elements (C, N, and P) in various Alaskan lakes (lake names on y-axis) that are located in one of three parks in Alaska (park names on right y groupings). B) A bar chart showing the abundance (in mg per L, x-axis) of the free elements (Cl, S, F, Br, Na, K, Ca, and Mg) in various Alaskan lakes (lake names on y-axis) that are located in one of three parks in Alaska (park names on right y groupings). The data are from a public chemistry data repository. Each bar represents the result of a single measurement of a single analyte, the identity of which is coded using color as shown in the color legend. Abbreviations: BELA - Bering Land Bridge National Preserve, GAAR - Gates Of The Arctic National Park &amp; Preserve, NOAT - Noatak National Preserve." width="100%" />
<p class="caption">(\#fig:unnamed-chunk-32)Figure 1: Carbon, nitrogen, and phosphorous in Alaskan lakes. A) A bar chart showing the abundance (in mg per L, x-axis) of the bound elements (C, N, and P) in various Alaskan lakes (lake names on y-axis) that are located in one of three parks in Alaska (park names on right y groupings). B) A bar chart showing the abundance (in mg per L, x-axis) of the free elements (Cl, S, F, Br, Na, K, Ca, and Mg) in various Alaskan lakes (lake names on y-axis) that are located in one of three parks in Alaska (park names on right y groupings). The data are from a public chemistry data repository. Each bar represents the result of a single measurement of a single analyte, the identity of which is coded using color as shown in the color legend. Abbreviations: BELA - Bering Land Bridge National Preserve, GAAR - Gates Of The Arctic National Park & Preserve, NOAT - Noatak National Preserve.</p>
</div>

### buildCaption {-}

To help you manage the suggestions above, please consider using the `buildCaption()` tool, which you can open using the command below. That command should open an interactive window with a checklist to help you quickly build quality captions.


``` r
buildCaption()
```

## {-}

## further reading {-}

- [Grammar extensions and insets with `ggpp`](https://docs.r4photobiology.info/ggpp/articles/grammar-extensions.html#geom_plot)  
  This article explains how to use the `ggpp` extension to add insets and annotations to `ggplot2` graphics in R. It introduces grammar extensions that allow you to insert subplots, highlight specific regions, and incorporate custom graphical elements in a composable and expressive way. Particularly useful for emphasizing detail or providing context within complex figures.

- [Patchwork: Simple plot layout with ggplot2](https://patchwork.data-imaginist.com/index.html)  
  `Patchwork` is an elegant and intuitive package for arranging multiple `ggplot2` plots into a single composite figure. With a minimal syntax that mirrors mathematical layout expressions, it allows users to combine plots vertically, horizontally, or in nested arrangements—ideal for creating figure panels for publications or presentations.

- [Cowplot: Versatile plot composition](https://wilkelab.org/cowplot/articles/plot_grid.html)  
  `Cowplot` is another popular package for composing multiple `ggplot2` plots. It offers more control and customization than `patchwork`, particularly for aligning plots, adjusting spacing, and embedding annotations. This makes it well-suited for fine-tuned figure design when preparing publication-quality graphics.


# question-driven report {-}

<!-- * Objective: Walk your reader through your results, drawing conclusions as you go. -->

<!-- * Tense: past tense and passive voice, because we are talking about what happened in the experiments, and reflecting with distance on the results.
 -->
<!-- * As you go: make notes of what should go into the introduction. -->

- A structured report with the following components:
  + Cover page: Title and author, Introductory paragraph
  + One section for each of three scientific question. Each section has one figure + caption, as well as one or more paragraphs explaining that figure and any conclusions that follow from it.
  + References

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/report_outline.jpg" width="100%" style="display: block; margin: auto;" />

## {-}

## brief structural example {-}

> **Chemical Pollutants in Minnesota Soils**\
>
>
>In order to better understand pollution in the state of Minnesota, this study focused on a detailed analyses of chemical measurements from soil samples from 300 sites around the state. The analyses consisted of a principal components analysis to determine which sites were similar to one another, thus answering the question of what sites exhibited similar pollutant profiles (Section 2.1). This first analysis was followed by statistical tests to see whether differences could be detected in the sites' chemistry, answering the question of whether there were any significant differences in the pollutant profiles at each site (Section 2.2).\
>
>
>*2.1 Principal Components Analysis*\
>
>
>To understand relationships between the sites from which soil chemistry was sampled, a principal components analysis was used. Each of the 20 different analytes, all of which contained halogen atoms, were included in the analysis. A scatter plot showing the position of each of the 300 samples in a space defined by dimesions 1 and 2 (which explain 54% and 35% of the total variance of the dataset, respectively), revealed that two major clusters are present, with a small number of outliers (Fig. 1). By color coding these two clusters according to whether the samples were collected from rural versus urban areas, it was possible to see that the first cluster was made out of almost exclusively samples from urban areas, while the second cluster was made up of almost entirely samples from rural areas. This suggested that variance in pollutant chemistry among the samples collected was assocaited with urban versus rural environments.\
>
>
>*2.2 Statistical Analyses*\
>
>
>Using the groupings that were identified via principal components analysis, statistical tests were conducted to determine if chemical abundances differed between groups. Tests for normality and homogeneity of variance (Shapiro and Levene tests) revealed that the data could not be assessed using ANOVA but instead required the use of a non-parametric test. Accordingly, the Kruskall-Wallis test followed by post-hoc Dunn tests were applied, which showed that the abundances of halogenated pollutants is significantly higher in urban versus rural areas (p = 0.0035, Fig. 2A). These direct observations are consistent with conclusions drawn by others in recent literature reviews focused on hydrocarbon compounds (Petrucci et al., 2018; Hendrix et al., 2019). *Thus, the new chemical analyses presented here demonstrate that the discrepancy in urban versus rural pollution is true not only for hydrocarbon compounds (as had been found previously), but also for halogenated compounds.* Together, these findings strongly suggest that either cities are a source of more pollution or that there is some other mechanism that concentrates pollution in cities.\
>

## structure {-}

(key: **number of suggested sentences**: *purpose*: "example")

* Title and Author
  + **1** Use about 75-140 characters (ideally no more than 125 characters). There are essentially two types of titles: descriptive titles and mechanistic titles. 
    + If your manuscript is exploratory research, consider using a descriptive title. For example: "Comparative analysis of carbon, sulfur, and phoshorous chemistry in six Alaskan lakes."
    + If your manuscript is hypothesis-driven research, consider using a mechanistic title. For example: "Dissolved organic carbon in Alaskan lakes is heavily influenced by water pH and temperature."
    + A good title should:
      + Be indicative of the content of the paper
      + Attract the interest of potential readers
      + Reflect whether the article is deascriptive or mechanistic
      + Include important keywords

* Introductory paragraph:
  + **1**: *State the aim of the report*: "The objective of this report is to..."
  + **3-4**: *Call out the subsections of the report according to methodology and scientific questions*: "We used method X to quantify property Y of our study subject and address the question "Does property Y vary based on the date the sample was collected?" (section 2.1)."

* Each subsection paragraph:
  <!-- + 1: *Statement of observation or lead-in*: "Based on the observation of X..." -->
  + **1**: *Purpose of the work described in this paragraph*: "In order to determine..."
  + **1**: *Review methods or experimental design specific to this subsection (if necessary)*
  + **4-5**: *Results of that method or experiment (i.e. data features)*
  + **1-2**: *Comparison of new results against those in literature (if possible)*
  + **1-2**: *Conclusion from the combined results or some other concluding remark* "Thus, analysis X revealed that..."
  <!-- + **1**: *Interpretation of the conclusion in a larger context (if possible / reasonable)* -->

<!-- * Introductory paragraph: -->
<!--   + **1**: *Review the aim of the paper*: "In order to understand…" -->
<!--   + **3-4**: *Use a methods summary to call out subsections*: "We used method X to quantify property Y of our study subject (section 2.1)" -->

<!-- * Each subsection paragraph: -->
<!--   <!-- + 1: *Statement of observation or lead-in*: "Based on the observation of X..." -->
<!--   + **1**: *Purpose of the work described in this paragraph*: "In order to determine..." -->
<!--   + **1**: *Review methods or experimental design specific to this subsection (if necessary)* -->
<!--   + **4-5**: *Results of that method or experiment (i.e. data features)* -->
<!--   + **1-2**: *Comparison of new results against those in literature (if possible)* -->
<!--   + **1-2**: *Conclusion from the combined results or some other concluding remark* "Thus, analysis X revealed that..." -->
<!--   <!-- + **1**: *Interpretation of the conclusion in a larger context (if possible / reasonable)* -->

## suggestions {-}

Is there an efficient way to write in the format outlined above? Yes. Follow the step-by-step instructions below:

### outline then draft paragraphs {-}

<img src="https://thebustalab.github.io/integrated_bioanalytics/images/res_disc_1.jpg" width="100%" style="display: block; margin: auto;" />

1. **Identify "data features" -> "conclusion" combinations.** Using your figures, make a list of all the potentially interesting features in your data, then pair each feature with a possible conclusiona it could lead to. Example:
  + "The GC-MS data presented here indicates that cities have higher levels of pollution than rural areas (Fig. 1)," (a data feature)
  + "suggesting that either cities are a source of more pollution or that there is some other mechanism that concentrates pollution in cities." (a conclusion)
  
2. **Perform targeted literature searches.** Expand your "data feature" -> "conclusion" combinations with "supplementary information" or "literature information". Example:
  + "The GC-MS data presented here indicates that cities have higher levels of pollution than rural areas." (data feature)
  + "These direct observations are consistent with conclusions drawn by others in recent literature reviews (Petrucci., 2018; Hendrix et al., 2019)" (literature information)
  + "Overall, this suggests that either cities are a source of more pollution or that there is some other mechanism that concentrates pollution in cities." (conclusion)
  
3. **Group "data feature" -> "supp/lit info" -> "conclusion" combinations into paragraphs.** Edit each conclusion so that it highlights what new contribution your data makes to the situation. Also consider whether any of the parargaphs now suggest the existence of mechanisms. Example (note the conclusion sentence in italics that highlights the new findings):
  + "The GC-MS data presented here indicates that cities have higher levels of pollution than rural areas (Fig. 1). These direct observations are consistent with meta-analyses of previously published observations (Supplemental Figure 1), as well as with conclusions drawn by others in recent literature reviews (So and so et al., 2018; The other person et al., 2019). *The new chemical analyses presented here thus confirm this is true for hydrocarbon compounds, and extend the observation to halogenated compounds in the atmosphere.* Together these findings strongly suggest that either cities are a source of more pollution or that there is some other mechanism that concentrates pollution in cities.

### order then edit paragraphs {-}

1. **Identify paragraph characteristics and group**
  + Consider whether any of your paragraphs are prerequisites for others and whether any paragraphs can be grouped according to topic.
  + Group paragraphs according to topic and prerequisite dependencies (putting prereq dependencies as close to eachother as possible.)

2. **Rearrange paragraph groups** Create the most natural flow. Consider:
  + Starting with group of paragraphs most relevant to the overall pitch/goal of the paper
  + Ending on the group of paragraphs that has the most future perspective
  + Ending in a strong suit (i.e. not something too speculative)
  + Consider putting orphaned paragraphs (or a shortened version of them) into the conclusion section.

3. **Edit transitions between groups.** Edit each paragraph, particularly its first and last sentences, to connect the paragraphs into a flowing document. Specifically, this means several things:
  + There should be no implicit cross-paragraph references (i.e. a new paragraph should not begin "The compound described above exhibited other interesting properties", rather, "3-hydroxycinnamic acid exhibited other interesting properties.").
  + There should be no abrupt jumps in subject between paragraphs, if there are consider breaking the discussion into subsections to help the reader identify logical resting points.
  + The discussion should not require the reader to go back and read its first half in order to understand its second half.

<!-- ### other thoughts -->

<!-- * Somewhere in the discussion, be sure to list out what any unsolved problems you faced are. -->

<!-- end -->

<!-- start Conclusions -->

<!-- # conclusion and introduction {-} -->

<!-- * Objective (conclusion): to convey a short statement of the take-home messages of your study. What are the most important things that you want the reader to remember from your study? -->

<!-- * Objective (introduciton): to prepare the reader by giving the reader sufficient background to understand the study as a whole. It therefore should only contain information pertinent to understanding the study and its broader significance.  -->

<!-- * Make sure that the scope of your introduction is in-line with the scope of the conclusion. That way, the reader will not be underwhelmed, nor will your work be undersold. -->

<!-- ## structure {-} -->

<!-- **Conclusion:** -->

<!-- * *One paragraph* -->
<!--   <!-- first sentence of the conclusion should be a restatment of the goal -->
<!--   + **2-3**: Summarize over-arching conclusions from each section of the paper (omit the details described in results or discussion) -->
<!--   + **2-3**: Based on a general description of findings, use pros and cons to argue for, if possible, alternative hypotheses. -->
<!--   + **1-2**: Suggest experiments to test these hypotheses. -->
<!--   + **1-2**: Describe future directions. -->
<!-- <!-- * Integrate over-arching conclusions to discuss new avenues, i.e. integrating chain length specificity and secondary functional group installation to discuss how both might affect physical properties of wax mixtures, then how these might have evolved. -->

<!-- **Introduction:** -->

<!-- ```{r fig.align='center', echo=FALSE, include=identical(knitr:::pandoc_to(), 'html'), results="markup"} -->
<!-- knitr:::include_graphics('https://thebustalab.github.io/integrated_bioanalytics/images/knowledge_gaps.jpeg', dpi = NA) -->
<!-- ``` -->

<!-- * *Paragraph 1: Introduce the topic* -->
<!--   + **1**: Introduce a topic and, ideally, an application of the research you will describe. Grab reader's attention. -->
<!--   + **1**: State why the topic is important. -->
<!--   + **1**: Describe what is known about the topic (at least, as pertains to the work at hand). -->
<!--   + **1**: Identify a gap in knowledge: "despite research in this area, here is what we don't know about the topic." -->
<!--   + **1**: List the negative things that will happen if we don't fill this gap in knowledge. -->

<!-- * *Paragraph 2: Provide background information* -->
<!--   + **3-5**: Describe, in moderate detail, the background information (concepts, literature) relevant to the study. -->
<!--   + **1**: End by saying how the details you just described relate to the application/topic described in the first paragraph. -->

<!-- * *Paragraph 3: Objectives of this study* -->
<!--   + **1**: State the objective of this study. -->
<!--   + **1**: Briefly describe what was done and the techniques or instruments used. -->
<!--   + **1-2**: For this project, briefly describe where you got the data, how you cleaned it up, if you merged multiple datasets, etc. -->
<!--   + **1**: (optional) State the major conclusion from the work and what it means for the application described in paragraph 1. -->

<!-- ## suggestions {-} -->

<!--   + If something is well-established, say so. -->
<!--   + Be clear about what is speculation. -->
<!--   + Last paragraph can mention objectives in list form. -->
<!--   + Last sentence can briefly mention methods (specific techniques or instruments) that were used. -->

<!-- end -->

<!-- start Abstract -->

<!-- # abstract and title {-} -->

<!-- ## abstract -->

<!-- ```{r fig.align='center', echo=FALSE, include=identical(knitr:::pandoc_to(), 'html'), results="markup"} -->
<!-- knitr:::include_graphics('https://thebustalab.github.io/integrated_bioanalytics/images/abstract_guide.jpeg', dpi = NA) -->
<!-- ``` -->

<!-- * **Structure** *One paragraph* Use about 200 - 500 words (ideally no more than 400 words) -->
<!--   + **1-2 sentences**: Introduction: Describe the topic, the motivation, and overall purpose of the research (Why is this research interesting and important? What gap in our knowledge does it fill?) -->
<!--   + **1-2 sentences**: Objective: Specific research objective, and potentially hypotheses/predictions, if any. -->
<!--   + **1-2 sentences**: Methods: Very concise overview of the methods used to address the research questions. -->
<!--   + **2-3 sentences**: Results/Discussion: Describe the major results (what you found) and interpretation of the results (what the results mean). -->
<!--   + **1-2 sentences**: Conclusions: Synthesizes the major contributions of the study into the context of the larger field to which the study belongs. What did we learn about the bigger picture of this field in general from doing this study? -->

<!-- * **Function: an abstract proves a short summary of the entire study.** The abstract should include the motivation or reason for conducting the study, what the research question or hypothesis was, how the experiments were conducted, what the results were, how the results are interpreted in light of the research question or hypothesis, and a concluding sentence about the general contribution or importance of the study. A good abstract should: -->
<!--   + Inform readers about the article’s content -->
<!--   + Summarize complex information in a clear, concise manner -->
<!--   + Help readers decide whether or not to read the article -->
<!--   + Used in conferences to summarize what the speaker will say during his/her presentation -->

<!-- ### other thoughts -->

<!-- * Whatever you identify as the strongest sentences in the main text, make sure those are reflected in the abstract -->

<!-- ## further reading {-} -->

<!-- * [Titles Guide](https://libguides.usc.edu/writingguide/title) -->

<!-- * [Abstract Guide] (https://www.cbs.umn.edu/sites/default/files/public/downloads/Annotated_Nature_abstract.pdf) -->

<!-- end -->
________________________________________________________________________________________________
________________________________________________________________________________________________
________________________________________________________________________________________________

# (PART) APPENDIX {-}

<!-- start links -->

# links {-}

## geoms {-}

[geoms and ggplot2 cheatsheet](https://thebustalab.github.io/integrated_bioanalytics/images/ggplot2_geoms.pdf)

## colors {-}

[ColorBrewer2](https://colorbrewer2.org/)

<!-- end -->

# datasets {-}

## alaska_lake_data {-}

The Alaska Lake Data was collected as part of a water quality monitoring initiative across various lakes in protected national parks, aimed at assessing the chemical composition and environmental conditions of these unique ecosystems. Researchers took water samples from several lakes, measuring key environmental parameters like water temperature and pH, along with analyzing the abundance of different chemical elements, such as carbon, nitrogen, and phosphorus. By comparing the concentrations of both bound and free elements, the study aimed to understand the health of these aquatic environments and the impact of natural and anthropogenic factors on the water chemistry. The dataset will be used to inform conservation strategies for maintaining the ecological balance in these sensitive regions.
- lake (Categorical): The name of the lake from which the water sample was collected. This refers to the sample.
- park (Categorical): The park or national park code where the lake is located. This is part of the sample identification.
water_temp (Continuous): The water temperature (in degrees Celsius) of the sample at the time of collection. This is an analyte describing an environmental condition of the sample.
- pH (Continuous): The pH value of the water, representing its acidity or alkalinity. This is an analyte providing an environmental characteristic of the sample.
- element (Categorical): The chemical element being measured in the water (e.g., C for carbon, N for nitrogen). This is an analyte.
- mg_per_L (Continuous): The concentration (in milligrams per liter) of the corresponding analyte from the element column, indicating the abundance of each analyte in the water sample.
- element_type (Categorical): Describes whether the element is in a "bound" or "free" state, providing context for the form of the analyte.

## algae_data {-}

This dataset was generated as part of a study investigating the biochemical composition of different algae strains under varying harvesting conditions. The goal of the research was to examine how different algae strains and harvesting regimes affect the abundance of various chemical species, particularly fatty acids and amino acids, which have potential applications in biofuel production and nutritional supplements. Replicates were performed to ensure consistency, and a wide range of chemical species was measured to provide insights into the algae's metabolic profile and its response to environmental or harvesting changes.

- replicate (Categorical): The replicate number of the sample for the experiment, indicating which iteration of the algae sample was analyzed.
- algae_strain (Categorical): The specific strain of algae used in the experiment (e.g., "Tsv1"). This refers to the strain from which each sample was collected.
- harvesting_regime (Categorical): The method or condition under which the algae sample was harvested (e.g., "Heavy" regime).
- chemical_species (Categorical): The type of chemical species or analyte measured in the algae sample, including various fatty acids (FAs) and amino acids (Aas).
- abundance (Continuous): The measured abundance of the chemical species or analyte in the algae sample, expressed in a continuous quantitative form (e.g., mg/L or similar units).

## beer_components {-}

This dataset captures the volatile compounds released from different ingredients like barley and corn, likely as part of a study on food aroma profiles. Researchers measured the abundance of specific analytes (such as 2-Methylpropanal) and classified them by chemical group (e.g., Aldehydes). The goal is to assess how different ingredients contribute to the overall aroma by linking each analyte to sensory descriptors, which include odor characteristics such as "Green," "Pungent," and "Malty." The dataset could be useful in food science research.

- ingredient (Categorical): The ingredient from which the analytes were measured (e.g., "barley," "corn"). This refers to the ingredient from which the sample was collected.
- replicate (Categorical): The replicate number of the sample for the experiment, indicating the repetition of the measurement for consistency.
- analyte (Categorical): The specific volatile compound or chemical measured in the ingredient (e.g., "2-Methylpropanal").
- analyte_class (Categorical): The chemical classification of the analyte (e.g., "Aldehyde").
- abundance (Continuous): The concentration of the analyte measured in the sample, likely in a quantitative unit such as mg/L.
- analyte_odor (Categorical): A sensory descriptor for the odors associated with the analyte, listed as a combination of descriptors (e.g., "Green; Pungent; Burnt; Malty; Toasted").

## hawaii_aquifers {-}

This dataset represents water quality measurements from various wells within an aquifer system, collected as part of a study on groundwater composition. Researchers measured the abundance of different dissolved elements and compounds, such as silica (SiO2) and chloride (Cl), from different wells in an aquifer. The dataset could be used to assess the chemical profile of groundwater and monitor any changes in water quality over time. Note the absence of geospatial data (latitude and longitude) for certain samples, and note that some samples come from the same well and aquifer but have different latitude and longitude coordinates.

- aquifer_code (Categorical): The code assigned to identify the aquifer system (e.g., "aquifer_1") that the sample came from.
- well_name (Categorical): The name of the well from which the water sample was collected (e.g., "Alewa_Heights_Spring").
- longitude (Continuous): The longitudinal coordinates of the well location from which the sample was take.
- latitude (Continuous): The latitudinal coordinates of the well location from which the sample was take.
- analyte (Categorical): The specific dissolved compound or element measured in the water sample (e.g., "SiO2," "Cl").
- abundance (Continuous): The concentration of the analyte in the water sample, expressed in a quantitative unit (mg/L).

## hops_components {-}

This dataset contains detailed information on various hop varieties used in brewing, including their country of origin, brewing usage (aroma or bittering), and the chemical composition of their essential oils and acids. The dataset serves as a resource for brewers to select hop varieties based on their aroma profiles and chemical content, such as alpha acids and essential oils like humulene and myrcene, which influence the bitterness, flavor, and aroma of beer. The goal of this dataset is to help optimize the hop selection process in brewing for specific flavor profiles and brewing techniques like dry hopping or bittering.

- hop_variety (Categorical): The specific variety of hop used (e.g., "Cascade," "Chinook").
- hop_origin (Categorical): The country of origin for the hop variety (e.g., "USA," "England").
- hop_brewing_usage (Categorical): The primary use of the hop in brewing, either for aroma or bittering, or techniques like dry hopping.
- hop_aroma (Categorical): The sensory description of the hop's aroma profile, which can include terms like "floral," "citrus," or "spicy."
- total_oil (Continuous): The total essential oil content of the hop, typically measured in milliliters per 100 grams.
- alpha_acids (Continuous): The percentage of alpha acids in the hop, which contribute to the bitterness of the beer.
- beta_acids (Continuous): The percentage of beta acids in the hop, which also contribute to the bitterness but degrade more slowly over time.
- humulene (Continuous): A compound contributing to the hop's woody, earthy aroma.
- myrcene (Continuous): A compound that contributes to the hop's citrus and floral aroma.
- humulone (Continuous): Another compound found in hop oils, related to bitterness and aroma.
- caryophyllene (Continuous): A compound contributing to spicy, peppery aromas.
- farnesene (Continuous): A compound contributing to green, woody, or fruity aromas.

<!-- start R FAQ -->
# r faq {-}

## Updating R and R Packages {-}

Close RStudio, open the plain R GUI, then run the following:

On Mac:


``` r
install.packages('remotes') #assuming it is not remotes installed
remotes::install_github('andreacirilloac/updateR')
updateR::updateR()
```

On PC:


``` r
install.packages("installr")
installr::updateR()
```

## ordering {-}

A list of numeric element has an inherent order to it: -inf -> +inf. A list of character element also has an inherent order to it: A -> Z, or if it's a mixed number and letter list (which is interpreted by R as a character list): 0 -> 9 -> A -> Z.

However, there are cases where we will want a list of character elements to have some order other than A -> Z. In these cases, we want to convert the list of character elements into a list of factor elements. Factors are lists of character elements that have an inherent order that is not A -> Z. For example, in the plot below, the y axis is not, perhaps, in the "correct" order:


``` r
ggplot(periodic_table) +
  geom_point(aes(y = group_number, x = atomic_mass_rounded))
```

<img src="index_files/figure-html/unnamed-chunk-38-1.png" width="100%" style="display: block; margin: auto;" />

How do we fix this? We need to convert the column `group_number` into a list of factors that have the correct order (see below). For this, we will use the command `factor`, which will accept an argument called `levels` in which we can define the order the the characters should be in:


``` r
periodic_table$group_number <- factor(
  periodic_table$group_number,
  levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "lanthanides", "actinides")
)

periodic_table
## # A tibble: 118 × 41
##    atomic_number element_name atomic_symbol group_number
##            <dbl> <chr>        <chr>         <fct>       
##  1             1 hydrogen     H             1           
##  2             2 helium       He            18          
##  3             3 lithium      Li            1           
##  4             4 beryllium    Be            2           
##  5             5 boron        B             13          
##  6             6 carbon       C             14          
##  7             7 nitrogen     N             15          
##  8             8 oxygen       O             16          
##  9             9 fluorine     F             17          
## 10            10 neon         Ne            18          
## # ℹ 108 more rows
## # ℹ 37 more variables: period <dbl>,
## #   atomic_mass_rounded <dbl>, melting_point_C <dbl>,
## #   boiling_point_C <dbl>, state_at_RT <chr>,
## #   density_g_per_mL <dbl>,
## #   electronegativity_pauling <dbl>,
## #   first_ionization_poten_eV <dbl>, …
```

Notice that now when we look at the type of data that is contained in the column `group_number` it says "<fct>". This is great! It means we have converted that column into a list of factors, instead of characters. Now what happens when we make our plot?


``` r
ggplot(periodic_table) +
  geom_point(aes(y = group_number, x = atomic_mass_rounded))
```

<img src="index_files/figure-html/unnamed-chunk-40-1.png" width="100%" style="display: block; margin: auto;" />

VICTORY!

## column manipulation {-}


How to select specific columns:


``` r
alaska_lake_data %>%
  select(water_temp, pH)
## # A tibble: 220 × 2
##    water_temp    pH
##         <dbl> <dbl>
##  1       6.46  7.69
##  2       6.46  7.69
##  3       6.46  7.69
##  4       6.46  7.69
##  5       6.46  7.69
##  6       6.46  7.69
##  7       6.46  7.69
##  8       6.46  7.69
##  9       6.46  7.69
## 10       6.46  7.69
## # ℹ 210 more rows
```

How to remove certain columns:

``` r
alaska_lake_data %>%
  select(!water_temp)
## # A tibble: 220 × 6
##    lake            park     pH element mg_per_L element_type
##    <chr>           <chr> <dbl> <chr>      <dbl> <chr>       
##  1 Devil_Mountain… BELA   7.69 C          3.4   bound       
##  2 Devil_Mountain… BELA   7.69 N          0.028 bound       
##  3 Devil_Mountain… BELA   7.69 P          0     bound       
##  4 Devil_Mountain… BELA   7.69 Cl        10.4   free        
##  5 Devil_Mountain… BELA   7.69 S          0.62  free        
##  6 Devil_Mountain… BELA   7.69 F          0.04  free        
##  7 Devil_Mountain… BELA   7.69 Br         0.02  free        
##  8 Devil_Mountain… BELA   7.69 Na         8.92  free        
##  9 Devil_Mountain… BELA   7.69 K          1.2   free        
## 10 Devil_Mountain… BELA   7.69 Ca         5.73  free        
## # ℹ 210 more rows
```

## user color palettes {-}

Suppose we want to create a specific color palette for each pack in `alaska_lake_data`. There are three unique parks:


``` r
unique(alaska_lake_data$park)
## [1] "BELA" "GAAR" "NOAT"
```

First we define the colors we want:


``` r
custom_colors_for_lakes <- c("#1a9850", "#ffffbf", "#d73027")
custom_colors_for_lakes
## [1] "#1a9850" "#ffffbf" "#d73027"
```

Then we name that vector according to which park we want to be which color:


``` r
names(custom_colors_for_lakes) <- c("GAAR", "NOAT", "BELA")
custom_colors_for_lakes
##      GAAR      NOAT      BELA 
## "#1a9850" "#ffffbf" "#d73027"
```

Now we feed that object to the `values` argument of scale_color_manual (or scale_fill_manual, if you want fill):


``` r
ggplot(alaska_lake_data) + 
  geom_point(aes(x = pH, y = water_temp, fill = park), size = 5, shape = 21, color = "black") +
  scale_fill_manual(values = custom_colors_for_lakes) +
  theme_classic()
```

<img src="index_files/figure-html/unnamed-chunk-46-1.png" width="100%" style="display: block; margin: auto;" />

# templates {-}

## matrix analyses

### basic `runMatrixAnalyses()` template


``` r

runMatrixAnalyses(   
  data = NULL,
  analysis = c("hclust", "pca", "pca_ord", "pca_dim"),
  columns_w_values_for_single_analyte = NULL,
  columns_w_sample_ID_info = NULL
)
```

### advanced `runMatrixAnalysis()` template


``` r

runMatrixAnalysis(
  data, # data to use for analysis
  analysis = c(
      "pca", "pca_ord", "pca_dim", # PCA
      "mca", "mca_ord", "mca_dim", # MCA (PCA on categorical data)
      "mds", "mds_ord", "mds_dim", # MDS
      "tsne", "dbscan", "kmeans", # Clustering
      "hclust", "hclust_phylo" # Hierarchical clustering
  ),
  parameters = NULL,
  column_w_names_of_multiple_analytes = NULL,
  column_w_values_for_multiple_analytes = NULL,
  columns_w_values_for_single_analyte = NULL,
  columns_w_additional_analyte_info = NULL,
  columns_w_sample_ID_info = NULL,
  transpose = FALSE, # default = FALSE, this chooses whether to transpose the data
  distance_method = c( # the distance metric to use in computing a distance matrix
    "euclidean", "maximum",
    "manhattan", "canberra",
    "binary", "minkowski",
    "coeff_unlike"
  ),
  agglomeration_method = c( # the clustering method to use in heirarchical clustering
      "ward.D2", "ward.D", "single", "complete",
      "average", # (= UPGMA)
      "mcquitty", # (= WPGMA)
      "median", # (= WPGMC)
      "centroid" # (= UPGMC)
  ),
  tree_method = c("nj"),
  unknown_sample_ID_info = NULL,
  components_to_return = 2, # how many principal components to return
  scale_variance = NULL, ## default = TRUE, except for hclust, then default = FALSE
  na_replacement = c("mean", "none", "zero", "drop"), # default = "mean", this chooses what to do with missing values
  output_format = c("wide", "long"), # default = "wide", this chooses whether to output a wide or long format
)
```

# loading analyzeGCMSdata {-}

This page explains how to load a simple application for integrating and analyzing GC-MS data. With the app, you can analyze .CDF.csv files. CDF.csv files contain essentially all the data from a GC-MS run, and can be exported from most GC-MS systems using the software provided by the manufacturer. To run the application, use the following guidelines:

1. Create a new folder on your hard drive and place your CDF.csv file into that folder. It doesn't matter what the name of that folder is, but it must not contain special characters (including a space ` ` in the name). For example, if my CDF.csv file is called "sorghum_bicolor.CDF.csv", then I might create a folder called `gc_data` on my hard drive, and place the "sorghum_bicolor.CDF.csv" file in that folder.

2. In R or RStudio, run the source command shown below. You will need to do this every time you re-open R or RStudio. This command will load the `analyzeGCMSdata` function into your R or RStudio environment. Note that you will need to be connected to the internet for this to work. The first time you run this command, a bunch of tools required for the analysis will be installed. It might take a while! Later runs of the command will not have this install portion - the app will just open straight away.


``` r
source("https://thebustalab.github.io/phylochemistry/gcms.R")
```

If the command was run successfully you should see something like:
<img src="https://thebustalab.github.io/integrated_bioanalytics/images/success.png" width="100%" style="display: block; margin: auto;" />
 \

3. In R or RStudio, run the `analyzeGCMSdata` command on the *folder* that contains your CDF.csv file. Do not run the command on the CDF.csv file itself, that will not work. For example, if my CDF.csv file is called "sorghum_bicolor.CDF.csv", and is inside the folder called `gc_data`, then I would run the following:
 \

If you are on a Mac, *use single forward slashes*. For example:

``` r
analyzeGCMSdata("/Volumes/My_Drive/gc_data")
```
 \

If you are on a PC, you may need to use double back slashes. For example:

``` r
analyzeGCMSdata("C:\\Users\\My_Profile\\gc_data")
```
 \

The first time you open your CDF.csv datafile, it may take a while to load. Once the new RShiny window opens, press shift+q to load the chromatogram(s).

# using analyzeGCMSdata {-}

## basic usage of analyteGCMSdata {-}

As a reference, below are the key commands used to operate the integration app. This is the information that is covered in the overview video.

To control the chromatogram window:

* shift + q = update
* shift + a = add selected peak
* shift + r = remove selected peak
* shift + z = save table

To control the mass spectrum window:

* shift+1 = extract mass spectra from highlighted chromatogram region, plot average mass spectrum in panel 1.
* shift+2 = refresh mass spectrum in panel 1. This is used for zooming in on a region of the mass spectrum that you have highlighted. A spectrum needs to first be extracted for this to be possible.
* shift+3 = extract mass spectra from highlighted chromatogram region, subtract their average from the mass spectrum in panel 1.
