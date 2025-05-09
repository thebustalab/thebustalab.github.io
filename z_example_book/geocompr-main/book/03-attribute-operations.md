# Attribute data operations {#attr}

## Prerequisites {-}

- This chapter requires the following packages to be installed and attached: 


```r
library(sf)      # vector data package introduced in Chapter 2
library(terra)   # raster data package introduced in Chapter 2
library(dplyr)   # tidyverse package for data frame manipulation
```

- It also relies on **spData**, which loads datasets used in the code examples of this chapter:


```r
library(spData)  # spatial data package introduced in Chapter 2
#> To access larger datasets in this package, install the spDataLarge
#> package with: `install.packages('spDataLarge',
#> repos='https://nowosad.github.io/drat/', type='source')`
```

## Introduction

Attribute data is non-spatial information associated with geographic (geometry) data.
A bus stop provides a simple example: its position would typically be represented by latitude and longitude coordinates (geometry data), in addition to its name.
The [Elephant & Castle / New Kent Road](https://www.openstreetmap.org/relation/6610626) stop in London, for example has coordinates of -0.098 degrees longitude and 51.495 degrees latitude which can be represented as `POINT (-0.098 51.495)` in the `sfc` representation described in Chapter \@ref(spatial-class).
Attributes such as the name *attribute*\index{attribute} of the POINT feature (to use Simple Features terminology) are the topic of this chapter.



Another example is the elevation value (attribute) for a specific grid cell in raster data.
Unlike the vector data model, the raster data model stores the coordinate of the grid cell indirectly, meaning the distinction between attribute and spatial information is less clear.
To illustrate the point, think of a pixel in the 3^rd^ row and the 4^th^ column of a raster matrix.
Its spatial location is defined by its index in the matrix: move from the origin four cells in the x direction (typically east and right on maps) and three cells in the y direction (typically south and down).
The raster's *resolution* defines the distance for each x- and y-step which is specified in a *header*.
The header is a vital component of raster datasets which specifies how pixels relate to geographic coordinates (see also Chapter \@ref(spatial-operations)).

This teaches how to manipulate geographic objects based on attributes such as the names of bus stops in a vector dataset and elevations of pixels in a raster dataset.
For vector data, this means techniques such as subsetting and aggregation (see Sections \@ref(vector-attribute-subsetting) and \@ref(vector-attribute-aggregation)).
Sections \@ref(vector-attribute-joining) and \@ref(vec-attr-creation) demonstrate how to join data onto simple feature objects using a shared ID and how to create new variables, respectively.
Each of these operations has a spatial equivalent:
the `[` operator in base R, for example, works equally for subsetting objects based on their attribute and spatial objects; you can also join attributes in two geographic datasets using spatial joins.
This is good news: skills developed in this chapter are cross-transferable.
Chapter \@ref(spatial-operations) extends the methods presented here to the spatial world.

After a deep dive into various types of *vector* attribute operations in the next section, *raster* attribute data operations are covered in Section \@ref(manipulating-raster-objects), which demonstrates how to create raster layers containing continuous and categorical attributes and extracting cell values from one or more layer (raster subsetting). 
Section \@ref(summarizing-raster-objects) provides an overview of 'global' raster operations which can be used to summarize entire raster datasets.

## Vector attribute manipulation

Geographic vector datasets are well supported in R thanks to the `sf` class, which extends base R's `data.frame`.
Like data frames, `sf` objects have one column per attribute variable (such as 'name') and one row per observation or *feature* (e.g., per bus station).
`sf` objects differ from basic data frames because they have a `geometry` column of class `sfc` which can contain a range of geographic entities (single and 'multi' point, line, and polygon features) per row.
This was described in Chapter \@ref(spatial-class), which demonstrated how *generic methods* such as `plot()` and `summary()` work with `sf` objects.
**sf** also provides generics that allow `sf` objects to behave like regular data frames, as shown by printing the class's methods:


```r
methods(class = "sf") # methods for sf objects, first 12 shown
```


```r
#>  [1] aggregate             cbind                 coerce               
#>  [4] initialize            merge                 plot                 
#>  [7] print                 rbind                 [                    
#> [10] [[<-                  $<-                   show                 
```



Many of these (`aggregate()`, `cbind()`, `merge()`, `rbind()` and `[`) are for manipulating data frames.
`rbind()`, for example, is binds rows two data frames together, one 'on top' of the other.
`$<-` creates new columns. 
A key feature of `sf` objects is that they store spatial and non-spatial data in the same way, as columns in a `data.frame`.

\BeginKnitrBlock{rmdnote}<div class="rmdnote">The geometry column of `sf` objects is typically called `geometry` or `geom` but any name can be used.
The following command, for example, creates a geometry column named g:
  
`st_sf(data.frame(n = world$name_long), g = world$geom)`

This enables geometries imported from spatial databases to have a variety of names such as `wkb_geometry` and `the_geom`.</div>\EndKnitrBlock{rmdnote}

`sf` objects can also extend the `tidyverse` classes for data frames, `tibble` and `tbl`.
\index{tidyverse (package)}.
Thus **sf** enables the full power of R's data analysis capabilities to be unleashed on geographic data, whether you use base R or tidyverse functions for data analysis.
\index{tibble}
(See [`Rdatatable/data.table#2273`](https://github.com/Rdatatable/data.table/issues/2273) for discussion of compatibility between `sf` objects and the fast `data.table` package.)
Before using these capabilities it is worth re-capping how to discover the basic properties of vector data objects.
Let's start by using base R functions to learn about the `world` dataset from the **spData** package:


```r
class(world) # it's an sf object and a (tidy) data frame
#> [1] "sf"         "tbl_df"     "tbl"        "data.frame"
dim(world)   # it is a 2 dimensional object, with 177 rows and 11 columns
#> [1] 177  11
```

`world` contains ten non-geographic columns (and one geometry list column) with almost 200 rows representing the world's countries.
The function `st_drop_geometry()` keeps only the attributes data of an `sf` object, in other words removing its geometry:


```r
world_df = st_drop_geometry(world)
class(world_df)
#> [1] "tbl_df"     "tbl"        "data.frame"
ncol(world_df)
#> [1] 10
```

Dropping the geometry column before working with attribute data can be useful; data manipulation processes can run faster when they work only on the attribute data and geometry columns are not always needed.
For most cases, however, it makes sense to keep the geometry column, explaining why the column is 'sticky' (it remains after most attribute operations unless specifically dropped).
Non-spatial data operations on `sf` objects only change an object's geometry when appropriate (e.g., by dissolving borders between adjacent polygons following aggregation).
Becoming skilled at geographic attribute data manipulation means becoming skilled at manipulating data frames.

For many applications, the tidyverse\index{tidyverse (package)} package **dplyr** offers an effective approach for working with data frames.
Tidyverse compatibility is an advantage of **sf** over its predecessor **sp**, but there are some pitfalls to avoid (see the supplementary `tidyverse-pitfalls` vignette at [geocompr.github.io](https://geocompr.github.io/geocompkg/articles/tidyverse-pitfalls.html) for details).

### Vector attribute subsetting

Base R subsetting methods include the operator `[` and the function `subset()`.
The key **dplyr** subsetting functions are  `filter()` and `slice()` for subsetting rows, and `select()` for subsetting columns.
Both approaches preserve the spatial components of attribute data in `sf` objects, while using the operator `$` or the **dplyr** function `pull()` to return a single attribute column as a vector will lose the attribute data, as we will see.
\index{attribute!subsetting}
This section focuses on subsetting `sf` data frames; for further details on subsetting vectors and non-geographic data frames we recommend reading section section [2.7](https://cran.r-project.org/doc/manuals/r-release/R-intro.html#Index-vectors) of An Introduction to R [@rcoreteam_introduction_2021] and Chapter [4](https://adv-r.hadley.nz/subsetting.html) of Advanced R Programming [@wickham_advanced_2019], respectively.

The `[` operator can subset both rows and columns. 
Indices placed inside square brackets placed directly after a data frame object name specify the elements to keep.
The command `object[i, j]` means 'return the rows represented by `i` and the columns represented by `j`, where `i` and `j` typically contain integers or `TRUE`s and `FALSE`s (indices can also be character strings, indicating row or column names).
`object[5, 1:3]`, for example, means 'return data containing the 5th row and columns 1 to 3: the result should be a data frame with only 1 row and 3 columns, and a fourth geometry column if it's an `sf` object.
Leaving `i` or `j` empty returns all rows or columns, so `world[1:5, ]` returns the first five rows and all 11 columns.
The examples below demonstrate subsetting with base R.
Guess the number of rows and columns in the `sf` data frames returned by each command and check the results on your own computer (see the end of the chapter for more exercises):


```r
world[1:6, ]    # subset rows by position
world[, 1:3]    # subset columns by position
world[1:6, 1:3] # subset rows and columns by position
world[, c("name_long", "pop")] # columns by name
world[, c(T, T, F, F, F, F, F, T, T, F, F)] # by logical indices
world[, 888] # an index representing a non-existent column
```



A demonstration of the utility of using `logical` vectors for subsetting is shown in the code chunk below.
This creates a new object, `small_countries`, containing nations whose surface area is smaller than 10,000 km^2^:


```r
i_small = world$area_km2 < 10000
summary(i_small) # a logical vector
#>    Mode   FALSE    TRUE 
#> logical     170       7
small_countries = world[i_small, ]
```

The intermediary `i_small` (short for index representing small countries) is a logical vector that can be used to subset the seven smallest countries in the `world` by surface area.
A more concise command, which omits the intermediary object, generates the same result:


```r
small_countries = world[world$area_km2 < 10000, ]
```

The base R function `subset()` provides another way to achieve the same result:


```r
small_countries = subset(world, area_km2 < 10000)
```

Base R functions are mature, stable and widely used, making them a rock solid choice, especially in contexts where reproducibility and reliability are key.
**dplyr** functions enable 'tidy' workflows which some people (the authors of this book included) find intuitive and productive for interactive data analysis, especially when combined with code editors such as RStudio that enable [auto-completion](https://support.rstudio.com/hc/en-us/articles/205273297-Code-Completion-in-the-RStudio-IDE) of column names.
Key functions for subsetting data frames (including `sf` data frames) with **dplyr** functions are demonstrated below.
<!-- The sentence below seems to be untrue based on the benchmark below. -->
<!-- `dplyr` is also faster than base R for some operations, due to its C++\index{C++} backend. -->
<!-- Something on dbplyr? I've never seen anyone use it regularly for spatial data 'in the wild' so leaving out the bit on integration with dbs for now (RL 2021-10) -->
<!-- The main **dplyr** subsetting functions are `select()`, `slice()`, `filter()` and `pull()`. -->



`select()` selects columns by name or position.
For example, you could select only two columns, `name_long` and `pop`, with the following command:


```r
world1 = dplyr::select(world, name_long, pop)
names(world1)
#> [1] "name_long" "pop"       "geom"
```

Note: as with the equivalent command in base R (`world[, c("name_long", "pop")]`), the sticky `geom` column remains.
`select()` also allows selecting a range of columns with the help of the `:` operator: 


```r
# all columns between name_long and pop (inclusive)
world2 = dplyr::select(world, name_long:pop)
```

You can remove specific columns with the `-` operator:


```r
# all columns except subregion and area_km2 (inclusive)
world3 = dplyr::select(world, -subregion, -area_km2)
```

Subset and rename columns at the same time with the `new_name = old_name` syntax:


```r
world4 = dplyr::select(world, name_long, population = pop)
```

It is worth noting that the command above is more concise than base R equivalent, which requires two lines of code:


```r
world5 = world[, c("name_long", "pop")] # subset columns by name
names(world5)[names(world5) == "pop"] = "population" # rename column manually
```

`select()` also works with 'helper functions' for more advanced subsetting operations, including `contains()`, `starts_with()` and `num_range()` (see the help page with `?select` for details).

Most **dplyr** verbs return a data frame, but you can extract a single column as a vector with `pull()`.
<!-- Note: I have commented out the statement below because it is not true for `sf` objects, it's a bit confusing that the behaviour differs between data frames and `sf` objects. -->
<!-- The subsetting operator in base R (see `?[`), by contrast, tries to return objects in the lowest possible dimension. -->
<!-- This means selecting a single column returns a vector in base R as demonstrated in code chunk below which returns a numeric vector representing the population of countries in the `world`: -->
You can get the same result in base R with the list subsetting operators `$` and `[[`, the three following commands return the same numeric vector:


```r
pull(world, pop)
world$pop
world[["pop"]]
```

<!-- Commenting out the following because it's confusing and covered better in other places (RL, 2021-10) -->
<!-- To turn off this behavior, set the `drop` argument to `FALSE`,  -->





`slice()` is the row-equivalent of `select()`.
The following code chunk, for example, selects rows 1 to 6:


```r
slice(world, 1:6)
```

`filter()` is **dplyr**'s equivalent of base R's `subset()` function.
It keeps only rows matching given criteria, e.g., only countries with and area below a certain threshold, or with a high average of life expectancy, as shown in the following examples:


```r
world7 = filter(world ,area_km2 < 10000) # countries with a small area
world7 = filter(world, lifeExp > 82)      # with high life expectancy
```

The standard set of comparison operators can be used in the `filter()` function, as illustrated in Table \@ref(tab:operators): 




Table: (\#tab:operators)Comparison operators that return Booleans (TRUE/FALSE).

|Symbol                        |Name                            |
|:-----------------------------|:-------------------------------|
|`==`                          |Equal to                        |
|`!=`                          |Not equal to                    |
|`>`, `<`                      |Greater/Less than               |
|`>=`, `<=`                    |Greater/Less than or equal      |
|`&`, <code>&#124;</code>, `!` |Logical operators: And, Or, Not |

### Chaining commands with pipes

Key to workflows using **dplyr** functions is the ['pipe'](http://r4ds.had.co.nz/pipes.html) operator `%>%` (or since R `4.1.0` the native pipe `|>`), which takes its name from the Unix pipe `|` [@grolemund_r_2016].
Pipes enable expressive code: the output of a previous function becomes the first argument of the next function, enabling *chaining*.
This is illustrated below, in which only countries from Asia are filtered from the `world` dataset, next the object is subset by columns (`name_long` and `continent`) and the first five rows (result not shown).


```r
world7 = world %>%
  filter(continent == "Asia") %>%
  dplyr::select(name_long, continent) %>%
  slice(1:5)
```

The above chunk shows how the pipe operator allows commands to be written in a clear order:
the above run from top to bottom (line-by-line) and left to right.
The alternative to `%>%` is nested function calls, which is harder to read:


```r
world8 = slice(
  dplyr::select(
    filter(world, continent == "Asia"),
    name_long, continent),
  1:5)
```

### Vector attribute aggregation

\index{attribute!aggregation}
\index{aggregation}
Aggregation involves summarizing data with one or more 'grouping variables', typically from columns in the data frame to be aggregated (geographic aggregation is covered in the next chapter).
An example of attribute aggregation is calculating the number of people per continent based on country-level data (one row per country).
The `world` dataset contains the necessary ingredients: the columns `pop` and `continent`, the population and the grouping variable, respectively.
The aim is to find the `sum()` of country populations for each continent, resulting in a smaller data frame (aggregation is a form of data reduction and can be a useful early step when working with large datasets).
This can be done with the base R function `aggregate()` as follows:


```r
world_agg1 = aggregate(pop ~ continent, FUN = sum, data = world, na.rm = TRUE)
class(world_agg1)
#> [1] "data.frame"
```

The result is a non-spatial data frame with six rows, one per continent, and two columns reporting the name and population of each continent (see Table \@ref(tab:continents) with results for the top 3 most populous continents).

`aggregate()` is a [generic function](https://adv-r.hadley.nz/s3.html#s3-methods) which means that it behaves differently depending on its inputs. 
**sf** provides the method `aggregate.sf()` which is activated automatically when `x` is an `sf` object and a `by` argument is provided:


```r
world_agg2 = aggregate(world["pop"], list(world$continent), FUN = sum, na.rm = TRUE)
class(world_agg2)
#> [1] "sf"         "data.frame"
nrow(world_agg2)
#> [1] 8
```

The resulting `world_agg2` object is a spatial object containing 8 features representing the continents of the world (and the open ocean).
`group_by() %>% summarize()` is the **dplyr** equivalent of `aggregate()`, with the variable name provided in the `group_by()` function specifying the grouping variable and information on what is to be summarized passed to the `summarize()` function, as shown below:


```r
world_agg3 = world %>%
  group_by(continent) %>%
  summarize(pop = sum(pop, na.rm = TRUE))
```

The approach may seem more complex but it has benefits: flexibility, readability, and control over the new column names.
This flexibility is illustrated in the command below, which calculates not only the population but also the area and number of countries in each continent:


```r
world_agg4  = world %>% 
  group_by(continent) %>%
  summarize(pop = sum(pop, na.rm = TRUE), `area (sqkm)` = sum(area_km2), n = n())
```

In the previous code chunk `pop`, `area (sqkm)` and `n` are column names in the result, and `sum()` and `n()` were the aggregating functions.
These aggregating functions return `sf` objects with rows representing continents and geometries containing the multiple polygons representing each land mass and associated islands (this works thanks to the geometric operation 'union', as explained in Section \@ref(geometry-unions)).

Let's combine what we have learned so far about **dplyr** functions, by chaining multiple commands to summarize attribute data about countries worldwide by continent.
The following command calculates population density (with `mutate()`), arranges continents by the number countries they contain (with `dplyr::arrange()`), and keeps only the 3 most populous continents (with `top_n()`), the result of which is presented in Table \@ref(tab:continents)):


```r
world_agg5 = world %>% 
  st_drop_geometry() %>%                      # drop the geometry for speed
  dplyr::select(pop, continent, area_km2) %>% # subset the columns of interest  
  group_by(continent) %>%                     # group by continent and summarize:
  summarize(Pop = sum(pop, na.rm = TRUE), Area = sum(area_km2), N = n()) %>%
  mutate(Density = round(Pop / Area)) %>%     # calculate population density
  top_n(n = 3, wt = Pop) %>%                  # keep only the top 3
  arrange(desc(N))                            # arrange in order of n. countries
```


Table: (\#tab:continents)The top 3 most populous continents ordered by population density (people per square km).

|continent |        Pop|     Area|  N| Density|
|:---------|----------:|--------:|--:|-------:|
|Africa    | 1154946633| 29946198| 51|      39|
|Asia      | 4311408059| 31252459| 47|     138|
|Europe    |  669036256| 23065219| 39|      29|

\BeginKnitrBlock{rmdnote}<div class="rmdnote">More details are provided in the help pages (which can be accessed via `?summarize` and `vignette(package = "dplyr")` and Chapter 5 of [R for Data Science](http://r4ds.had.co.nz/transform.html#grouped-summaries-with-summarize). </div>\EndKnitrBlock{rmdnote}

###  Vector attribute joining

Combining data from different sources is a common task in data preparation. 
Joins do this by combining tables based on a shared 'key' variable.
**dplyr** has multiple join functions including `left_join()` and `inner_join()` --- see `vignette("two-table")` for a full list.
These function names follow conventions used in the database language [SQL](http://r4ds.had.co.nz/relational-data.html) [@grolemund_r_2016, Chapter 13]; using them to join non-spatial datasets to `sf` objects is the focus of this section.
**dplyr** join functions work the same on data frames and `sf` objects, the only important difference being the `geometry` list column.
The result of data joins can be either an `sf` or `data.frame` object.
The most common type of attribute join on spatial data takes an `sf` object as the first argument and adds columns to it from a `data.frame` specified as the second argument.
\index{join}
\index{attribute!join}

To demonstrate joins, we will combine data on coffee production with the `world` dataset.
The coffee data is in a data frame called `coffee_data` from the **spData** package (see `?coffee_data` for details).
It has 3 columns:
`name_long` names major coffee-producing nations and `coffee_production_2016` and `coffee_production_2017` contain estimated values for coffee production in units of 60-kg bags in each year.
A 'left join', which preserves the first dataset, merges `world` with `coffee_data`:


```r
world_coffee = left_join(world, coffee_data)
#> Joining, by = "name_long"
class(world_coffee)
#> [1] "sf"         "tbl_df"     "tbl"        "data.frame"
```

Because the input datasets share a 'key variable' (`name_long`) the join worked without using the `by` argument (see `?left_join` for details).
The result is an `sf` object identical to the original `world` object but with two new variables (with column indices 11 and 12) on coffee production.
This can be plotted as a map, as illustrated in Figure \@ref(fig:coffeemap), generated with the `plot()` function below:


```r
names(world_coffee)
#>  [1] "iso_a2"                 "name_long"              "continent"             
#>  [4] "region_un"              "subregion"              "type"                  
#>  [7] "area_km2"               "pop"                    "lifeExp"               
#> [10] "gdpPercap"              "geom"                   "coffee_production_2016"
#> [13] "coffee_production_2017"
plot(world_coffee["coffee_production_2017"])
```

<div class="figure" style="text-align: center">
<img src="03-attribute-operations_files/figure-html/coffeemap-1.png" alt="World coffee production (thousand 60-kg bags) by country, 2017. Source: International Coffee Organization." width="100%" />
<p class="caption">(\#fig:coffeemap)World coffee production (thousand 60-kg bags) by country, 2017. Source: International Coffee Organization.</p>
</div>

For joining to work, a 'key variable' must be supplied in both datasets.
By default **dplyr** uses all variables with matching names.
In this case, both `world_coffee` and `world` objects contained a variable called `name_long`, explaining the message `Joining, by = "name_long"`.
In the majority of cases where variable names are not the same, you have two options:

1. Rename the key variable in one of the objects so they match.
2. Use the `by` argument to specify the joining variables.

The latter approach is demonstrated below on a renamed version of `coffee_data`:


```r
coffee_renamed = rename(coffee_data, nm = name_long)
world_coffee2 = left_join(world, coffee_renamed, by = c(name_long = "nm"))
```



Note that the name in the original object is kept, meaning that `world_coffee` and the new object `world_coffee2` are identical.
Another feature of the result is that it has the same number of rows as the original dataset.
Although there are only 47 rows of data in `coffee_data`, all 177 country records are kept intact in `world_coffee` and `world_coffee2`:
rows in the original dataset with no match are assigned `NA` values for the new coffee production variables.
What if we only want to keep countries that have a match in the key variable?
In that case an inner join can be used:


```r
world_coffee_inner = inner_join(world, coffee_data)
#> Joining, by = "name_long"
nrow(world_coffee_inner)
#> [1] 45
```

Note that the result of `inner_join()` has only 45 rows compared with 47 in `coffee_data`.
What happened to the remaining rows?
We can identify the rows that did not match using the `setdiff()` function as follows:


```r
setdiff(coffee_data$name_long, world$name_long)
#> [1] "Congo, Dem. Rep. of" "Others"
```

The result shows that `Others` accounts for one row not present in the `world` dataset and that the name of the `Democratic Republic of the Congo` accounts for the other:
it has been abbreviated, causing the join to miss it.
The following command uses a string matching (regex) function from the **stringr** package to confirm what `Congo, Dem. Rep. of` should be:


```r
(drc = stringr::str_subset(world$name_long, "Dem*.+Congo"))
#> [1] "Democratic Republic of the Congo"
```





To fix this issue, we will create a new version of `coffee_data` and update the name.
`inner_join()`ing the updated data frame returns a result with all 46 coffee-producing nations:


```r
coffee_data$name_long[grepl("Congo,", coffee_data$name_long)] = drc
world_coffee_match = inner_join(world, coffee_data)
#> Joining, by = "name_long"
nrow(world_coffee_match)
#> [1] 46
```

It is also possible to join in the other direction: starting with a non-spatial dataset and adding variables from a simple features object.
This is demonstrated below, which starts with the `coffee_data` object and adds variables from the original `world` dataset.
In contrast with the previous joins, the result is *not* another simple feature object, but a data frame in the form of a **tidyverse** tibble:
the output of a join tends to match its first argument:


```r
coffee_world = left_join(coffee_data, world)
#> Joining, by = "name_long"
class(coffee_world)
#> [1] "tbl_df"     "tbl"        "data.frame"
```

\BeginKnitrBlock{rmdnote}<div class="rmdnote">In most cases, the geometry column is only useful in an `sf` object.
The geometry column can only be used for creating maps and spatial operations if R 'knows' it is a spatial object, defined by a spatial package such as **sf**.
Fortunately, non-spatial data frames with a geometry list column (like `coffee_world`) can be coerced into an `sf` object as follows: `st_as_sf(coffee_world)`. </div>\EndKnitrBlock{rmdnote}

This section covers the majority of joining use cases.
For more information, we recommend @grolemund_r_2016, the [join vignette](https://geocompr.github.io/geocompkg/articles/join.html) in the **geocompkg** package that accompanies this book, and documentation of the **data.table** package.^[
**data.table** is a high-performance data processing package.
Its application to geographic data is covered in a blog post hosted at r-spatial.org/r/2017/11/13/perp-performance.html.
]
Another type of join is a spatial join, covered in the next chapter (Section \@ref(spatial-joining)).

### Creating attributes and removing spatial information {#vec-attr-creation}

Often, we would like to create a new column based on already existing columns.
For example, we want to calculate population density for each country.
For this we need to divide a population column, here `pop`, by an area column, here `area_km2` with unit area in square kilometers.
Using base R, we can type:


```r
world_new = world # do not overwrite our original data
world_new$pop_dens = world_new$pop / world_new$area_km2
```

Alternatively, we can use one of **dplyr** functions - `mutate()` or `transmute()`.
`mutate()` adds new columns at the penultimate position in the `sf` object (the last one is reserved for the geometry):


```r
world %>% 
  mutate(pop_dens = pop / area_km2)
```

The difference between `mutate()` and `transmute()` is that the latter drops all other existing columns (except for the sticky geometry column):


```r
world %>% 
  transmute(pop_dens = pop / area_km2)
```

`unite()` from the **tidyr** package (which provides many useful functions for reshaping datasets, including `pivot_longer()`) pastes together existing columns.
For example, we want to combine the `continent` and `region_un` columns into a new column named `con_reg`.
Additionally, we can define a separator (here: a colon `:`) which defines how the values of the input columns should be joined, and if the original columns should be removed (here: `TRUE`):


```r
world_unite = world %>%
  unite("con_reg", continent:region_un, sep = ":", remove = TRUE)
```

The `separate()` function does the opposite of `unite()`: it splits one column into multiple columns using either a regular expression or character positions.
This function also comes from the **tidyr** package.


```r
world_separate = world_unite %>% 
  separate(con_reg, c("continent", "region_un"), sep = ":")
```



The **dplyr** function `rename()` and the base R function `setNames()` are useful for renaming columns.
The first replaces an old name with a new one.
The following command, for example, renames the lengthy `name_long` column to simply `name`:


```r
world %>% 
  rename(name = name_long)
```

`setNames()` changes all column names at once, and requires a character vector with a name matching each column.
This is illustrated below, which outputs the same `world` object, but with very short names: 




```r
new_names = c("i", "n", "c", "r", "s", "t", "a", "p", "l", "gP", "geom")
world %>% 
  setNames(new_names)
```

It is important to note that attribute data operations preserve the geometry of the simple features.
As mentioned at the outset of the chapter, it can be useful to remove the geometry.
To do this, you have to explicitly remove it.
Hence, an approach such as `select(world, -geom)` will be unsuccessful and you should instead use `st_drop_geometry()`.^[
`st_geometry(world_st) = NULL` also works to remove the geometry from `world`, but overwrites the original object.
]


```r
world_data = world %>% st_drop_geometry()
class(world_data)
#> [1] "tbl_df"     "tbl"        "data.frame"
```

## Manipulating raster objects
<!--jn-->

In contrast to the vector data model underlying simple features (which represents points, lines and polygons as discrete entities in space), raster data represent continuous surfaces.
This section shows how raster objects work by creating them *from scratch*, building on Section \@ref(an-introduction-to-terra).
Because of their unique structure, subsetting and other operations on raster datasets work in a different way, as demonstrated in Section \@ref(raster-subsetting).
\index{raster!manipulation}

The following code recreates the raster dataset used in Section \@ref(raster-classes), the result of which is illustrated in Figure \@ref(fig:cont-raster).
This demonstrates how the `rast()` function works to create an example raster named `elev` (representing elevations).


```r
elev = rast(nrows = 6, ncols = 6, resolution = 0.5, 
            xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5,
            vals = 1:36)
```

The result is a raster object with 6 rows and 6 columns (specified by the `nrow` and `ncol` arguments), and a minimum and maximum spatial extent in x and y direction (`xmin`, `xmax`, `ymin`, `ymax`).
The `vals` argument sets the values that each cell contains: numeric data ranging from 1 to 36 in this case.
Raster objects can also contain categorical values of class `logical` or `factor` variables in R.
The following code creates a raster representing grain sizes (Figure \@ref(fig:cont-raster)):


```r
grain_order = c("clay", "silt", "sand")
grain_char = sample(grain_order, 36, replace = TRUE)
grain_fact = factor(grain_char, levels = grain_order)
grain = rast(nrows = 6, ncols = 6, resolution = 0.5, 
             xmin = -1.5, xmax = 1.5, ymin = -1.5, ymax = 1.5,
             vals = grain_fact)
```



The raster object stores the corresponding look-up table or "Raster Attribute Table" (RAT) as a list of data frames, which can be viewed with `cats(grain)` (see `?cats()` for more information).
Each element of this list is a layer of the raster.
It is also possible to use the function `levels()` for retrieving and adding new or replacing existing factor levels:


```r
levels(grain)[[1]] = c(levels(grain)[[1]], wetness = c("wet", "moist", "dry"))
levels(grain)
#> [[1]]
#> [1] "clay"  "silt"  "sand"  "wet"   "moist" "dry"
```

<div class="figure" style="text-align: center">
<img src="03-attribute-operations_files/figure-html/cont-raster-1.png" alt="Raster datasets with numeric (left) and categorical values (right)." width="100%" />
<p class="caption">(\#fig:cont-raster)Raster datasets with numeric (left) and categorical values (right).</p>
</div>

\BeginKnitrBlock{rmdnote}<div class="rmdnote">Categorical raster objects can also store information about the colors associated with each value using a color table.
The color table is a data frame with three (red, green, blue) or four (alpha) columns, where each row relates to one value.
Color tables in **terra** can be viewed or set with the `coltab()` function (see `?coltab`).
Importantly, saving a raster object with a color table to a file (e.g., GeoTIFF) will also save the color information.</div>\EndKnitrBlock{rmdnote}

### Raster subsetting

Raster subsetting is done with the base R operator `[`, which accepts a variety of inputs:
\index{raster!subsetting}

- Row-column indexing
- Cell IDs
- Coordinates (see Section \@ref(spatial-raster-subsetting))
- Another spatial object (see Section \@ref(spatial-raster-subsetting))

Here, we only show the first two options since these can be considered non-spatial operations.
If we need a spatial object to subset another or the output is a spatial object, we refer to this as spatial subsetting.
Therefore, the latter two options will be shown in the next chapter (see Section \@ref(spatial-raster-subsetting)).

The first two subsetting options are demonstrated in the commands below ---
both return the value of the top left pixel in the raster object `elev` (results not shown):


```r
# row 1, column 1
elev[1, 1]
# cell ID 1
elev[1]
```

Subsetting of multi-layered raster objects will return the cell value(s) for each layer.
For example, `c(elev, grain)[1]` returns a data frame with one row and two columns --- one for each layer.
To extract all values or complete rows, you can also use `values()`.

Cell values can be modified by overwriting existing values in conjunction with a subsetting operation.
The following code chunk, for example, sets the upper left cell of `elev` to 0 (results not shown):


```r
elev[1, 1] = 0
elev[]
```

Leaving the square brackets empty is a shortcut version of `values()` for retrieving all values of a raster.
Multiple cells can also be modified in this way:


```r
elev[1, c(1, 2)] = 0
```

Replacing values of multilayered rasters can be done with a matrix with as many columns as layers and rows as replaceable cells (results not shown):


```r
two_layers = c(grain, elev) 
two_layers[1] = cbind(c(0), c(4))
two_layers[]
```

### Summarizing raster objects

**terra** contains functions for extracting descriptive statistics\index{statistics} for entire rasters.
Printing a raster object to the console by typing its name returns minimum and maximum values of a raster.
`summary()` provides common descriptive statistics\index{statistics} -- minimum, maximum, quartiles and number of `NA`s for continuous rasters and a number of cells of each class for categorical rasters.
Further summary operations such as the standard deviation (see below) or custom summary statistics can be calculated with `global()`. 
\index{raster!summarizing}


```r
global(elev, sd)
```

\BeginKnitrBlock{rmdnote}<div class="rmdnote">If you provide the `summary()` and `global()` functions with a multi-layered raster object, they will summarize each layer separately, as can be illustrated by running: `summary(c(elev, grain))`.</div>\EndKnitrBlock{rmdnote}

Additionally, the `freq()` function allows to get the frequency table of categorical values.

Raster value statistics can be visualized in a variety of ways.
Specific functions such as `boxplot()`, `density()`, `hist()` and `pairs()` work also with raster objects, as demonstrated in the histogram created with the command below (not shown):


```r
hist(elev)
```

In case the desired visualization function does not work with raster objects, one can extract the raster data to be plotted with the help of `values()` (Section \@ref(raster-subsetting)).
\index{raster!values}

Descriptive raster statistics belong to the so-called global raster operations.
These and other typical raster processing operations are part of the map algebra scheme, which are covered in the next chapter (Section \@ref(map-algebra)).

<div class="rmdnote">
<p>Some function names clash between packages (e.g., a function with the name <code>extract()</code> exist in both <strong>terra</strong> and <strong>tidyr</strong> packages). In addition to not loading packages by referring to functions verbosely (e.g., <code>tidyr::extract()</code>), another way to prevent function names clashes is by unloading the offending package with <code>detach()</code>. The following command, for example, unloads the <strong>terra</strong> package (this can also be done in the <em>package</em> tab which resides by default in the right-bottom pane in RStudio): <code>detach("package:terra", unload = TRUE, force = TRUE)</code>. The <code>force</code> argument makes sure that the package will be detached even if other packages depend on it. This, however, may lead to a restricted usability of packages depending on the detached package, and is therefore not recommended.</p>
</div>

## Exercises


For these exercises we will use the `us_states` and `us_states_df` datasets from the **spData** package.
You must have attached the package, and other packages used in the attribute operations chapter (**sf**, **dplyr**, **terra**) with commands such as `library(spData)` before attempting these exercises:

```r
library(sf)
library(dplyr)
library(terra)
library(spData)
data(us_states)
data(us_states_df)
```

`us_states` is a spatial object (of class `sf`), containing geometry and a few attributes (including name, region, area, and population) of states within the contiguous United States.
`us_states_df` is a data frame (of class `data.frame`) containing the name and additional variables (including median income and poverty level, for the years 2010 and 2015) of US states, including Alaska, Hawaii and Puerto Rico.
The data comes from the United States Census Bureau, and is documented in `?us_states` and `?us_states_df`.

E1. Create a new object called `us_states_name` that contains only the `NAME` column from the `us_states` object using either base R (`[`) or tidyverse (`select()`) syntax.
What is the class of the new object and what makes it geographic?



E2. Select columns from the `us_states` object which contain population data.
Obtain the same result using a different command (bonus: try to find three ways of obtaining the same result).
Hint: try to use helper functions, such as `contains` or `starts_with` from **dplyr** (see `?contains`).

E3. Find all states with the following characteristics (bonus find *and* plot them):

- Belong to the Midwest region.
- Belong to the West region, have an area below 250,000 km^2^ *and* in 2015 a population greater than 5,000,000 residents (hint: you may need to use the function `units::set_units()` or `as.numeric()`).
- Belong to the South region, had an area larger than 150,000 km^2^ or a total population in 2015 larger than 7,000,000 residents.

E4. What was the total population in 2015 in the `us_states` dataset?
What was the minimum and maximum total population in 2015?

E5. How many states are there in each region?

E6. What was the minimum and maximum total population in 2015 in each region?
What was the total population in 2015 in each region?

E7. Add variables from `us_states_df` to `us_states`, and create a new object called `us_states_stats`.
What function did you use and why?
Which variable is the key in both datasets?
What is the class of the new object?

E8. `us_states_df` has two more rows than `us_states`.
How can you find them? (hint: try to use the `dplyr::anti_join()` function)

E9. What was the population density in 2015 in each state?
What was the population density in 2010 in each state?

E10. How much has population density changed between 2010 and 2015 in each state?
Calculate the change in percentages and map them.

E11. Change the columns' names in `us_states` to lowercase. (Hint: helper functions - `tolower()` and `colnames()` may help.)

E12. Using `us_states` and `us_states_df` create a new object called `us_states_sel`.
The new object should have only two variables - `median_income_15` and `geometry`.
Change the name of the `median_income_15` column to `Income`.

E13. Calculate the change in the number of residents living below the poverty level between 2010 and 2015 for each state. (Hint: See ?us_states_df for documentation on the poverty level columns.)
Bonus: Calculate the change in the *percentage* of residents living below the poverty level in each state.

E14. What was the minimum, average and maximum state's number of people living below the poverty line in 2015 for each region?
Bonus: What is the region with the largest increase in people living below the poverty line?

E15. Create a raster from scratch with nine rows and columns and a resolution of 0.5 decimal degrees (WGS84).
Fill it with random numbers.
Extract the values of the four corner cells. 



E16. What is the most common class of our example raster `grain` (hint: `modal`)?



E17. Plot the histogram and the boxplot of the `dem.tif` file from the **spDataLarge** package (`system.file("raster/dem.tif", package = "spDataLarge")`). 
