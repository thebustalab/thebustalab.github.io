# Reprojecting geographic data {#reproj-geo-data}

## Prerequisites {-}

- This chapter requires the following packages (**lwgeom** is also used, but does not need to be attached):


```r
library(sf)
library(terra)
library(dplyr)
library(spData)
library(spDataLarge)
#> Warning: package 'spDataLarge' was built under R version 4.1.2
```

## Introduction {#reproj-intro}

Section \@ref(crs-intro) introduced coordinate reference systems (CRSs) and demonstrated their importance.
This chapter goes further, highlighting specific issues that can arise due to ignoring CRSs, and demonstrating how to **set** coordinate systems and *transform* geographic data from one CRS to another.
\index{CRS!geographic} 
\index{CRS!projected} 
As illustrated in Figure \@ref(fig:vectorplots) from that earlier chapter, there are two types of CRSs: *geographic* ('lon/lat', with units in degrees longitude and latitude) and *projected* (typically with units of meters from a datum).
In many projects there is no need to worry about, let alone convert between, different CRSs.
It is important to know if your data is in a projected or geographic coordinate system, and the consequences of this for geometry operations.
However, if you know the CRS of your data and the consequences for geometry operations (covered in the next section), CRS should *just work*: knowledge of CRSs is often most important when things go wrong.
Having a clearly defined project CRS that all project data is in (or is converted into), plus understanding how and why to use different CRSs, can ensure that things don't go wrong. 

This chapter teaches the fundamentals of CRSs, demonstrates the consequences of using different CRSs (including what can go wrong), and how to 'reproject' datasets from one coordinate system to another.
The next section introduces CRSs in R, followed by Section \@ref(crs-in-r) which shows how to get and set CRSs associated with spatial objects. 
Section \@ref(geom-proj) demonstrates the importance of knowing what CRS your data is in with reference to a worked example of creating buffers.
The questions of when to reproject and which CRS to use are covered in Section \@ref(whenproject) and Section \@ref(which-crs), respectively.
Reprojecting vector and raster objects is covered in sections \@ref(reproj-vec-geom) and \@ref(reproj-ras).
Modifying map projections is covered in Section \@ref(mapproj).


## CRSs in R {#crs-in-r}

\index{CRS!EPSG}
\index{CRS!WKT2}
\index{CRS!proj4string}
Spatial R packages support a wide range of CRSs and they use the long-established [PROJ](https://proj.org) library.
Two recommend ways to describe CRSs in R are (a) Spatial Reference System Identifier (SRID) or (b) well-known text (known as WKT2^[
Several WKT dialects were created to describe CRSs, including ESRI WKT, GDAL WKT1, and the current WKT2:2018 [@lott_geographic_2015]]) definitions.
Both of these approaches have advantages and disadvantages. 

A SRID is a unique value used to identify coordinate reference system definitions in a form of *AUTHORITY:CODE*.
The most popular registry of SRIDs is *EPSG*, however, other registries, such as *ESRI* or *OGR*, exist.
For example, *EPSG:4326* represents the latitude/longitude WGS84 CRS, and *ESRI:54030* - Robinson projection.
SRIDs are usually short and therefore easier to remember. 
Each SRID is associated with a well-known text (WKT2) definition of the coordinate reference system. 

A WKT2 describes coordinate reference systems (CRSs) and coordinates operations between them in the form of well-known text strings.
It is exhaustive, detailed, and precise, allowing for unambiguous CRSs storage and transformations.
It consists of all information about any given CRS, including its datum and ellipsoid, prime meridian, projection, units, etc.
This feature also makes the WKT2 approach more complicated and usually too complex to be manually defined.

The `wkt` component stands for '**w**ell-**k**nown **t**ext representation of coordinate reference systems'.
The language was developed as an open standard by the Open Geospatial Commission (OGC) "for the description of coordinate operations" and is related to the WKT representation of geometries, which was also developed by the OGC and is used when printing vector geometries, as outlined in Section \@ref(geometry).
The full the WKT CRS format specification, the latest version of which was published in 2019 as an internationally agreed standard (ISO 19162:2019), is available in a 132 page document published at [docs.opengeospatial.org](http://docs.opengeospatial.org/is/18-010r7/18-010r7.html).

In the past, the `proj4string` definitions, was the standard way to specify coordinate operations and store CRSs.
These string representations, built on a key=value form (e.g, `+proj=longlat +datum=WGS84 +no_defs`), are, however, currently discouraged in most cases.
PROJ version 6 and further still allows to use `proj4string`s to define coordinate operations, but some `proj4string` keys are no longer supported or are not advisable to use (e.g., `+nadgrids`, `+towgs84`, `+k`, `+init=epsg:`) and only three datums (i.e., WGS84, NAD83, and NAD27) can be directly set in `proj4string`.
Importantly, `proj4string`s are not used to store CRSs anymore.
Longer explanations on the recent changes in the PROJ library and why `proj4string` was replaced by `WKT2` can be found in @bivand_progress_2021, Chapter 2 of @pebesma_spatial_2022, and [blog post by Floris Vanderhaeghe](https://inbo.github.io/tutorials/tutorials/spatial_crs_coding/).

## Querying and setting coordinate systems {#crs-setting}

Let's look at how CRSs are stored in R spatial objects and how they can be set.
For this, we need to read-in a vector dataset:


```r
vector_filepath = system.file("shapes/world.gpkg", package = "spData")
new_vector = read_sf(vector_filepath)
```

Our new object, `new_vector`, is a polygon representing a world map data (`?spData::world`).
In **sf** the CRS of an object can be retrieved using `st_crs()`.


```r
st_crs(new_vector) # get CRS
#> Coordinate Reference System:
#>   User input: WGS 84 
#>   wkt:
#>   ...
```

The output is a list containing two main components: 1) `User input` (in this case `EPSG:4326`) and 2) `wkt`.
The `User input` component is simply the text that the user entered to describe the CRS, with `EPSG:` added as a prefix if the CRS was given as an EPSG code.
`crs = 4326` is understood by `sf` as `crs = "EPSG:4326"` and can be used, although we prefer the more explicit character string to prevent ambiguity.

The `input` element is quite flexible, and depending on the input file or user input, can contain SRID representation (e.g., `"EPSG:4326"`), CRS's name (e.g., `"WGS84"`), or even `proj4string` definition.
The `wkt` element stores the WKT2 representation, which is used when saving the object to a file or doing any coordinate operations.
Above, we can see that the `new_vector` object has the WGS84 ellipsoid, uses the Greenwich prime meridian, and the latitude and longitude axis order.
In this case, we also have some additional elements, such as `USAGE` explaining the area suitable for the use of this CRS, and `ID` pointing to the CRS's SRID - `"EPSG:4326"`.

The `st_crs` function also has one helpful feature -- we can retrieve some additional information about the used CRS. 
For example, try to run:

- `st_crs(new_vector)$IsGeographic` to check is the CRS is geographic or not
- `st_crs(new_vector)$units_gdal` to find out the CRS units
- `st_crs(new_vector)$srid` extracts its SRID (when available)
- `st_crs(new_vector)$proj4string` extracts the `proj4string` representation

In cases when a coordinate reference system (CRS) is missing or the wrong CRS is set, the `st_set_crs()` function can be used:


```r
new_vector = st_set_crs(new_vector, "EPSG:4326") # set CRS
```

The second argument in the above function could be either SRID (`"EPSG:4326"` in the example), complete WKT2 representation, `proj4string`, or CRS extracted from the existing object with `st_crs()`.

The `crs()` function can be used to access CRS information from a `SpatRaster` object (note the use of the `cat()` function to print it nicely): 


```r
raster_filepath = system.file("raster/srtm.tif", package = "spDataLarge")
my_rast = rast(raster_filepath)
cat(crs(my_rast)) # get CRS
#> GEOGCRS["WGS 84",
#>     DATUM["World Geodetic System 1984",
#>         ELLIPSOID["WGS 84",6378137,298.257223563,
#>             LENGTHUNIT["metre",1]]],
#>     PRIMEM["Greenwich",0,
#>         ANGLEUNIT["degree",0.0174532925199433]],
#>     CS[ellipsoidal,2],
#>         AXIS["geodetic latitude (Lat)",north,
#>             ORDER[1],
#>             ANGLEUNIT["degree",0.0174532925199433]],
#>         AXIS["geodetic longitude (Lon)",east,
#>             ORDER[2],
#>             ANGLEUNIT["degree",0.0174532925199433]],
#>     ID["EPSG",4326]]
```

The output is the WKT2 representation of CRS. 
The same function, `crs()`, is can be also used to set a CRS for raster objects.


```r
crs(my_rast) = "EPSG:26912" # set CRS
```

Here, we can use either SRID, complete WKT2 representation, `proj4string`, or CRS extracted from other existing object with `crs()`.

Importantly, the `st_crs()` and `crs()` functions do not alter coordinates' values or geometries.
Their role is only to set a metadata information about the object CRS.

In some cases the CRS of a geographic object is unknown, as is the case in the `london` dataset created in the code chunk below, building on the example of London introduced in Section \@ref(vector-data):


```r
london = data.frame(lon = -0.1, lat = 51.5) %>% 
  st_as_sf(coords = c("lon", "lat"))
st_is_longlat(london)
#> [1] NA
```

The output `NA` shows that `sf` does not know what the CRS is and is unwilling to guess (`NA` literally means 'not available').
Unless a CRS is manually specified or is loaded from a source that has CRS metadata, `sf` does not make any explicit assumptions about which coordinate systems, other than to say "I don't know".
This behavior makes sense given the diversity of available CRSs but differs from some approaches, such as the GeoJSON file format specification, which makes the simplifying assumption that all coordinates have a lon/lat CRS: `EPSG:4326`.

A CRS can be added to `sf` objects in three main ways:

- By assigning the CRS to a pre-existing object, e.g. with `st_crs(london) = "EPSG:4326"`.
- By passing a CRS to the `crs` argument in `sf` functions that create geometry objects such as `st_as_sf(... crs = "EPSG:4326")`. The same argument can also be used to set the CRS when creating raster datasets (e.g., `rast(crs = "EPSG:4326")`).
- With the `st_set_crs()`, which returns a version of the data that has a new CRS, an approach that is demonstrated in the following code chunk.


```r
london_geo = st_set_crs(london, "EPSG:4326")
st_is_longlat(london_geo)
#> [1] TRUE
```

<!-- The following example demonstrates how to add CRS metadata to raster datasets. -->
<!-- Todo: add this -->

Datasets without a specified CRS can cause problems: all geographic coordinates have a coordinate system and software can only make good decisions around plotting and and geometry operations if it knows what type of CRS it is working with.

## Geometry operations on projected and unprojected data {#geom-proj}

If no CRS has been set, `sf` uses the GEOS geometry library for many operations.
GEOS is not well suited to lon/lat CRSs, as we will see later in this chapter.
If a CRS has been set, `sf` will use either GEOS or the S2 *spherical geometry engine* depending on the type of CRS.
<!-- Todo: add s2 section -->
<!--jn: s2 section is still missing from the book-->
Since `sf` version 1.0.0, R's ability to work with geographic vector datasets that have lon/lat CRSs has improved substantially, thanks to its integration with S2 introduced in Section \@ref(s2).

To demonstrate the importance of CRSs, we will in this section create a buffer of 100 km around the `london` object created in the previous section.
We will also create a deliberately faulty buffer with a 'distance' of 1 degree, which is roughly equivalent to 100 km (1 degree is about 111 km at the equator).
Before diving into the code, it may be worth skipping briefly ahead to peek at Figure \@ref(fig:crs-buf) to get a visual handle on the outputs that you should be able to reproduce by following the code chunks below.

The first stage is to create three buffers around the `london` and `london_geo` objects created above with boundary distances of 1 degree and 100 km  (or 100,000 m, which can be expressed as `1e5` in scientific notation) from central London:


```r
london_buff_no_crs = st_buffer(london, dist = 1)   # incorrect: no CRS
london_buff_s2 = st_buffer(london_geo, dist = 1e5) # silent use of s2
london_buff_s2_100_cells = st_buffer(london_geo, dist = 1e5, max_cells = 100) 
```

In the first line above, `sf` assumes that the input is projected and generates a result that has a buffer in units of degrees, which is problematic, as we will see.
In the second line, `sf` silently uses the spherical geometry engine S2, introduced in Chapter \@ref(spatial-class), to calculate the extent of the buffer using the default value of `max_cells = 1000` --- set to `100` in line three --- the consequences which will become apparent shortly (see `?s2::s2_buffer_cells` for details).
To highlight the impact of `sf`'s use of the S2 geometry engine for unprojected (geographic) coordinate systems, we will temporarily disable it with the command `sf_use_s2()` (which is on, `TRUE`, by default), in the code chunk below.
Like `london_buff_no_crs`, the new `london_geo` object is a geographic abomination: it has units of degrees, which makes no sense in the vast majority of cases:


```r
sf::sf_use_s2(FALSE)
#> Spherical geometry (s2) switched off
london_buff_lonlat = st_buffer(london_geo, dist = 1) # incorrect result
#> Warning in st_buffer.sfc(st_geometry(x), dist, nQuadSegs, endCapStyle =
#> endCapStyle, : st_buffer does not correctly buffer longitude/latitude data
#> dist is assumed to be in decimal degrees (arc_degrees).
sf::sf_use_s2(TRUE)
#> Spherical geometry (s2) switched on
```

The warning message above hints at issues with performing planar geometry operations on lon/lat data. 
When spherical geometry operations are turned off, with the command `sf::sf_use_s2(FALSE)`, buffers (and other geometric operations) may result in worthless outputs because they use units of latitude and longitude, a poor substitute for proper units of distances such as meters.

\BeginKnitrBlock{rmdnote}<div class="rmdnote">The distance between two lines of longitude, called meridians, is around 111 km at the equator (execute `geosphere::distGeo(c(0, 0), c(1, 0))` to find the precise distance).
This shrinks to zero at the poles.
At the latitude of London, for example, meridians are less than 70 km apart (challenge: execute code that verifies this).
<!-- `geosphere::distGeo(c(0, 51.5), c(1, 51.5))` -->
Lines of latitude, by contrast, are equidistant from each other irrespective of latitude: they are always around 111 km apart, including at the equator and near the poles (see Figures \@ref(fig:crs-buf) to \@ref(fig:wintriproj)).</div>\EndKnitrBlock{rmdnote}

Do not interpret the warning about the geographic (`longitude/latitude`) CRS as "the CRS should not be set": it almost always should be!
It is better understood as a suggestion to *reproject* the data onto a projected CRS.
This suggestion does not always need to be heeded: performing spatial and geometric operations makes little or no difference in some cases (e.g., spatial subsetting).
But for operations involving distances such as buffering, the only way to ensure a good result (without using spherical geometry engines) is to create a projected copy of the data and run the operation on that.
<!--toDo:rl-->
<!-- jn: idea -- maybe it would be add a table somewhere in the book showing which operations are impacted by s2? -->
This is done in the code chunk below:


```r
london_proj = data.frame(x = 530000, y = 180000) %>% 
  st_as_sf(coords = 1:2, crs = "EPSG:27700")
```

The result is a new object that is identical to `london`, but reprojected onto a suitable CRS (the British National Grid, which has an EPSG code of 27700 in this case) that has units of meters.
We can verify that the CRS has changed using `st_crs()` as follows (some of the output has been replaced by `...`):


```r
st_crs(london_proj)
#> Coordinate Reference System:
#>   User input: EPSG:27700 
#>   wkt:
#> PROJCRS["OSGB36 / British National Grid",
#> ...
```

<!--toDo:rl-->
<!-- jn: the next paragraph need to be updated! -->
Notable components of this CRS description include the EPSG code (`EPSG: 27700`), the projection ([transverse Mercator](https://en.wikipedia.org/wiki/Transverse_Mercator_projection), `+proj=tmerc`), the origin (`+lat_0=49 +lon_0=-2`) and units (`+units=m`).^[
For a short description of the most relevant projection parameters and related concepts, see the fourth lecture by Jochen Albrecht hosted at
http://www.geography.hunter.cuny.edu/~jochen/GTECH361/lectures/ and information at https://proj.org/usage/projections.html.
Other great resources on projections are spatialreference.org and progonos.com/furuti/MapProj.
]
The fact that the units of the CRS are meters (rather than degrees) tells us that this is a projected CRS: `st_is_longlat(london_proj)` now returns `FALSE` and geometry operations on `london_proj` will work without a warning, meaning buffers can be produced from it using proper units of distance.
The following line of code creates a buffer around *projected* data of exactly 100 km:


```r
london_buff_projected = st_buffer(london_proj, 1e5)
```

The geometries of the three `london_buff*` objects that *have* a specified CRS created above (`london_buff_s2`, `london_buff_lonlat` and `london_buff_projected`) created in the preceding code chunks are illustrated in Figure \@ref(fig:crs-buf).



<div class="figure" style="text-align: center">
<img src="07-reproj_files/figure-html/crs-buf-1.png" alt="Buffers around London showing results created with the S2 spherical geometry engine on lon/lat data (left), projected data (middle) and lon/lat data without using spherical geometry (right). The left plot illustrates the result of buffering unprojected data with sf, which calls Google's S2 spherical geometry engine by default with `max_cells = 1000` (thin line). The thick 'blocky' line illustrates the result of the same operation with `max_cells = 100`." width="100%" />
<p class="caption">(\#fig:crs-buf)Buffers around London showing results created with the S2 spherical geometry engine on lon/lat data (left), projected data (middle) and lon/lat data without using spherical geometry (right). The left plot illustrates the result of buffering unprojected data with sf, which calls Google's S2 spherical geometry engine by default with `max_cells = 1000` (thin line). The thick 'blocky' line illustrates the result of the same operation with `max_cells = 100`.</p>
</div>

It is clear from Figure \@ref(fig:crs-buf) that buffers based on `s2` and properly projected CRSs are not 'squashed', meaning that every part of the buffer boundary is equidistant to London.
The results that are generated from lon/lat CRSs when `s2` is *not* used, either because the input lacks a CRS or because `sf_use_s2()` is turned off, are heavily distorted, with the result elongated in the north-south axis, highlighting the dangers of using algorithms that assume projected data on lon/lat inputs (as GEOS does).
The results generated using S2 are also distorted, however, although less dramatically.
Both buffer boundaries in Figure \@ref(fig:crs-buf) (left) are jagged, although this may only be apparent or relevant when for the thick boundary representing a buffer created with the `s2` argument `max_cells` set to 100.
<!--toDo:rl-->
<!--jn: maybe it is worth to emphasize that the differences are due to the use of S2 vs GEOS-->
<!--jn: you mention S2 a lot in this section, but not GEOS...-->
The less is that results obtained from lon/lat data via `s2` will be different from results obtained from using projected data, although these differences reduce as the value of `max_cells` increases: the 'right' value for this argument may depend on many factors and the default value 1000 is a reasonable default, balancing speed of computation against resolution of results, in many cases.
In situations where curved boundaries are advantageous, transforming to a projected CRS before buffering (or performing other geometry operations) may be appropriate.

The importance of CRSs (primarily whether they are projected or geographic) and the impacts of `sf`'s default setting to use S2 for buffers on lon/lat data is clear from the example above.
The subsequent sections go into more depth, exploring which CRS to use when projected CRSs *are* needed and the details of reprojecting vector and raster objects.

## When to reproject? {#whenproject}

\index{CRS!reprojection} 
The previous section showed how to set the CRS manually, with `st_set_crs(london, "EPSG:4326")`.
In real world applications, however, CRSs are usually set automatically when data is read-in.
In many projects the main CRS-related task is to *transform* objects, from one CRS into another.
But when should data be transformed? 
And into which CRS?
There are no clear-cut answers to these questions and CRS selection always involves trade-offs [@maling_coordinate_1992].
However, there are some general principles provided in this section that can help you decide. 

First it's worth considering *when to transform*.
<!--toDo:rl-->
<!--not longer valid-->
In some cases transformation to a projected CRS is essential, such as when using geometric functions such as `st_buffer()`, as Figure \@ref(fig:crs-buf) showed.
Conversely, publishing data online with the **leaflet** package may require a geographic CRS.
Another case is when two objects with different CRSs must be compared or combined, as shown when we try to find the distance between two objects with different CRSs:


```r
st_distance(london_geo, london_proj)
# > Error: st_crs(x) == st_crs(y) is not TRUE
```

To make the `london` and `london_proj` objects geographically comparable one of them must be transformed into the CRS of the other.
But which CRS to use?
The answer depends on context: many projects, especially those involving web mapping, require outputs in EPSG:4326, in which case it is worth transforming the projected object.
If, however, the project requires planar geometry operations rather than spherical geometry operations engine (e.g. to create buffers with smooth edges), it may be worth transforming data with a geographic CRS into an equivalent object with a projected CRS, such as the British National Grid (EPSG:27700).
That is the subject of Section \@ref(reproj-vec-geom).

## Which CRS to use? {#which-crs}

\index{CRS!reprojection} 
\index{projection!World Geodetic System}
The question of *which CRS* is tricky, and there is rarely a 'right' answer:
"There exist no all-purpose projections, all involve distortion when far from the center of the specified frame" [@bivand_applied_2013].
Additionally, you should not be attached just to one projection for every task.
It is possible to use one projection for some part of the analysis, another projection for a different part, and even some other for visualization.
Always try to pick the CRS that serves your goal best!

When selecting **geographic CRSs**, the answer is often [WGS84](https://en.wikipedia.org/wiki/World_Geodetic_System#A_new_World_Geodetic_System:_WGS_84).
It is used not only for web mapping, but also because GPS datasets and thousands of raster and vector datasets are provided in this CRS by default.
WGS84 is the most common CRS in the world, so it is worth knowing its EPSG code: 4326.
This 'magic number' can be used to convert objects with unusual projected CRSs into something that is widely understood.

What about when a **projected CRS** is required?
In some cases, it is not something that we are free to decide:
"often the choice of projection is made by a public mapping agency" [@bivand_applied_2013].
This means that when working with local data sources, it is likely preferable to work with the CRS in which the data was provided, to ensure compatibility, even if the official CRS is not the most accurate.
The example of London was easy to answer because (a) the British National Grid (with its associated EPSG code 27700) is well known and (b) the original dataset (`london`) already had that CRS.

\index{UTM} 
A commonly used default is Universal Transverse Mercator ([UTM](https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system)), a set of CRSs that divides the Earth into 60 longitudinal wedges and 20 latitudinal segments.
The transverse Mercator projection used by UTM CRSs is conformal but distorts areas and distances with increasing severity with distance from the center of the UTM zone.
Documentation from the GIS software Manifold therefore suggests restricting the longitudinal extent of projects using UTM zones to 6 degrees from the central meridian (source: [manifold.net](http://www.manifold.net/doc/mfd9/universal_transverse_mercator_projection.htm)).
Therefore, we recommend using UTM only when your focus is on preserving angles for relatively small area!

Almost every place on Earth has a UTM code, such as "60H" which refers to northern New Zealand where R was invented.
UTM EPSG codes run sequentially from 32601 to 32660 for northern hemisphere locations and from 32701 to 32760 for southern hemisphere locations.



To show how the system works, let's create a function, `lonlat2UTM()` to calculate the EPSG code associated with any point on the planet as [follows](https://stackoverflow.com/a/9188972/): 


```r
lonlat2UTM = function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {
    utm + 32600
  } else{
    utm + 32700
  }
}
```

The following command uses this function to identify the UTM zone and associated EPSG code for Auckland and London:




```r
lonlat2UTM(c(174.7, -36.9))
#> [1] 32760
lonlat2UTM(st_coordinates(london))
#> [1] 32630
```

Currently, we also have tools helping us to select a proper CRS, which includes the **crssuggest** package <!--add ref or docs-->.
The main function in this package, `suggest_crs()`, takes a spatial object with geographic CRS and returns a list of possible projected CRSs that could be used for the given area.^[This package also allows to figure out the true CRS of the data without any CRS information attached.]
Important note: while this package is helpful in many situations, you need to be aware of the properties of the recommended CRS before you apply it.

\index{CRS!custom} 
In cases where an appropriate CRS is not immediately clear, the choice of CRS should depend on the properties that are most important to preserve in the subsequent maps and analysis.
All CRSs are either equal-area, equidistant, conformal (with shapes remaining unchanged), or some combination of compromises of those (section \@ref(projected-coordinate-reference-systems)).
Custom CRSs with local parameters can be created for a region of interest and multiple CRSs can be used in projects when no single CRS suits all tasks.
'Geodesic calculations' can provide a fall-back if no CRSs are appropriate (see [proj.org/geodesic.html](https://proj.org/geodesic.html)).
Regardless of the projected CRS used, the results may not be accurate for geometries covering hundreds of kilometers.

\index{CRS!custom}
When deciding on a custom CRS, we recommend the following:^[
<!--toDo:rl-->
<!-- jn:I we can assume who is the "anonymous reviewer", can we ask him/her to use his/her name? -->
Many thanks to an anonymous reviewer whose comments formed the basis of this advice.
]

\index{projection!Lambert azimuthal equal-area}
\index{projection!Azimuthal equidistant}
\index{projection!Lambert conformal conic}
\index{projection!Stereographic}
\index{projection!Universal Transverse Mercator}

- A Lambert azimuthal equal-area ([LAEA](https://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection)) projection for a custom local projection (set latitude and longitude of origin to the center of the study area), which is an equal-area projection at all locations but distorts shapes beyond thousands of kilometers
- Azimuthal equidistant ([AEQD](https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection)) projections for a specifically accurate straight-line distance between a point and the center point of the local projection
- Lambert conformal conic ([LCC](https://en.wikipedia.org/wiki/Lambert_conformal_conic_projection)) projections for regions covering thousands of kilometers, with the cone set to keep distance and area properties reasonable between the secant lines
- Stereographic ([STERE](https://en.wikipedia.org/wiki/Stereographic_projection)) projections for polar regions, but taking care not to rely on area and distance calculations thousands of kilometers from the center

One possible approach to automatically select a projected CRS specific to a local dataset is to create an azimuthal equidistant ([AEQD](https://en.wikipedia.org/wiki/Azimuthal_equidistant_projection)) projection for the center-point of the study area.
This involves creating a custom CRS (with no EPSG code) with units of meters based on the center point of a dataset.
This approach should be used with caution: no other datasets will be compatible with the custom CRS created and results may not be accurate when used on extensive datasets covering hundreds of kilometers.

The principles outlined in this section apply equally to vector and raster datasets.
Some features of CRS transformation however are unique to each geographic data model.
We will cover the particularities of vector data transformation in Section \@ref(reproj-vec-geom) and those of raster transformation in Section \@ref(reproj-ras).
<!--toDo:jn-->
<!-- add reference to the modifying section -->

## Reprojecting vector geometries {#reproj-vec-geom}

<!--toDo:rl-->
<!--jn: idea adding info about custom piplines?-->

\index{CRS!reprojection} 
\index{vector!reprojection} 
Chapter \@ref(spatial-class) demonstrated how vector geometries are made-up of points, and how points form the basis of more complex objects such as lines and polygons.
Reprojecting vectors thus consists of transforming the coordinates of these points, which form the vertices of lines and polygons.

Section \@ref(whenproject) contains an example in which at least one `sf` object must be transformed into an equivalent object with a different CRS to calculate the distance between two objects.



```r
london2 = st_transform(london_geo, "EPSG:27700")
```

Now that a transformed version of `london` has been created, using the **sf** function `st_transform()`, the distance between the two representations of London can be found.
It may come as a surprise that `london` and `london2` are just over 2 km apart!^[
The difference in location between the two points is not due to imperfections in the transforming operation (which is in fact very accurate) but the low precision of the manually-created coordinates that created `london` and `london_proj`.
Also surprising may be that the result is provided in a matrix with units of meters.
This is because `st_distance()` can provide distances between many features and because the CRS has units of meters.
Use `as.numeric()` to coerce the result into a regular number.
]


```r
st_distance(london2, london_proj)
#> Units: [m]
#>      [,1]
#> [1,] 2018
```


This is demonstrated below with reference to `cycle_hire_osm`, an `sf` object from **spData** that represents 'docking stations' where you can hire bicycles in London.
The CRS of `sf` objects can be queried --- and as we learned in Section \@ref(reproj-intro) set --- with the function `st_crs()`.
The output is printed as multiple lines of text containing information about the coordinate system:


```r
st_crs(cycle_hire_osm)
#> Coordinate Reference System:
#>   User input: EPSG:4326 
#>   wkt:
#> GEOGCS["WGS 84",
#>     DATUM["WGS_1984",
#>         SPHEROID["WGS 84",6378137,298.257223563,
#>             AUTHORITY["EPSG","7030"]],
#>         AUTHORITY["EPSG","6326"]],
#>     PRIMEM["Greenwich",0,
#>         AUTHORITY["EPSG","8901"]],
#>     UNIT["degree",0.0174532925199433,
#>         AUTHORITY["EPSG","9122"]],
#>     AUTHORITY["EPSG","4326"]]
```

As we saw in Section \@ref(crs-setting), the main CRS components, `User input` and `wkt`, are printed as a single entity, the output of `st_crs()` is in fact a named list of class `crs` with two elements, single character strings named `input` and `wkt`:


```
#> [1] "crs"
#> [1] "input" "wkt"
```

Additional elements can be retrieved with the `$` operator, including `Name`, `proj4string` and `epsg` (see [`?st_crs`](https://r-spatial.github.io/sf/reference/st_crs.html) and the CRS and tranformation tutorial on the GDAL [website](https://gdal.org/tutorials/osr_api_tut.html#querying-coordinate-reference-system) for details):


```r
crs_lnd$Name
#> [1] "WGS 84"
crs_lnd$proj4string
#> [1] "+proj=longlat +datum=WGS84 +no_defs"
crs_lnd$epsg
#> [1] 4326
```


<!--toDo:rl-->
<!--not longer valid-->
<!-- This duality of CRS objects means that they can be set either using an EPSG code or a `proj4string`. -->
<!-- This means that `st_crs("+proj=longlat +datum=WGS84 +no_defs")` is equivalent to `st_crs(4326)`, although not all `proj4string`s have an associated EPSG code. -->
<!-- Both elements of the CRS are changed by transforming the object to a projected CRS: -->



<!--toDo:rl-->
<!--not longer valid-->
<!-- The resulting object has a new CRS with an EPSG code 27700. -->
<!-- But how to find out more details about this EPSG code, or any code? -->
<!-- One option is to search for it online. -->
<!-- Another option is to use a function from the **rgdal** package to find the name of the CRS: -->



<!--toDo:rl-->
<!--not longer valid-->
<!-- The result shows that the EPSG code 27700 represents the British National Grid, a result that could have been found by searching online for "[EPSG 27700](https://www.google.com/search?q=CRS+27700)". -->
<!-- But what about the `proj4string` element? -->
<!-- `proj4string`s are text strings that describe the CRS. -->
<!-- They can be seen as formulas for converting a projected point into a point on the surface of the Earth and can be accessed from `crs` objects as follows (see [proj.org/](https://proj.org/) for further details of what the output means): -->



\BeginKnitrBlock{rmdnote}<div class="rmdnote">Printing a spatial object in the console automatically returns its coordinate reference system.
To access and modify it explicitly, use the `st_crs` function, for example, `st_crs(cycle_hire_osm)`.</div>\EndKnitrBlock{rmdnote}


## Reprojecting raster geometries {#reproj-ras}

\index{raster!reprojection} 
\index{raster!warping} 
\index{raster!transformation} 
\index{raster!resampling} 
The projection concepts described in the previous section apply equally to rasters.
However, there are important differences in reprojection of vectors and rasters:
transforming a vector object involves changing the coordinates of every vertex but this does not apply to raster data.
Rasters are composed of rectangular cells of the same size (expressed by map units, such as degrees or meters), so it is usually impracticable to transform coordinates of pixels separately.
Raster reprojection involves creating a new raster object, often with a different number of columns and rows than the original.
The attributes must subsequently be re-estimated, allowing the new pixels to be 'filled' with appropriate values.
In other words, raster reprojection can be thought of as two separate spatial operations: a vector reprojection of the raster extent to another CRS (Section \@ref(reproj-vec-geom)), and computation of new pixel values through resampling (Section \@ref(resampling)).
Thus in most cases when both raster and vector data are used, it is better to avoid reprojecting rasters and reproject vectors instead.

\BeginKnitrBlock{rmdnote}<div class="rmdnote">Reprojection of the regular rasters is also known as warping. 
Additionally, there is a second similar operation called "transformation".
Instead of resampling all of the values, it leaves all values intact but recomputes new coordinates for every raster cell, changing the grid geometry.
For example, it could convert the input raster (a regular grid) into a curvilinear grid.
The transformation operation can be performed in R using [the **stars** package](https://r-spatial.github.io/stars/articles/stars5.html).</div>\EndKnitrBlock{rmdnote}



The raster reprojection process is done with `project()` from the **terra** package.
Like the `st_transform()` function demonstrated in the previous section, `project()` takes a geographic object (a raster dataset in this case) and some CRS representation as the second argument.
On a side note -- the second argument can also be an existing raster object with a different CRS.

Let's take a look at two examples of raster transformation: using categorical and continuous data.
Land cover data are usually represented by categorical maps.
The `nlcd.tif` file provides information for a small area in Utah, USA obtained from [National Land Cover Database 2011](https://www.mrlc.gov/data/nlcd-2011-land-cover-conus) in the NAD83 / UTM zone 12N CRS.


```r
cat_raster = rast(system.file("raster/nlcd.tif", package = "spDataLarge"))
crs(cat_raster)
#> [1] "PROJCRS[\"NAD83 / UTM zone 12N\",\n    BASEGEOGCRS[\"NAD83\",\n        DATUM[\"North American Datum 1983\",\n            ELLIPSOID[\"GRS 1980\",6378137,298.257222101,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4269]],\n    CONVERSION[\"UTM zone 12N\",\n        METHOD[\"Transverse Mercator\",\n            ID[\"EPSG\",9807]],\n        PARAMETER[\"Latitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",-111,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"Scale factor at natural origin\",0.9996,\n            SCALEUNIT[\"unity\",1],\n            ID[\"EPSG\",8805]],\n        PARAMETER[\"False easting\",500000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]],\n    USAGE[\n        SCOPE[\"Engineering survey, topographic mapping.\"],\n        AREA[\"North America - between 114°W and 108°W - onshore and offshore. Canada - Alberta; Northwest Territories; Nunavut; Saskatchewan.  United States (USA) - Arizona; Colorado; Idaho; Montana; New Mexico; Utah; Wyoming.\"],\n        BBOX[31.33,-114,84,-108]],\n    ID[\"EPSG\",26912]]"
```

In this region, 8 land cover classes were distinguished (a full list of NLCD2011 land cover classes can be found at [mrlc.gov](https://www.mrlc.gov/data/legends/national-land-cover-database-2011-nlcd2011-legend)):


```r
unique(cat_raster)
#>       levels
#> 1      Water
#> 2  Developed
#> 3     Barren
#> 4     Forest
#> 5  Shrubland
#> 6 Herbaceous
#> 7 Cultivated
#> 8   Wetlands
```

When reprojecting categorical rasters, the estimated values must be the same as those of the original.
This could be done using the nearest neighbor method (`near`), which sets each new cell value to the value of the nearest cell (center) of the input raster.
An example is reprojecting `cat_raster` to WGS84, a geographic CRS well suited for web mapping.
The first step is to obtain the PROJ definition of this CRS, which can be done, for example using the [http://spatialreference.org](http://spatialreference.org/ref/epsg/wgs-84/) webpage. 
The final step is to reproject the raster with the `project()` function which, in the case of categorical data, uses the nearest neighbor method (`near`):


```r
cat_raster_wgs84 = project(cat_raster, "EPSG:4326", method = "near")
```

Many properties of the new object differ from the previous one, including the number of columns and rows (and therefore number of cells), resolution (transformed from meters into degrees), and extent, as illustrated in Table \@ref(tab:catraster) (note that the number of categories increases from 8 to 9 because of the addition of `NA` values, not because a new category has been created --- the land cover classes are preserved).


Table: (\#tab:catraster)Key attributes in the original ('cat\_raster') and projected ('cat\_raster\_wgs84') categorical raster datasets.

|CRS   | nrow| ncol|   ncell| resolution| unique_categories|
|:-----|----:|----:|-------:|----------:|-----------------:|
|NAD83 | 1359| 1073| 1458207|    31.5275|                 8|
|WGS84 | 1246| 1244| 1550024|     0.0003|                 9|

Reprojecting numeric rasters (with `numeric` or in this case `integer` values) follows an almost identical procedure.
This is demonstrated below with `srtm.tif` in **spDataLarge** from [the Shuttle Radar Topography Mission (SRTM)](https://www2.jpl.nasa.gov/srtm/), which represents height in meters above sea level (elevation) with the WGS84 CRS:


```r
con_raster = rast(system.file("raster/srtm.tif", package = "spDataLarge"))
crs(con_raster)
#> [1] "GEOGCRS[\"WGS 84\",\n    DATUM[\"World Geodetic System 1984\",\n        ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n            LENGTHUNIT[\"metre\",1]]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433]],\n    CS[ellipsoidal,2],\n        AXIS[\"geodetic latitude (Lat)\",north,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        AXIS[\"geodetic longitude (Lon)\",east,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n    ID[\"EPSG\",4326]]"
```

We will reproject this dataset into a projected CRS, but *not* with the nearest neighbor method which is appropriate for categorical data.
Instead, we will use the bilinear method which computes the output cell value based on the four nearest cells in the original raster.^[Other methods mentioned in Section \@ref(resampling) also can be used here.]
The values in the projected dataset are the distance-weighted average of the values from these four cells:
the closer the input cell is to the center of the output cell, the greater its weight.
The following commands create a text string representing WGS 84 / UTM zone 12N, and reproject the raster into this CRS, using the `bilinear` method:


```r
con_raster_ea = project(con_raster, "EPSG:32612", method = "bilinear")
crs(con_raster_ea)
#> [1] "PROJCRS[\"WGS 84 / UTM zone 12N\",\n    BASEGEOGCRS[\"WGS 84\",\n        DATUM[\"World Geodetic System 1984\",\n            ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n                LENGTHUNIT[\"metre\",1]]],\n        PRIMEM[\"Greenwich\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        ID[\"EPSG\",4326]],\n    CONVERSION[\"UTM zone 12N\",\n        METHOD[\"Transverse Mercator\",\n            ID[\"EPSG\",9807]],\n        PARAMETER[\"Latitude of natural origin\",0,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8801]],\n        PARAMETER[\"Longitude of natural origin\",-111,\n            ANGLEUNIT[\"degree\",0.0174532925199433],\n            ID[\"EPSG\",8802]],\n        PARAMETER[\"Scale factor at natural origin\",0.9996,\n            SCALEUNIT[\"unity\",1],\n            ID[\"EPSG\",8805]],\n        PARAMETER[\"False easting\",500000,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8806]],\n        PARAMETER[\"False northing\",0,\n            LENGTHUNIT[\"metre\",1],\n            ID[\"EPSG\",8807]]],\n    CS[Cartesian,2],\n        AXIS[\"(E)\",east,\n            ORDER[1],\n            LENGTHUNIT[\"metre\",1]],\n        AXIS[\"(N)\",north,\n            ORDER[2],\n            LENGTHUNIT[\"metre\",1]],\n    USAGE[\n        SCOPE[\"Engineering survey, topographic mapping.\"],\n        AREA[\"Between 114°W and 108°W, northern hemisphere between equator and 84°N, onshore and offshore. Canada - Alberta; Northwest Territories (NWT); Nunavut; Saskatchewan. Mexico. United States (USA).\"],\n        BBOX[0,-114,84,-108]],\n    ID[\"EPSG\",32612]]"
```

Raster reprojection on numeric variables also leads to small changes to values and spatial properties, such as the number of cells, resolution, and extent.
These changes are demonstrated in Table \@ref(tab:rastercrs)^[
Another minor change, that is not represented in Table \@ref(tab:rastercrs), is that the class of the values in the new projected raster dataset is `numeric`.
This is because the `bilinear` method works with continuous data and the results are rarely coerced into whole integer values.
This can have implications for file sizes when raster datasets are saved.
]:


Table: (\#tab:rastercrs)Key attributes in the original ('con\_raster') and projected ('con\_raster\_ea') continuous raster datasets.

|CRS          | nrow| ncol|  ncell| resolution| mean|
|:------------|----:|----:|------:|----------:|----:|
|WGS84        |  457|  465| 212505|     0.0008| 1843|
|UTM zone 12N |  515|  422| 217330|    83.5334| 1842|

\BeginKnitrBlock{rmdnote}<div class="rmdnote">Of course, the limitations of 2D Earth projections apply as much to vector as to raster data.
At best we can comply with two out of three spatial properties (distance, area, direction).
Therefore, the task at hand determines which projection to choose. 
For instance, if we are interested in a density (points per grid cell or inhabitants per grid cell) we should use an equal-area projection (see also Chapter \@ref(location)).</div>\EndKnitrBlock{rmdnote}

## Custom map projections {#mapproj}

<!--     Custom CRSs are also ideally specified as WKT2 -->
<!--     https://epsg.io/ -->
<!-- the two below websites are not up-to-date -->
<!--     https://spatialreference.org/ref/epsg/ -->
<!--     https://epsg.org/home.html -->

<!-- https://projectionwizard.org/ -->

<!--toDo:jn-->
<!--not longer valid-->
<!-- proj4strings still can be used - explain how to use them and when -->
<!-- however, focus on wkt2 customization here! -->
<!-- also, consider moving this section to the bottom of the chapter and show some raster examples -->

<!-- \index{CRS!proj4string}  -->
<!-- Established CRSs captured by EPSG codes are well-suited for many applications. -->
<!-- However in some cases it is desirable to create a new CRS, using a custom `proj4string`. -->
<!-- This system allows a very wide range of projections to be created, as we'll see in some of the custom map projections in this section. -->

<!-- A long and growing list of projections has been developed and many of these can be set with the `+proj=` element of `proj4string`s.^[ -->
<!-- The Wikipedia page 'List of map projections' has 70+ projections and illustrations. -->
<!-- ] -->

<!-- When mapping the world while preserving area relationships, the Mollweide projection is a good choice [@jenny_guide_2017] (Figure \@ref(fig:mollproj)). -->
<!-- To use this projection, we need to specify it using the `proj4string` element, `"+proj=moll"`, in the `st_transform` function: -->





On the other hand, when mapping the world, it is often desirable to have as little distortion as possible for all spatial properties (area, direction, distance).
One of the most popular projections to achieve this is the Winkel tripel projection (Figure \@ref(fig:wintriproj)).^[
This projection is used, among others, by the National Geographic Society.
]
`st_transform_proj()` from the **lwgeom** package allows for coordinate transformations to a wide range of CRSs, including the Winkel tripel projection:


```r
world_wintri = lwgeom::st_transform_proj(world, crs = "+proj=wintri")
```

<div class="figure" style="text-align: center">
<img src="07-reproj_files/figure-html/wintriproj-1.png" alt="Winkel tripel projection of the world." width="100%" />
<p class="caption">(\#fig:wintriproj)Winkel tripel projection of the world.</p>
</div>





<!-- Moreover, PROJ parameters can be modified in most CRS definitions. -->
<!-- The below code transforms the coordinates to the Lambert azimuthal equal-area projection centered on longitude and latitude of `0` (Figure \@ref(fig:laeaproj1)). -->


<!-- plot(world_laea1$geom) -->
<!-- plot(world_laea1$geom, graticule = TRUE) -->



<!-- We can change the PROJ parameters, for example the center of the projection, using the `+lon_0` and `+lat_0` parameters.  -->
<!-- The code below gives the map centered on New York City (Figure \@ref(fig:laeaproj2)). -->





More information on CRS modifications can be found in the [Using PROJ](https://proj.org/usage/index.html) documentation.

There is more to learn about CRSs.
An excellent resource in this area, also implemented in R, is the website R Spatial.
Chapter 6 from this free online book is recommended reading --- see: [rspatial.org/terra/spatial/6-crs.html](https://rspatial.org/terra/spatial/6-crs.html)

## Exercises


E1. Create a new object called `nz_wgs` by transforming `nz` object into the WGS84 CRS.

- Create an object of class `crs` for both and use this to query their CRSs.
- With reference to the bounding box of each object, what units does each CRS use?
- Remove the CRS from `nz_wgs` and plot the result: what is wrong with this map of New Zealand and why?



E2. Transform the `world` dataset to the transverse Mercator projection (`"+proj=tmerc"`) and plot the result.
What has changed and why?
Try to transform it back into WGS 84 and plot the new object.
Why does the new object differ from the original one?



E3. Transform the continuous raster (`con_raster`) into NAD83 / UTM zone 12N using the nearest neighbor interpolation method.
What has changed?
How does it influence the results?



E4. Transform the categorical raster (`cat_raster`) into WGS 84 using the bilinear interpolation method.
What has changed?
How does it influence the results?



<!--toDo:jn-->
<!--improve/replace/modify the following q-->
<!-- E5. Create your own `proj4string`.  -->
<!-- It should have the Lambert Azimuthal Equal Area (`laea`) projection, the WGS84 ellipsoid, the longitude of projection center of 95 degrees west, the latitude of projection center of 60 degrees north, and its units should be in meters. -->
<!-- Next, subset Canada from the `world` object and transform it into the new projection.  -->
<!-- Plot and compare a map before and after the transformation. -->

<!-- ```{r 06-reproj-40} -->
<!-- new_p4s = "+proj=laea +ellps=WGS84 +lon_0=-95 +lat_0=60 +units=m" -->
<!-- canada = dplyr::filter(world, name_long == "Canada") -->
<!-- new_canada = st_transform(canada, new_p4s) -->
<!-- par(mfrow = c(1, 2)) -->
<!-- plot(st_geometry(canada), graticule = TRUE, axes = TRUE) -->
<!-- plot(st_geometry(new_canada), graticule = TRUE, axes = TRUE) -->
<!-- ``` -->
