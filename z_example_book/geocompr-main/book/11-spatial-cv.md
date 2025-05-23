# Statistical learning {#spatial-cv}

## Prerequisites {-}

This chapter assumes proficiency with geographic data analysis\index{geographic data analysis}, for example gained by studying the contents and working-through the exercises in Chapters \@ref(spatial-class) to \@ref(reproj-geo-data).
A familiarity with generalized linear models (GLM)\index{GLM} and machine learning\index{machine learning} is highly recommended [for example from @zuur_mixed_2009;@james_introduction_2013].

The chapter uses the following packages:^[
Package **kernlab**, **pROC**, **RSAGA**\index{RSAGA (package)} and **spDataLarge** must also be installed although these do not need to be attached.
]


```r
library(sf)
library(raster)
library(dplyr)
```


```r
library(mlr)
library(parallelMap)
```


Required data will be attached in due course.

## Introduction {#intro-cv1}

Statistical learning\index{statistical learning} is concerned with the use of statistical and computational models for identifying patterns in data and predicting from these patterns.
Due to its origins, statistical learning\index{statistical learning} is one of R's\index{R} great strengths (see Section \@ref(software-for-geocomputation)).^[
Applying statistical techniques to geographic data has been an active topic of research for many decades in the fields of Geostatistics, Spatial Statistics and point pattern analysis [@diggle_modelbased_2007; @gelfand_handbook_2010; @baddeley_spatial_2015].
]
Statistical learning\index{statistical learning} combines methods from statistics\index{statistics} and machine learning\index{machine learning} and its methods can be categorized into supervised and unsupervised techniques.
Both are increasingly used in disciplines ranging from physics, biology and ecology to geography and economics [@james_introduction_2013].

This chapter focuses on supervised techniques in which there is a training dataset, as opposed to unsupervised techniques such as clustering\index{clustering}.
Response variables can be binary (such as landslide occurrence), categorical (land use), integer (species richness count) or numeric (soil acidity measured in pH).
Supervised techniques model the relationship between such responses --- which are known for a sample of observations --- and one or more predictors.

<!-- For this we can use techniques from the field of statistics or from the field of machine learning.
Which to use depends on the primary aim: statistical inference or prediction.
Statistical regression techniques are especially useful if the aim is statistical inference.
These techniques also allow predictions of unseen data points but this is usually only of secondary interest to statisticians.
Statistical inference, on the other hand, refers among others to a predictor's significance, its importance for a specific model, its relationship with the response and the uncertainties associated with the estimated coefficients.
To trust the p-values and standard errors of such models we need to perform a thorough model validation testing if one or several of the underlying model assumptions (heterogeneity, independence, etc.) have been violated [@zuur_mixed_2009].
By contrast, statistical inference is impossible with machine learning [@james_introduction_2013].
-->
<!-- The primary aim of machine learning is to make good predictions, whereas the field of statistics is more focussed on the underlying theory [e.g. @zuur_mixed_2009] -->

The primary aim of much machine learning\index{machine learning} research is to make good predictions, as opposed to statistical/Bayesian inference, which is good at helping to understand underlying mechanisms and uncertainties in the data [see @krainski_advanced_2018].
Machine learning thrives in the age of 'big data'\index{big data} because its methods make few assumptions about input variables and can handle huge datasets.
Machine learning is conducive to tasks such as the prediction of future customer behavior, recommendation services (music, movies, what to buy next), face recognition, autonomous driving, text classification and predictive maintenance (infrastructure, industry).

<!-- ^[In this case we do not have too worry too much about possible model misspecifications since we explicitly do not want to do statistical inference.] -->

This chapter is based on a case study: the (spatial) prediction of landslides.
This application links to the applied nature of geocomputation, defined in Chapter \@ref(intro), and illustrates how machine learning\index{machine learning} borrows from the field of statistics\index{statistics} when the sole aim is prediction.
Therefore, this chapter first introduces modeling and cross-validation\index{cross-validation} concepts with the help of a Generalized Linear Model \index{GLM} [@zuur_mixed_2009].
Building on this, the chapter implements a more typical machine learning\index{machine learning} algorithm\index{algorithm}, namely a Support Vector Machine (SVM)\index{SVM}.
The models' **predictive performance** will be assessed using spatial cross-validation (CV)\index{cross-validation!spatial CV}, which accounts for the fact that geographic data is special.

CV\index{cross-validation} determines a model's ability to generalize to new data, by splitting a dataset (repeatedly) into training and test sets.
It uses the training data to fit the model, and checks its performance when predicting against the test data.
CV helps to detect overfitting\index{overfitting} since models that predict the training data too closely (noise) will tend to perform poorly on the test data.

Randomly splitting spatial data can lead to training points that are neighbors in space with test points.
Due to spatial autocorrelation\index{autocorrelation!spatial}, test and training datasets would not be independent in this scenario, with the consequence that CV\index{cross-validation} fails to detect a possible overfitting\index{overfitting}.
Spatial CV\index{cross-validation!spatial CV} alleviates this problem and is the **central** theme in this chapter.

## Case study: Landslide susceptibility {#case-landslide}

This case study is based on a dataset of landslide locations in Southern Ecuador, illustrated in Figure \@ref(fig:lsl-map) and described in detail in @muenchow_geomorphic_2012.
A subset of the dataset used in that paper is provided in the **RSAGA**\index{RSAGA (package)} package, which can be loaded as follows:


```r
data("landslides", package = "RSAGA")
```

This should load three objects: a `data.frame` named `landslides`, a `list` named `dem`, and an `sf` object named `study_area`.
`landslides` contains a factor column `lslpts` where `TRUE` corresponds to an observed landslide 'initiation point', with the coordinates stored in columns `x` and `y`.^[
The landslide initiation point is located in the scarp of a landslide polygon. See @muenchow_geomorphic_2012 for further details.
]

There are 175 landslide points and 1360 non-landslide, as shown by `summary(landslides)`.
The 1360 non-landslide points were sampled randomly from the study area, with the restriction that they must fall outside a small buffer around the landslide polygons.

To make the number of landslide and non-landslide points balanced, let us sample 175 from the 1360 non-landslide points.^[The `landslides` dataset has been used in classes and summer schools.
To show how predictive performance of different algorithms changes with an unbalanced and highly spatially autocorrelated response variable, 1360 non-landslide points were randomly selected, i.e., many more absences than presences.
However, especially a logistic regression\index{regression!logistic} with a log-link, as used in this chapter, expects roughly the same number of presences and absences in the response.]


```r
# select non-landslide points
non_pts = filter(landslides, lslpts == FALSE)
# select landslide points
lsl_pts = filter(landslides, lslpts == TRUE)
# randomly select 175 non-landslide points
set.seed(11042018)
non_pts_sub = sample_n(non_pts, size = nrow(lsl_pts))
# create smaller landslide dataset (lsl)
lsl = bind_rows(non_pts_sub, lsl_pts)
```

`dem` is a digital elevation model\index{digital elevation model} consisting of two elements:
`dem$header`, a `list` which represents a raster 'header'\index{raster!header} (see Section \@ref(raster-data)), and `dem$data`, a matrix with the altitude of each pixel.
`dem` can be converted into a `raster` object with:

<!-- Idea: could create a function to do this -->


```r
dem = raster(
  dem$data, 
  crs = dem$header$proj4string,
  xmn = dem$header$xllcorner, 
  xmx = dem$header$xllcorner + dem$header$ncols * dem$header$cellsize,
  ymn = dem$header$yllcorner,
  ymx = dem$header$yllcorner + dem$header$nrows * dem$header$cellsize
  )
```


<div class="figure" style="text-align: center">
<img src="figures/lsl-map-1.png" alt="Landslide initiation points (red) and points unaffected by landsliding (blue) in Southern Ecuador." width="100%" />
<p class="caption">(\#fig:lsl-map)Landslide initiation points (red) and points unaffected by landsliding (blue) in Southern Ecuador.</p>
</div>
\index{hillshade}

To model landslide susceptibility, we need some predictors.
Terrain attributes are frequently associated with landsliding [@muenchow_geomorphic_2012], and these can be computed from the digital elevation model (`dem`) using R-GIS bridges (see Chapter \@ref(gis)).
We leave it as an exercise to the reader to compute the following terrain attribute rasters and extract the corresponding values to our landslide/non-landslide data frame (see exercises; we also provide the resulting data frame via the **spDataLarge** package, see further below):

- `slope`: slope angle (°).
- `cplan`: plan curvature (rad m^−1^) expressing the convergence or divergence of a slope and thus water flow.
- `cprof`: profile curvature (rad m^-1^) as a measure of flow acceleration, also known as downslope change in slope angle.
- `elev`: elevation (m a.s.l.) as the representation of different altitudinal zones of vegetation and precipitation in the study area.
- `log10_carea`: the decadic logarithm of the catchment area (log10 m^2^) representing the amount of water flowing towards a location.

Data containing the landslide points, with the corresponding terrain attributes, is provided in the **spDataLarge** package, along with the terrain attribute raster stack from which the values were extracted.
Hence, if you have not computed the predictors yourself, attach the corresponding data before running the code of the remaining chapter:


```r
# attach landslide points with terrain attributes
data("lsl", package = "spDataLarge")
# attach terrain attribute raster stack
data("ta", package = "spDataLarge")
```

The first three rows of `lsl`, rounded to two significant digits, can be found in Table \@ref(tab:lslsummary).

<table class="table" style="margin-left: auto; margin-right: auto;">
<caption>(\#tab:lslsummary)Structure of the lsl dataset.</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> x </th>
   <th style="text-align:right;"> y </th>
   <th style="text-align:left;"> lslpts </th>
   <th style="text-align:right;"> slope </th>
   <th style="text-align:right;"> cplan </th>
   <th style="text-align:right;"> cprof </th>
   <th style="text-align:right;"> elev </th>
   <th style="text-align:right;"> log10_carea </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 715078 </td>
   <td style="text-align:right;"> 9558647 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 37 </td>
   <td style="text-align:right;"> 0.021 </td>
   <td style="text-align:right;"> 0.009 </td>
   <td style="text-align:right;"> 2500 </td>
   <td style="text-align:right;"> 2.6 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 713748 </td>
   <td style="text-align:right;"> 9558047 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 42 </td>
   <td style="text-align:right;"> -0.024 </td>
   <td style="text-align:right;"> 0.007 </td>
   <td style="text-align:right;"> 2500 </td>
   <td style="text-align:right;"> 3.1 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 712508 </td>
   <td style="text-align:right;"> 9558887 </td>
   <td style="text-align:left;"> FALSE </td>
   <td style="text-align:right;"> 20 </td>
   <td style="text-align:right;"> 0.039 </td>
   <td style="text-align:right;"> 0.015 </td>
   <td style="text-align:right;"> 2100 </td>
   <td style="text-align:right;"> 2.3 </td>
  </tr>
</tbody>
</table>

## Conventional modeling approach in R {#conventional-model}

Before introducing the **mlr**\index{mlr (package)} package, an umbrella-package providing a unified interface to dozens of learning algorithms (Section \@ref(spatial-cv-with-mlr)), it is worth taking a look at the conventional modeling interface in R\index{R}.
This introduction to supervised statistical learning\index{statistical learning} provides the basis for doing spatial CV\index{cross-validation!spatial CV}, and contributes to a better grasp on the **mlr**\index{mlr (package)} approach presented subsequently.

Supervised learning involves predicting a response variable as a function of predictors (Section \@ref(intro-cv)). 
In R\index{R}, modeling functions are usually specified using formulas (see `?formula` and the detailed [Formulas in R Tutorial](https://www.datacamp.com/community/tutorials/r-formula-tutorial) for details of R formulas).
The following command specifies and runs a generalized linear model\index{GLM}:


```r
fit = glm(lslpts ~ slope + cplan + cprof + elev + log10_carea,
          family = binomial(),
          data = lsl)
```

It is worth understanding each of the three input arguments:

- A formula, which specifies landslide occurrence (`lslpts`) as a function of the predictors
- A family, which specifies the type of model, in this case `binomial` because the response is binary (see `?family`)
- The data frame which contains the response and the predictors

The results of this model can be printed as follows (`summary(fit)` provides a more detailed account of the results):


```r
class(fit)
#> [1] "glm" "lm"
fit
#> 
#> Call:  glm(formula = lslpts ~ slope + cplan + cprof + elev + log10_carea, 
#>     family = binomial(), data = lsl)
#> 
#> Coefficients:
#> (Intercept)        slope        cplan        cprof         elev  log10_carea  
#>    1.97e+00     9.30e-02    -2.57e+01    -1.43e+01     2.41e-05    -2.12e+00  
#> 
#> Degrees of Freedom: 349 Total (i.e. Null);  344 Residual
#> Null Deviance:	    485 
#> Residual Deviance: 361 	AIC: 373
```

The model object `fit`, of class `glm`, contains the coefficients defining the fitted relationship between response and predictors.
It can also be used for prediction.
This is done with the generic `predict()` method, which in this case calls the function `predict.glm()`.
Setting `type` to `response` returns the predicted probabilities (of landslide occurrence) for each observation in `lsl`, as illustrated below (see `?predict.glm`):


```r
pred_glm = predict(object = fit, type = "response")
head(pred_glm)
#>      1      2      3      4      5      6 
#> 0.3327 0.4755 0.0995 0.1480 0.3486 0.6766
```

Spatial predictions can be made by applying the coefficients to the predictor rasters. 
This can be done manually or with `raster::predict()`.
In addition to a model object (`fit`), this function also expects a raster stack with the predictors named as in the model's input data frame (Figure \@ref(fig:lsl-susc)).


```r
# making the prediction
pred = raster::predict(ta, model = fit, type = "response")
```

<div class="figure" style="text-align: center">
<img src="figures/lsl-susc-1.png" alt="Spatial prediction of landslide susceptibility using a GLM." width="100%" />
<p class="caption">(\#fig:lsl-susc)Spatial prediction of landslide susceptibility using a GLM.</p>
</div>

Here, when making predictions we neglect spatial autocorrelation\index{autocorrelation!spatial} since we assume that on average the predictive accuracy remains the same with or without spatial autocorrelation structures.
However, it is possible to include spatial autocorrelation\index{autocorrelation!spatial} structures into models [@zuur_mixed_2009;@blangiardo_spatial_2015;@zuur_beginners_2017] as well as into predictions [kriging approaches,  see, e.g., @goovaerts_geostatistics_1997;@hengl_practical_2007;@bivand_applied_2013].
This is, however, beyond the scope of this book.
<!--
Nevertheless, we give the interested reader some pointers where to look it up:

1. The predictions of regression kriging combines the predictions of a regression with the kriging of the regression's residuals [@bivand_applied_2013]. 
2. One can also add a spatial correlation (dependency) structure to a generalized least squares model  [`nlme::gls()`; @zuur_mixed_2009; @zuur_beginners_2017].  
3. Finally, there are mixed-effect modeling approaches.
Basically, a random effect imposes a dependency structure on the response variable which in turn allows for observations of one class to be more similar to each other than to those of another class [@zuur_mixed_2009]. 
Classes can be, for example, bee hives, owl nests, vegetation transects or an altitudinal stratification.
This mixed modeling approach assumes normal and independent distributed random intercepts.^[Note that for spatial predictions one would usually use the population intercept.]
This can even be extended by using a random intercept that is normal and spatially dependent.
For this, however, you will have to resort most likely to Bayesian modeling approaches since frequentist software tools are rather limited in this respect especially for more complex models [@blangiardo_spatial_2015; @zuur_beginners_2017]. 
-->

Spatial prediction maps are one very important outcome of a model.
Even more important is how good the underlying model is at making them since a prediction map is useless if the model's predictive performance is bad.
The most popular measure to assess the predictive performance of a binomial model is the Area Under the Receiver Operator Characteristic Curve (AUROC)\index{AUROC}.
This is a value between 0.5 and 1.0, with 0.5 indicating a model that is no better than random and 1.0 indicating perfect prediction of the two classes. 
Thus, the higher the AUROC\index{AUROC}, the better the model's predictive power.
The following code chunk computes the AUROC\index{AUROC} value of the model with `roc()`, which takes the response and the predicted values as inputs. 
`auc()` returns the area under the curve.


```r
pROC::auc(pROC::roc(lsl$lslpts, fitted(fit)))
#> Area under the curve: 0.826
```

An AUROC\index{AUROC} value of
<!--  -->
0.83
represents a good fit.
However, this is an overoptimistic estimation since we have computed it on the complete dataset. 
To derive a biased-reduced assessment, we have to use cross-validation\index{cross-validation} and in the case of spatial data should make use of spatial CV\index{cross-validation!spatial CV}.

## Introduction to (spatial) cross-validation {#intro-cv} 

Cross-validation\index{cross-validation} belongs to the family of resampling methods\index{resampling} [@james_introduction_2013].
The basic idea is to split (repeatedly) a dataset into training and test sets whereby the training data is used to fit a model which then is applied to the test set.
Comparing the predicted values with the known response values from the test set (using a performance measure such as the AUROC\index{AUROC} in the binomial case) gives a bias-reduced assessment of the model's capability to generalize the learned relationship to independent data.
For example, a 100-repeated 5-fold cross-validation means to randomly split the data into five partitions (folds) with each fold being used once as a test set (see upper row of Figure \@ref(fig:partitioning)). 
This guarantees that each observation is used once in one of the test sets, and requires the fitting of five models.
Subsequently, this procedure is repeated 100 times.
Of course, the data splitting will differ in each repetition.
<!--if the error is calc. on the fold-level. most often its calc. on the repetition level. maybe worth noting.
talk about this in person
-->
Overall, this sums up to 500 models, whereas the mean performance measure (AUROC\index{AUROC}) of all models is the model's overall predictive power.

However, geographic data is special.
As we will see in Chapter \@ref(transport), the 'first law' of geography states that points close to each other are, generally, more similar than points further away [@miller_tobler_2004].
This means these points are not statistically independent because training and test points in conventional CV\index{cross-validation} are often too close to each other (see first row of Figure \@ref(fig:partitioning)).
'Training' observations near the 'test' observations can provide a kind of 'sneak preview':
information that should be unavailable to the training dataset.
<!-- "folds" only for the repetition split, "partitions" or "subsets" for splitting within a fold
talk about this in person
-->
To alleviate this problem 'spatial partitioning' is used to split the observations into spatially disjointed subsets (using the observations' coordinates in a *k*-means clustering\index{clustering!kmeans}; @brenning_spatial_2012; second row of Figure \@ref(fig:partitioning)).
This partitioning strategy is the **only** difference between spatial and conventional CV.
As a result, spatial CV leads to a bias-reduced assessment of a model's predictive performance, and hence helps to avoid overfitting\index{overfitting}.
<!-- Alex suggested to remove this: 
It is important to note that spatial CV reduces the bias introduced by spatial autocorrelation but does not completely remove it. 
This is because there are still a few points in the test and training data which are still neighbors (@brenning_spatial_2012; see second row of \@ref(fig:partitioning)).
-->

<div class="figure" style="text-align: center">
<img src="figures/13_partitioning.png" alt="Spatial visualization of selected test and training observations for cross-validation of one repetition. Random (upper row) and spatial partitioning (lower row)." width="100%" />
<p class="caption">(\#fig:partitioning)Spatial visualization of selected test and training observations for cross-validation of one repetition. Random (upper row) and spatial partitioning (lower row).</p>
</div>

## Spatial CV with **mlr**
\index{mlr (package)}
There are dozens of packages for statistical learning\index{statistical learning}, as described for example in the [CRAN machine learning task view](https://CRAN.R-project.org/view=MachineLearning).
Getting acquainted with each of these packages, including how to undertake cross-validation and hyperparameter\index{hyperparameter} tuning, can be a time-consuming process.
Comparing model results from different packages can be even more laborious.
The **mlr** package was developed to address these issues.
It acts as a 'meta-package', providing a unified interface to popular supervised and unsupervised statistical learning techniques including classification, regression\index{regression}, survival analysis and clustering\index{clustering} [@bischl_mlr:_2016].
<!-- ^[ -->
<!-- We will solely focus on supervised learning techniques in this chapter. -->
<!-- ] -->
The standardized **mlr** interface is based on eight 'building blocks'.
As illustrated in Figure \@ref(fig:building-blocks), these have a clear order.

<div class="figure" style="text-align: center">
<img src="figures/13_ml_abstraction_crop.png" alt="Basic building blocks of the mlr package. Source: http://bit.ly/2tcb2b7. (Permission to reuse this figure was kindly granted.)" width="100%" />
<p class="caption">(\#fig:building-blocks)Basic building blocks of the mlr package. Source: http://bit.ly/2tcb2b7. (Permission to reuse this figure was kindly granted.)</p>
</div>

The **mlr** modeling process consists of three main stages.
First, a **task** specifies the data (including response and predictor variables) and the model type (such as regression\index{regression} or classification\index{classification}).
Second, a **learner** defines the specific learning algorithm that is applied to the created task.
Third, the **resampling** approach assesses the predictive performance of the model, i.e., its ability to generalize to new data (see also Section \@ref(intro-cv)).

### Generalized linear model {#glm}

To implement a GLM\index{GLM} in **mlr**\index{mlr (package)}, we must create a **task** containing the landslide data.
Since the response is binary (two-category variable), we create a classification\index{classification} task with `makeClassifTask()` (for regression\index{regression} tasks, use `makeRegrTask()`, see `?makeRegrTask` for other task types).
The first essential argument of these `make*()` functions is `data`.
The `target` argument expects the name of a response variable and `positive` determines which of the two factor levels of the response variable indicate the landslide initiation point (in our case this is `TRUE`).
All other variables of the `lsl` dataset will serve as predictors except for the coordinates (see the result of `getTaskFormula(task)` for the model formula).
For spatial CV, the `coordinates` parameter is used (see Section \@ref(intro-cv) and Figure \@ref(fig:partitioning)) which expects the coordinates as a xy data frame.


```r
library(mlr)
# coordinates needed for the spatial partitioning
coords = lsl[, c("x", "y")]
# select response and predictors to use in the modeling
data = dplyr::select(lsl, -x, -y)
# create task
task = makeClassifTask(data = data, target = "lslpts",
                       positive = "TRUE", coordinates = coords)
```

`makeLearner()` determines the statistical learning\index{statistical learning} method to use.
All classification\index{classification} **learners** start with `classif.` and all regression\index{regression} learners with `regr.` (see `?makeLearners` for details). 
`listLearners()` helps to find out about all available learners and from which package **mlr** imports them (Table \@ref(tab:lrns)). 
For a specific task, we can run:


```r
listLearners(task, warn.missing.packages = FALSE) %>%
  dplyr::select(class, name, short.name, package) %>%
  head()
```

<table>
<caption>(\#tab:lrns)Sample of available learners for binomial tasks in the mlr package.</caption>
 <thead>
  <tr>
   <th style="text-align:left;"> Class </th>
   <th style="text-align:left;"> Name </th>
   <th style="text-align:left;"> Short name </th>
   <th style="text-align:left;"> Package </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> classif.adaboostm1 </td>
   <td style="text-align:left;"> ada Boosting M1 </td>
   <td style="text-align:left;"> adaboostm1 </td>
   <td style="text-align:left;"> RWeka </td>
  </tr>
  <tr>
   <td style="text-align:left;"> classif.binomial </td>
   <td style="text-align:left;"> Binomial Regression </td>
   <td style="text-align:left;"> binomial </td>
   <td style="text-align:left;"> stats </td>
  </tr>
  <tr>
   <td style="text-align:left;"> classif.featureless </td>
   <td style="text-align:left;"> Featureless classifier </td>
   <td style="text-align:left;"> featureless </td>
   <td style="text-align:left;"> mlr </td>
  </tr>
  <tr>
   <td style="text-align:left;"> classif.fnn </td>
   <td style="text-align:left;"> Fast k-Nearest Neighbour </td>
   <td style="text-align:left;"> fnn </td>
   <td style="text-align:left;"> FNN </td>
  </tr>
  <tr>
   <td style="text-align:left;"> classif.gausspr </td>
   <td style="text-align:left;"> Gaussian Processes </td>
   <td style="text-align:left;"> gausspr </td>
   <td style="text-align:left;"> kernlab </td>
  </tr>
  <tr>
   <td style="text-align:left;"> classif.IBk </td>
   <td style="text-align:left;"> k-Nearest Neighbours </td>
   <td style="text-align:left;"> ibk </td>
   <td style="text-align:left;"> RWeka </td>
  </tr>
</tbody>
</table>

This yields all learners able to model two-class problems (landslide yes or no).
We opt for the binomial classification\index{classification} method used in Section \@ref(conventional-model) and implemented as `classif.binomial` in **mlr**.
Additionally, we must specify the link-function, `logit` in this case, which is also the default of the `binomial()` function.
`predict.type` determines the type of the prediction with `prob` resulting in the predicted probability for landslide occurrence between 0 and 1 (this corresponds to `type = response` in `predict.glm`).


```r
lrn = makeLearner(cl = "classif.binomial",
                  link = "logit",
                  predict.type = "prob",
                  fix.factors.prediction = TRUE)
```

To find out from which package the specified learner is taken and how to access the corresponding help pages, we can run:


```r
getLearnerPackages(lrn)
helpLearner(lrn)
```

<!--
Having specified a learner and a task, we can train our model which basically executes the `glm()` command in the background for our task. 


```r
mod = train(learner = lrn, task = task)
mlr_fit = getLearnerModel(mod)
```



`getLearnerModel()` extracts the used model which shows that **mlr** passed all specified parameters to the `glm` function in the background as also proved by following code:


```r
fit = glm(lslpts ~ ., family = binomial(link = "logit"), data = data)
identical(fit$coefficients, mlr_fit$coefficients)
```
-->

The set-up steps for modeling with **mlr**\index{mlr (package)} may seem tedious. 
But remember, this single interface provides access to the 150+ learners shown by `listLearners()`; it would be far more tedious to learn the interface for each learner!
Further advantages are simple parallelization of resampling techniques and the ability to tune machine learning hyperparameters\index{hyperparameter} (see Section \@ref(svm)).
Most importantly, (spatial) resampling in **mlr** is straightforward, requiring only two more steps: specifying a resampling method and running it.
We will use a 100-repeated 5-fold spatial CV\index{cross-validation!spatial CV}: five partitions will be chosen based on the provided coordinates in our `task` and the partitioning will be repeated 100 times:[^13]

[^13]: 

    Note that package **sperrorest** initially implemented spatial cross-validation in R [@brenning_spatial_2012].
    In the meantime, its functionality was integrated into the **mlr** package which is the reason why we are using **mlr** [@schratz_performance_nodate].The **caret** package is another umbrella-package [@kuhn_applied_2013] for streamlined modeling in R; however, so far it does not provide spatial CV which is why we refrain from using it for spatial data.



```r
perf_level = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 100)
```

To execute the spatial resampling, we run `resample()` using the specified learner, task, resampling strategy and of course the performance measure, here the AUROC\index{AUROC}.
This takes some time (around 10 seconds on a modern laptop) because it computes the AUROC for 500 models. 
Setting a seed ensures the reproducibility\index{reproducibility} of the obtained result and will ensure the same spatial partitioning when re-running the code.

<!-- I just thought it might be worth showing the differences between an error on the fold level and repetition level but aggregating to the rep level is not a one-liner in mlr -->



```r
set.seed(012348)
sp_cv = mlr::resample(learner = lrn, task = task,
                      resampling = perf_level, 
                      measures = mlr::auc)
```



<!-- sp_cv and conv_cv have been saved in spatialcv.Rdata. I needed to run the modeling outside of the book since knitr sets its own seed and I am not sure if this actually helps to make sure that the same partitions are used in the cv.
I really don't understand why I have to load spatialcv.Rdata here a third time...-->

The output of the preceding code chunk is a bias-reduced assessment of the model's predictive performance, as illustrated in the following code chunk  (required input data is saved in the file `spatialcv.Rdata` in the book's GitHub repo):




```r
# summary statistics of the 500 models
summary(sp_cv$measures.test$auc)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   0.686   0.757   0.789   0.780   0.795   0.861
# mean AUROC of the 500 models
mean(sp_cv$measures.test$auc)
#> [1] 0.78
```

To put these results in perspective, let us compare them with AUROC\index{AUROC} values from a 100-repeated 5-fold non-spatial cross-validation (Figure \@ref(fig:boxplot-cv); the code for the non-spatial cross-validation\index{cross-validation} is not shown here but will be explored in the exercise section).
As expected, the spatially cross-validated result yields lower AUROC values on average than the conventional cross-validation approach, underlining the over-optimistic predictive performance due to spatial autocorrelation\index{autocorrelation!spatial} of the latter.

<div class="figure" style="text-align: center">
<img src="11-spatial-cv_files/figure-html/boxplot-cv-1.png" alt="Boxplot showing the difference in AUROC values between spatial and conventional 100-repeated 5-fold cross-validation." width="100%" />
<p class="caption">(\#fig:boxplot-cv)Boxplot showing the difference in AUROC values between spatial and conventional 100-repeated 5-fold cross-validation.</p>
</div>

### Spatial tuning of machine-learning hyperparameters {#svm}

Section \@ref(intro-cv) introduced machine learning\index{machine learning} as part of statistical learning\index{statistical learning}.
To recap, we adhere to the following definition of machine learning by [Jason Brownlee](https://machinelearningmastery.com/linear-regression-for-machine-learning/):

> Machine learning, more specifically the field of predictive modeling, is primarily concerned with minimizing the error of a model or making the most accurate predictions possible, at the expense of explainability.
In applied machine learning we will borrow, reuse and steal algorithms from many different fields, including statistics and use them towards these ends.

In Section \@ref(glm) a GLM was used to predict landslide susceptibility.
This section introduces support vector machines (SVM)\index{SVM} for the same purpose.
Random forest\index{random forest} models might be more popular than SVMs; however, the positive effect of tuning hyperparameters\index{hyperparameter} on model performance is much more pronounced in the case of SVMs [@probst_hyperparameters_2018].
Since (spatial) hyperparameter tuning is the major aim of this section, we will use an SVM.
For those wishing to apply a random forest model, we recommend to read this chapter, and then proceed to Chapter \@ref(eco) in which we will apply the currently covered concepts and techniques to make spatial predictions based on a random forest model.

SVMs\index{SVM} search for the best possible 'hyperplanes' to separate classes (in a classification\index{classification} case) and estimate 'kernels' with specific hyperparameters to allow for non-linear boundaries between classes [@james_introduction_2013].
Hyperparameters\index{hyperparameter} should not be confused with coefficients of parametric models, which are sometimes also referred to as parameters.^[
For a detailed description of the difference between coefficients and hyperparameters, see the 'machine mastery' blog post on the subject.
<!-- For a more detailed description of the difference between coefficients and hyperparameters, see the [machine mastery blog](https://machinelearningmastery.com/difference-between-a-parameter-and-a-hyperparameter/). -->
]
Coefficients can be estimated from the data, while hyperparameters are set before the learning begins.
Optimal hyperparameters are usually determined within a defined range with the help of cross-validation methods.
This is called hyperparameter tuning.

Some SVM implementations such as that provided by **kernlab** allow hyperparameters to be tuned automatically, usually based on random sampling (see upper row of Figure \@ref(fig:partitioning)).
This works for non-spatial data but is of less use for spatial data where 'spatial tuning' should be undertaken.

Before defining spatial tuning, we will set up the **mlr**\index{mlr (package)} building blocks, introduced in Section \@ref(glm), for the SVM.
The classification\index{classification} task remains the same, hence we can simply reuse the `task` object created in Section \@ref(glm).
Learners implementing SVM can be found using `listLearners()` as follows:


```r
lrns = listLearners(task, warn.missing.packages = FALSE)
filter(lrns, grepl("svm", class)) %>% 
  dplyr::select(class, name, short.name, package)
#>            class                                 name short.name package
#> 6   classif.ksvm              Support Vector Machines       ksvm kernlab
#> 9  classif.lssvm Least Squares Support Vector Machine      lssvm kernlab
#> 17   classif.svm     Support Vector Machines (libsvm)        svm   e1071
```

Of the options illustrated above, we will use `ksvm()` from the **kernlab** package [@karatzoglou_kernlab_2004].
To allow for non-linear relationships, we use the popular radial basis function (or Gaussian) kernel which is also the default of `ksvm()`.


```r
lrn_ksvm = makeLearner("classif.ksvm",
                        predict.type = "prob",
                        kernel = "rbfdot")
```

The next stage is to specify a resampling strategy.
Again we will use a 100-repeated 5-fold spatial CV\index{cross-validation!spatial CV}.

<!-- Instead of saying "outer resampling" we concluded to use "performance estimation level" and "tuning level" (inner) in our paper
# this is also what is shown in the nested CV figure so it would be more consistent -->


```r
# performance estimation level
perf_level = makeResampleDesc(method = "SpRepCV", folds = 5, reps = 100)
```

Note that this is the exact same code as used for the GLM\index{GLM} in Section \@ref(glm); we have simply repeated it here as a reminder.

So far, the process has been identical to that described in Section \@ref(glm).
The next step is new, however: to tune the hyperparameters\index{hyperparameter}.
Using the same data for the performance assessment and the tuning would potentially lead to overoptimistic results [@cawley_overfitting_2010].
This can be avoided using nested spatial CV\index{cross-validation!spatial CV}.

<div class="figure" style="text-align: center">
<img src="figures/13_cv.png" alt="Schematic of hyperparameter tuning and performance estimation levels in CV. (Figure was taken from Schratz et al. (2018). Permission to reuse it  was kindly granted.)" width="100%" />
<p class="caption">(\#fig:inner-outer)Schematic of hyperparameter tuning and performance estimation levels in CV. (Figure was taken from Schratz et al. (2018). Permission to reuse it  was kindly granted.)</p>
</div>

This means that we split each fold again into five spatially disjoint subfolds which are used to determine the optimal hyperparameters\index{hyperparameter} (`tune_level` object in the code chunk below; see Figure \@ref(fig:inner-outer) for a visual representation).
To find the optimal hyperparameter combination, we fit 50 models (`ctrl` object in the code chunk below) in each of these subfolds with randomly selected values for the hyperparameters C and Sigma.
The random selection of values C and Sigma is additionally restricted to a predefined tuning space (`ps` object).
The range of the tuning space was chosen with values recommended in the literature [@schratz_performance_nodate].

<!--
Questions Pat:
- why not using e1071 svm -> inner hyperparameter tuning also possible I guess...
## Because kernlab has more kernel options. Other than that there is no argument
- explanation correct?
## If you mean the paragraph above, yes
- trafo-function?
## is just a different approach of writing the limits. You could also directly write 2^{-15}. Makes it easier to see the limits at the first glance. Personal preference though
- 125,000 models
-->

<!--
talk in person (see also exercises):
- can I compare the mean AUROC of the GLM and the SVM when using the same seed? Or is seeding not strictly necessary? I mean, ok, the partitions vary a bit but overall...
-->


```r
# five spatially disjoint partitions
tune_level = makeResampleDesc("SpCV", iters = 5)
# use 50 randomly selected hyperparameters
ctrl = makeTuneControlRandom(maxit = 50)
# define the outer limits of the randomly selected hyperparameters
ps = makeParamSet(
  makeNumericParam("C", lower = -12, upper = 15, trafo = function(x) 2^x),
  makeNumericParam("sigma", lower = -15, upper = 6, trafo = function(x) 2^x)
  )
```

The next stage is to modify the learner `lrn_ksvm` in accordance with all the characteristics defining the hyperparameter tuning with `makeTuneWrapper()`.


```r
wrapped_lrn_ksvm = makeTuneWrapper(learner = lrn_ksvm, 
                                   resampling = tune_level,
                                   par.set = ps,
                                   control = ctrl, 
                                   show.info = TRUE,
                                   measures = mlr::auc)
```

The **mlr** is now set-up to fit 250 models to determine optimal hyperparameters for one fold.
Repeating this for each fold, we end up with 1250 (250 \* 5) models for each repetition.
Repeated 100 times means fitting a total of 125,000 models to identify optimal hyperparameters (Figure \@ref(fig:partitioning)).
These are used in the performance estimation, which requires the fitting of another 500 models (5 folds \* 100 repetitions; see Figure \@ref(fig:partitioning)). 
To make the performance estimation processing chain even clearer, let us write down the commands we have given to the computer:

1. Performance level (upper left part of Figure \@ref(fig:inner-outer)): split the dataset into five spatially disjoint (outer) subfolds.
1. Tuning level (lower left part of Figure \@ref(fig:inner-outer)): use the first fold of the performance level and split it again spatially into five (inner) subfolds for the hyperparameter tuning. 
Use the 50 randomly selected hyperparameters\index{hyperparameter} in each of these inner subfolds, i.e., fit 250 models.
1. Performance estimation: Use the best hyperparameter combination from the previous step (tuning level) and apply it to the first outer fold in the performance level to estimate the performance (AUROC\index{AUROC}).
1. Repeat steps 2 and 3 for the remaining four outer folds.
1. Repeat steps 2 to 4, 100 times.

The process of hyperparameter tuning and performance estimation is computationally intensive.
Model runtime can be reduced with parallelization, which can be done in a number of ways, depending on the operating system.
<!-- "cloud development is done on linux servers" is somehow a strange read that I cannot relate really. Maybe sth like: "Parallelilaztion and cloud-computing are most often done on Linux operating systems nowadays. This has some reasons, one of them that directly affects us is that only on Linux systems we can set a parallel seed in R that makes the parallel processes reproducible [this is still an assumption, I will check on that!]"-->

Before starting the parallelization\index{parallelization}, we ensure that the processing continues even if one of the models throws an error by setting `on.learner.error` to `warn`.
This avoids the process stopping just because of one failed model, which is desirable on large model runs.
To inspect the failed models once the processing is completed, we dump them:


```r
configureMlr(on.learner.error = "warn", on.error.dump = TRUE)
```

To start the parallelization\index{parallelization}, we set the `mode` to `multicore` which will use `mclapply()` in the background on a single machine in the case of a Unix-based operating system.^[
See `?parallelStart` for further modes and github.com/berndbischl/parallelMap for more on the unified interface to popular parallelization back-ends.
]
Equivalenty, `parallelStartSocket()` enables parallelization under Windows. 
`level` defines the level at which to enable parallelization, with `mlr.tuneParams` determining that the hyperparameter tuning level should be parallelized (see lower left part of Figure \@ref(fig:inner-outer), `?parallelGetRegisteredLevels`, and the **mlr** [parallelization tutorial](https://mlr-org.github.io/mlr-tutorial/release/html/parallelization/index.html#parallelization-levels) for details).
We will use half of the available cores (set with the `cpus` parameter), a setting that allows possible other users to work on the same high performance computing cluster in case one is used (which was the case when we ran the code).
<!-- the partitions are created before the parallelization by the normal set.seed() call. mc.set.seed makes sure that the randomly chosen hyperparameters for the tuning are reproducible. These will first set within the parallelization.-->
Setting `mc.set.seed` to `TRUE` ensures that the randomly chosen hyperparameters\index{hyperparameter} during the tuning can be reproduced when running the code again.
Unfortunately, `mc.set.seed` is only available under Unix-based systems.


```r
library(parallelMap)
if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
parallelStart(mode = "multicore", 
              # parallelize the hyperparameter tuning level
              level = "mlr.tuneParams", 
              # just use half of the available cores
              cpus = round(parallel::detectCores() / 2),
              mc.set.seed = TRUE)
}

if (Sys.info()["sysname"] == "Windows") {
  parallelStartSocket(level = "mlr.tuneParams",
                      cpus =  round(parallel::detectCores() / 2))
}
```

Now we are set up for computing the nested spatial CV.
Using a seed allows us to recreate the exact same spatial partitions when re-running the code.
Specifying the `resample()` parameters follows the exact same procedure as presented when using a GLM\index{GLM}, the only difference being the `extract` argument.
This allows the extraction of the hyperparameter\index{hyperparameter} tuning results which is important if we plan follow-up analyses on the tuning.
After the processing, it is good practice to explicitly stop the parallelization\index{parallelization} with `parallelStop()`.
Finally, we save the output object (`result`) to disk in case we would like to use it another R session.
Before running the subsequent code, be aware that it is time-consuming:
the 125,500 models took ~1/2hr on a server using 24 cores (see below).


```r
set.seed(12345)
result = mlr::resample(learner = wrapped_lrn_ksvm, 
                       task = task,
                       resampling = perf_level,
                       extract = getTuneResult,
                       measures = mlr::auc)
# stop parallelization
parallelStop()
# save your result, e.g.:
# saveRDS(result, "svm_sp_sp_rbf_50it.rds")
```

In case you do not want to run the code locally, we have saved a subset of the [results](https://github.com/Robinlovelace/geocompr/blob/main/extdata/spatial_cv_result.rds) in the book's GitHub repo.
They can be loaded as follows:


```r
result = readRDS("extdata/spatial_cv_result.rds")
```

Note that runtime depends on many aspects: CPU speed, the selected algorithm, the selected number of cores and the dataset.


```r
# Exploring the results
# runtime in minutes
round(result$runtime / 60, 2)
#> [1] 37.4
```

<!-- interestingly, running the code with 4 cores takes 41.67 minutes -->

Even more important than the runtime is the final aggregated AUROC\index{AUROC}: the model's ability to discriminate the two classes. 


```r
# final aggregated AUROC 
result$aggr
#> auc.test.mean 
#>         0.758
# same as
mean(result$measures.test$auc)
#> [1] 0.758
```

It appears that the GLM\index{GLM} (aggregated AUROC\index{AUROC} was 0.78) is slightly better than the SVM\index{SVM} in this specific case.
However, using more than 50 iterations in the random search would probably yield hyperparameters\index{hyperparameter} that result in models with a better AUROC [@schratz_performance_nodate].
On the other hand, increasing the number of random search iterations would also increase the total number of models and thus runtime.

The estimated optimal hyperparameters for each fold at the performance estimation level can also be viewed.
The following command shows the best hyperparameter combination of the first fold of the first iteration (recall this results from the first 5 \* 50 model runs):


```r
# winning hyperparameters of tuning step, 
# i.e. the best combination out of 50 * 5 models
result$extract[[1]]$x
#> $C
#> [1] 0.458
#> 
#> $sigma
#> [1] 0.023
```

The estimated hyperparameters have been used for the first fold in the first iteration of the performance estimation level which resulted in the following AUROC\index{AUROC} value:


```r
result$measures.test[1, ]
#>   iter   auc
#> 1    1 0.799
```

<!-- # maybe add a figure (boxplot) showing the differences between tuning and no tuning?-->

So far spatial CV\index{cross-validation!spatial CV} has been used to assess the ability of learning algorithms to generalize to unseen data.
For spatial predictions, one would tune the hyperparameters\index{hyperparameter} on the complete dataset.
This will be covered in Chapter \@ref(eco).

## Conclusions

Resampling methods are an important part of a data scientist's toolbox [@james_introduction_2013]. 
This chapter used cross-validation\index{cross-validation} to assess predictive performance of various models.
As described in Section \@ref(intro-cv), observations with spatial coordinates may not be statistically independent due to spatial autocorrelation\index{autocorrelation!spatial}, violating a fundamental assumption of cross-validation.
Spatial CV\index{cross-validation!spatial CV} addresses this issue by reducing bias introduced by spatial autocorrelation\index{autocorrelation!spatial}. 

The **mlr**\index{mlr (package)} package facilitates (spatial) resampling\index{resampling} techniques in combination with the most popular statistical learning\index{statistical learning} techniques including linear regression\index{regression!linear}, semi-parametric models such as generalized additive models\index{generalized additive model} and machine learning\index{machine learning} techniques such as random forests\index{random forest}, SVMs\index{SVM}, and boosted regression trees [@bischl_mlr:_2016;@schratz_performance_nodate].
Machine learning algorithms often require hyperparameter\index{hyperparameter} inputs, the optimal 'tuning' of which can require thousands of model runs which require large computational resources, consuming much time, RAM and/or cores.
**mlr** tackles this issue by enabling parallelization\index{parallelization}.

Machine learning overall, and its use to understand spatial data, is a large field and this chapter has provided the basics, but there is more to learn.
We recommend the following resources in this direction:

- The **mlr** tutorials on [Machine Learning in R](https://mlr-org.github.io/mlr-tutorial/release/html/) and [Handling of spatial Data](https://mlr-org.github.io/mlr-tutorial/release/html/handling_of_spatial_data/index.html)
- An academic paper on hyperparameter\index{hyperparameter} tuning [@schratz_performance_nodate]
- In case of spatio-temporal data, one should account for spatial\index{autocorrelation!spatial} and temporal\index{autocorrelation!temporal} autocorrelation when doing CV\index{cross-validation} [@meyer_improving_2018]

## Exercises

1. Compute the following terrain attributes from the `dem` datasets loaded with `data("landslides", package = "RSAGA")` with the help of R-GIS bridges (see Chapter \@ref(gis)):
    - Slope
    - Plan curvature
    - Profile curvature
    - Catchment area
1. Extract the values from the corresponding output rasters to the `landslides` data frame (`data(landslides, package = "RSAGA"`) by adding new variables called `slope`, `cplan`, `cprof`, `elev` and `log_carea`. Keep all landslide initiation points and 175 randomly selected non-landslide points (see Section \@ref(case-landslide) for details).
1. Use the derived terrain attribute rasters in combination with a GLM to make a spatial prediction map similar to that shown in Figure \@ref(fig:lsl-susc).
Running `data("study_mask", package = "spDataLarge")` attaches a mask of the study area.
1. Compute a 100-repeated 5-fold non-spatial cross-validation and spatial CV based on the GLM learner and compare the AUROC values from both resampling\index{resampling} strategies with the help of boxplots (see Figure \@ref(fig:boxplot-cv)).
Hint: You need to specify a non-spatial task and a non-spatial resampling\index{resampling} strategy.
<!-- @Patrick: talk in person; but I think this step is not necessary since spatial and non-spatial partitions must be different -->
<!-- Before running the spatial cross-validation for both tasks set a seed to make sure that both use the same partitions which in turn guarantees comparability.-->
1. Model landslide susceptibility using a quadratic discriminant analysis [QDA, @james_introduction_2013].
Assess the predictive performance (AUROC) of the QDA. 
What is the difference between the spatially cross-validated mean AUROC value of the QDA and the GLM?
<!-- so I think, setting a seed here makes sure that the same spatial partitions are used for both models, right?-->
Hint: Before running the spatial cross-validation for both learners, set a seed to make sure that both use the same spatial partitions which in turn guarantees comparability.
1. Run the SVM without tuning the hyperparameters.
Use the `rbfdot` kernel with $\sigma$ = 1 and *C* = 1. 
Leaving the hyperparameters unspecified in **kernlab**'s `ksvm()` would otherwise initialize an automatic non-spatial hyperparameter tuning.
For a discussion on the need for (spatial) tuning of hyperparameters, please refer to @schratz_performance_nodate.
