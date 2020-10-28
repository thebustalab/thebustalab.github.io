non versus parametric
https://statisticsbyjim.com/hypothesis-testing/nonparametric-parametric-tests/

source("https://thebustalab.github.io/R_For_Chemists/custom_functions/chem.R")


hawaii_aquifers

K_data <- hawaii_aquifers %>%
	dplyr::filter(analyte == "K")

## Single mean

	K_data_1_2 <- K_data %>%
		dplyr::filter(aquifer_code %in% c("aquifer_1", "aquifer_6"))

	K_data_1_2

	ggplot(K_data_1_2, aes(x = aquifer_code, y = abundance)) +
		geom_boxplot() +
		geom_point()

	## Is the data normally distributed with similar variance?

		## Visually?

			K_data_1_2 %>%
			 	group_by(aquifer_code) %>%
				ggplot(aes(x = abundance)) + 
					geom_histogram(bins = 30) +
					facet_wrap(~aquifer_code) +
					geom_density(aes(y = ..density..*10), colour = "blue")

		## Statistically?

			# Look within each group by grouping and shapiro test
				K_data_1_2 %>%
					group_by(aquifer_code) %>% 
					shapiro_test(abundance)

			# Look across groups by not grouping and running levene
				K_data_1_2 %>%
					levene_test(abundance ~ aquifer_code)


	## If assumptions are met, use the more powerful t-test
		K_data_1_2 %>%
			t_test(abundance ~ aquifer_code)

	## If assumptions are not met, use the less powerful wilcox test
		K_data_1_2 %>%
			wilcox_test(abundance ~ aquifer_code)

## Multiple means
	
	## What is the number of measurements in each group?

		K_data %>%
			group_by(aquifer_code) %>%
			summarize(count = n())

	## Is each group normally distributed with similar variance?

		## Visually?

			K_data %>%
			 	group_by(aquifer_code) %>%
				ggplot(aes(x = abundance)) + 
					geom_histogram(bins = 30) +
					facet_wrap(~aquifer_code) +
					geom_density(aes(y = ..density..*10), colour = "blue")

		## Statistically?

			# Here we want to look at each aquifer individually - is there a normal distribution?
				K_data %>%
					group_by(aquifer_code) %>% 
					shapiro_test(abundance)

			# Here we want to look across aquifers - are variances equal?
				K_data %>%
					levene_test(abundance ~ aquifer_code)


			# If yes, ANOVA:
				K_data %>%
					anova_test(abundance ~ aquifer_code)

				# If significant, tukey
					K_data %>%
						tukey_hsd(abundance ~ aquifer_code) %>%
						p_groups()

			# If no, less powerful kruskal test:
				K_data %>%
					kruskal_test(abundance ~ aquifer_code)

				# If significant, dunn
					K_data %>%
						dunn_test(abundance ~ aquifer_code) %>%
						p_groups()