Downloading occurrence data from GBIF
================
Achyut K Banerjee
August 22, 2021

**Premise:** To download occurrence data for **a large number of
species**, the following code has been used. The example in this code is
for the ILORA database version 1.0.

**Data availability:** The associated data is available in
\[C:/Users/Achyut/Desktop/web/general\]

*Install and load package*

``` r
install.packages("rgbif")
# fill in your gbif.org credentials 
user <- "banerjeeachyut31" # your gbif.org username 
pwd <- "**********" # your gbif.org password
email <- "banerjeeachyut31@gmail.com" # your email
```

*Load libraries*

``` r
library(dplyr)
library(purrr)
library(readr)  
library(magrittr) # for %T>% pipe
library(rgbif) # for occ_download
library(taxize) # for get_gbifid_
```

*Read data*

``` r
data<-read.csv("C:/Users/Achyut/Desktop/web/sdm/sdm6/species_list.csv")
```

*Get data*

``` r
gbif_taxon_keys <- 
  readr::read_csv("species_list.csv") %>% 
  pull("Taxon name") %>% 
  taxize::get_gbifid_(method="backbone") %>% # match names to the GBIF backbone to get taxonkeys
  imap(~ .x %>% mutate(original_sciname = .y)) %>% # add original name back into data.frame
  bind_rows() %T>% # combine all data.frames into one
  readr::write_tsv(path = "all_matches.tsv") %>% # save as side effect for you to inspect if you want
  filter(matchtype == "EXACT" & status == "ACCEPTED") %>% # get only accepted and matched names
  filter(kingdom == "Plantae") %>% # remove anything that might have matched to a non-plant
  pull(usagekey) # get the gbif taxonkeys
```

*Download data*

``` r
occ_download(
  pred_in("taxonKey", gbif_taxon_keys),
  format = "DWCA",
  pred("country", "IN"),
  pred("hasCoordinate", TRUE),
  pred("hasGeospatialIssue", FALSE),
  user=user,pwd=pwd,email=email
)
```

**END**

**References**

-   \[Read the blog here\]
    (<https://data-blog.gbif.org/post/downloading-long-species-lists-on-gbif/>)
