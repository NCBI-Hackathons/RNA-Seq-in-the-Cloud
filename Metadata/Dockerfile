FROM jupyter/datascience-notebook:7f1482f5a136

#maintainer "Ariella Gladstein <aglad@med.unc.edu>"
#organization "University of North Carolina at Chapel Hill"
#department "Genetics"
#date "13 March 2019"
#application "Data science applications"


RUN conda install --quiet --yes 'simplegeneric'

# Install R packages
RUN wget https://cran.r-project.org/src/contrib/rjson_0.2.20.tar.gz
RUN R CMD INSTALL rjson_0.2.20.tar.gz


# Switch back from root
USER $NB_UID