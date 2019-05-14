# Pipeline to find samples based on their metadata

The interactive pipeline is accessible through the Jupyter Notebooks 
- https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Metadata/case_control_search.ipynb
- https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Metadata/series_search.ipynb

## How to interact with the Jupyter Notebooks

### Accessing Notebooks via Docker
The Notebook can be run via a Docker container so you do not have to worry about dependencies!  
*Note: If you do not have Docker installed, see the Docker documentation 
https://docs.docker.com/install/*  

On your local machine, simply run:
```
sudo docker run -it -v $(pwd):/home/jovyan/work --rm -p 8888:8888 ncbihackathons/metadata jupyter-lab
```
You must be in the `RNA-Seq-in-the-Cloud` or `RNA-Seq-in-the-Cloud/Metadata` directory.

You will then see something like:
```
[I 23:49:14.625 LabApp] Writing notebook server cookie secret to /home/jovyan/.local/share/jupyter/runtime/notebook_cookie_secret
[I 23:49:15.165 LabApp] JupyterLab extension loaded from /opt/conda/lib/python3.6/site-packages/jupyterlab
[I 23:49:15.165 LabApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
[W 23:49:15.167 LabApp] JupyterLab server extension not enabled, manually loading...
[I 23:49:15.174 LabApp] JupyterLab extension loaded from /opt/conda/lib/python3.6/site-packages/jupyterlab
[I 23:49:15.174 LabApp] JupyterLab application directory is /opt/conda/share/jupyter/lab
[I 23:49:15.175 LabApp] Serving notebooks from local directory: /home/jovyan
[I 23:49:15.176 LabApp] The Jupyter Notebook is running at:
[I 23:49:15.176 LabApp] http://(a8e01d3931bb or 127.0.0.1):8888/?token=53c5dcc8cab1af007ff3a7cabf41201e65771f1846ec75d4
[I 23:49:15.176 LabApp] Use Control-C to stop this server and shut down all kernels (twice to skip confirmation).
[C 23:49:15.177 LabApp]

    Copy/paste this URL into your browser when you connect for the first time,
    to login with a token:
        http://(a8e01d3931bb or 127.0.0.1):8888/?token=53c5dcc8cab1af007ff3a7cabf41201e65771f1846ec75d4
```
Then, paste the URL and token into your browser (we recommend using `127.0.0.1`), e.g.:  
`http://127.0.0.1:8888/?token=53c5dcc8cab1af007ff3a7cabf41201e65771f1846ec75d4`.  
*Note: If you run this Docker image via ssh on an HPC, server, or virtual machine, paste the IP of the machine it is running on, not `127.0.0.1`.*

In your browser, Jupyter Lab will open and you will see something like:  
![Alt text](https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Metadata/jupyterlab-home.png?raw=true "Title")

You can then navigate to the file you want to interact with:  
![Alt text](https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Metadata/jupyterlab-notebook.png?raw=true "Title")

### Accessing Notebooks via Singularity
If you do not have root access on your machine, you can also run the same Docker image via Singularity (https://www.sylabs.io/docs/)!
```
singularity exec docker://ncbihackathons/metadata jupyter-lab
```

### Accessing Notebook without Docker or Singulaity
__Not Recommended__

It is not necessary to access the notebook via the Docker image.
If you opt to not use the Docker image, you will need to install:
- Jupyter Notebooks (or Jupyter Lab)
- All the dependencies for the Notebook

## Additional information on the functions
- The Python functions can be found in https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Metadata/functions.py
- The R code can be found in 
  - https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Metadata/Metadata_plot.R
  - https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Metadata/Metadata_table.R
  - https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Metadata/Metadata_piecharts.R
