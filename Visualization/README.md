# Visualization of RNAseq data

[Mybinder](https://hub.gke.mybinder.org/user/ncbi-hackathons-eq-in-the-cloud-i1a58q64/notebooks/Visualization/Visualizations.ipynb) version of the github.





PCA
![PCA](https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Visualization/3D_pca.PNG)

Heatmap
![Heatmap](https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Visualization/interactive_heatmap.PNG)

Volcano Plot
![Vol plt](https://github.com/NCBI-Hackathons/RNA-Seq-in-the-Cloud/blob/master/Visualization/basic_heatmap.PNG)


### Running with Docker

Build the Docker image (This step can be skipped when the Docker image has been pushed to Docker Hub)  
*must be in `Visualization` directory*
```bash
sudo docker build -t ncbihackathons/viz-notebook .
```

Run Jupyter Lab via Docker container:
```bash
sudo docker run -it -v $(pwd):/home/jovyan/work --rm -p 8888:8888 ncbihackathons/viz-notebook jupyter-lab
```

If you want to enter bash in the container:
```
sudo docker run -it -v $(pwd):/home/jovyan/work --rm -p 8888:8888 ncbihackathons/viz-notebook /bin/bash
```
To run PCA dash app:
```bash
python PCA_dash_app.py -cd /path/to/results/table.txt -a /path/to/attributes/table.txt
```
