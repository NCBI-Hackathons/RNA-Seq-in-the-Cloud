FROM jupyter/minimal-notebook:7db1bd2a7511
#COPY requirements.txt /src/
#RUN pip install --upgrade -r /src/requirements.txt
RUN pip install \
    numpy \
    pandas \
    seaborn \
    sklearn \
    plotly \
    matplotlib \
    dash


COPY pca_dash_app.py /src/pca_dash_app.py
EXPOSE 8080
EXPOSE 8888

CMD ["python", "/src/pca_dash_app.py"]
