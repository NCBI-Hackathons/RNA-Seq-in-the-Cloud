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



COPY app.py /src/app.py
EXPOSE 8080
EXPOSE 8888
CMD ["python", "/src/app.py"]
