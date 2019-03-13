FROM python:3-alpine
COPY requirements.txt /src/
RUN pip install --upgrade -r /src/requirements.txt

COPY pca_dash_app.py /src/pca_dash_app.py
EXPOSE 8080
CMD ["python", "/src/pca_dash_app.py"]
