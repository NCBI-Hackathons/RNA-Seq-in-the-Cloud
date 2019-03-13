FROM python:3-alpine
COPY requirements.txt /src/
RUN pip install --upgrade -r /src/requirements.txt
    
COPY app.py /src/app.py
EXPOSE 8080
CMD ["python", "/src/app.py"]
