# Deployment 

The Dockerfile in this repository can be used to create a production-ready image of the application.
The image can be built using the following command:

```bash
cd deployment
docker build -t <image-name> .
```
or a prebuilt one can be pulled from ghcr.io:

```bash
docker pull ghcr.io/rinikerlab/lightweight-registration:main
```


