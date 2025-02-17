# Deployment 

The Dockerfile in this repository can be used to create a production-ready image of the application.
The image can be built using the following command from the main lwreg directory:

```bash
cd deployment
docker build -f deployment/Dockerfile -t <image-name> .
```
or a prebuilt one can be pulled from ghcr.io:

```bash
docker pull ghcr.io/rinikerlab/lightweight-registration:main
```

## Sample session using the command-line interface

Build the image and start a container:

```bash
% docker build -f deployment/Dockerfile -t lwreg_1 .
% docker run --name lwreg_demo  -d -it lwreg_1 bash

```

And now you can use the running container:
```bash
% docker exec lwreg_demo lwreg initdb --confirm=yes
% docker exec lwreg_demo lwreg register --smiles="CC"
1
% docker exec lwreg_demo lwreg register --smiles="CCO"
2
% docker exec lwreg_demo lwreg register --smiles="CCN"
3
% docker exec lwreg_demo lwreg query --smiles="CCN"
3
% docker exec lwreg_demo lwreg retrieve --id=3
(3, '\n     RDKit          2D\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.000000 0.000000 0.000000 0\nM  V30 2 C 1.299038 0.750000 0.000000 0\nM  V30 3 N 2.598076 -0.000000 0.000000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 1 2 3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n', 'mol')
% docker exec lwreg_demo lwreg query --smiles="CCOC"
not found
```
