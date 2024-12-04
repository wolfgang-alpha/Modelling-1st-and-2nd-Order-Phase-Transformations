#!/bin/bash 

docker run --rm --init -p 8888:8888 -v "$(pwd)":/root/shared dolfinx/lab:v0.7.2 ./shared/ch-limit.ipynb
