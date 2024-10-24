#!/bin/bash
############### By: WMAlex12#################### 
# Este archivo ejecutable te permite realizar las descargas desde la base
# de datos de MGnify, utilizando como indicador de selección el estudio
# (studies).Se realizó un Script de R el cual permite hacer la descarga
# Utilizando la API, por lo que puedes acceder a ella y hacer todas las
# modificaciones que creas necesarias par ala descarga de tus datos. 

echo "Introduce el número de acceso:"
read study

# Crear o sobrescribir el archivo con el número de acceso
echo -n "$study" > study.txt

echo "Ejecuntando el Script"
Rscript MGnify_DownloadDB.R
