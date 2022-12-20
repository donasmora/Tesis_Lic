# Este archivo funciona como un preprocesamiento de los datos analizados en la tesis.

import pandas as pd
from func import *

# importas datos
columnas = ["Short_name","Date","JD","Filter","Mag","Mag_error"]
data = pd.read_table("Multicolor_Blazar.dat",sep="\s+",usecols=columnas)
# Los datos utilizados en la tesis no son públicos, sin embargo,
# se prentende hacerlos públicos como parte de un artículo a futuro

datos = data.copy() # Para no modificar la importacion original
# se arreglan renglones con nombres de fuente con errores
datos = datos.replace('2c273','3c273')
datos = datos.replace('3c270','3c279')
datos = datos.replace('3c245','3c345')
# se quitan fuentes con porblemas y renglones que estorban
nombres_a_quitar = ['oj049','FILTER','%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%','3c345']
# se quitan fuentes que tienen datos de menos de 4 filtros
names = datos['Short_name'].unique()
for name in names:
    proofs = extraer_por_nombre(name,datos)
    filtros = proofs['Filter'].unique()
    if len(filtros)<4:
        nombres_a_quitar.append(name)
indices = []
for name in nombres_a_quitar:
    for i in range(len(datos['Short_name'])):
        if datos['Short_name'][i] == name:
            indices.append(i)
new_datos = datos.drop(indices)
# se definen nuevas columnas para el análisis
new_datos['JD'] = new_datos['JD'].map(lambda x: float(x))
new_datos['MJD'] = new_datos['JD'].map(lambda x: x-2400000.5)
new_datos['Flux'] = [mag_to_flux(new_datos['Mag'][i],new_datos['Mag_error'][i],new_datos['Filter'][i])[0] for i in new_datos.index]
new_datos['Flux_err'] = [mag_to_flux(new_datos['Mag'][i],new_datos['Mag_error'][i],new_datos['Filter'][i])[1] for i in new_datos.index]
# se quitan los outliers usando los algoritmos creados para dicho fin
names = new_datos['Short_name'].unique()
cotas = [0,59332,59360,59400]
filtros = ['I','R','V','B']
for name in names:
    proofs = extraer_por_nombre(name,new_datos)
    for i in range(len(cotas)-1):
        for filtro in filtros:
            proofs2 = extraer_por_temporada2(cotas[i],cotas[i+1],extraer_por_filtro(filtro,proofs))
            if len(proofs2['MJD']) > 1:
                new_datos = new_datos.drop(outliers(proofs2['MJD'],proofs2['Mag']))
            else:
                continue
# indices cuyos outliers no se detectan con el algoritmo
# pero son considerados outliers pues no hay observaciones continuas que 
# contribuyan a los saltos presentados por estas observaciones
outliers_extra = [616,373,628,516,154,551,21,278,510,511,434]
new_datos = new_datos.drop(outliers_extra)