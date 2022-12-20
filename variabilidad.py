import numpy as np
import pandas as pd
from scipy.stats import chi2
from datos_completos import *
from sed import *
import math
import spectral_index as si

def chi2_mag(magnitudes,errores):
    chi2_f = 0
    M_mean = sum(magnitudes)/len(magnitudes)
    for i in range(len(magnitudes)):
        chi2_f += (magnitudes[i]-M_mean)**2 / (errores[i]**2)
    return chi2_f

# se puede cambiar la significancia al valor que se desee
def test_chi2(magnitudes,errores,significance=0.01):
    mag = np.array(magnitudes)
    err = np.array(errores)
    dof = len(mag)-1
    chi2_calc = chi2_mag(mag,err)
    if chi2_calc > chi2_critic_values[str(significance)][dof-1]:
        return True
    else:
        return False

def amplitud_var_mag(magnitudes,errores):
    mag = np.array(magnitudes)
    err = np.array(errores)
    A_max = max(mag)
    A_min = min(mag)
    error_mean = sum(err)/len(err)
    return round(((A_max-A_min)**2-2*error_mean**2)**(1/2),2)

### Escala de variacion con el otro formalismo

def escala_variabilidad(dias,flujos):
    taus = []
    dates = np.array(dias)
    flux = np.array(flujos)
    for i in range(len(dates)-1):
        aux = np.array([flux[i],flux[i+1]])
        F_2 = max(aux)
        F_1 = min(aux)
        if F_2 == F_1:
            continue
        tau = (dates[i+1]-dates[i])/(np.log(F_2/F_1))
        taus.append(tau)
    if min(taus) == 0:
        taus.remove(0)
    return min(taus)


### Escala de variaicion minima con doubling-timescale

def min_doubling_timescale(dias,flujos):
    scales = []
    dates = np.array(dias)
    flux = np.array(flujos)
    inf = float('inf')
    for i in range(len(flujos)-1):
        F_2 = max(flux[i],flux[i+1])
        F_1 = min(flux[i],flux[i+1])
        if F_2 == F_1:
            continue
        scale = (dates[i+1]-dates[i])/np.log2(F_2/F_1)
        scales.append(scale)
    if min(scales) == 0:
        scales.remove(0)
    return min(scales)

def resultados_variabilidad(name):
    print(nombres_largos[name])
    for filtro in filtros:
        print("Filtro:",filtro)
        for i in range(len(cotas)-1):
            proof = extraer_por_temporada2(
                cotas[i],
                cotas[i+1],
                extraer_por_filtro(filtro,extraer_por_nombre(name,new_datos))
            )
            if len(proof['MJD']) > 1:
                if test_chi2(proof['Mag'],proof['Mag_error']):
                    proof = proof.sort_values(by=['MJD'])
                    resultados = escala_variabilidad(proof['MJD'],proof['Flux'])
                    doubling_scale = min_doubling_timescale(proof['MJD'], proof['Flux'])
                    cadena = "& "+str(amplitud_var_mag(proof['Mag'],proof['Mag_error']))+\
                    " & "+str(round(resultados,5))+\
                    " & "+ str(round(doubling_scale,5))
                else:
                    cadena = ''
                inicio = math.floor(min(np.array(proof['MJD'])))
                final = math.ceil(max(np.array(proof['MJD'])))
                if inicio == final-1:
                    epoca = str(inicio)
                else:
                    epoca = str(inicio)+'-'+str(final)
                print("&",epoca,
                    "&",test_chi2(proof['Mag'],proof['Mag_error']),cadena
                     )
                    
            else:
                continue

def chi2_ind_esp(indices,errores):
    chi2_ind = 0
    ind_mean = sum(indices)/len(indices)
    for i in range(len(indices)):
        chi2_ind += (indices[i]-ind_mean)**2 / (errores[i]**2)
    return chi2_ind

def amplitud_var_ind(indices,errores):
    index = np.array(indices)
    err = np.array(errores)
    ind_max = max(index)
    ind_min = min(index)
    error_mean = sum(err)/len(err)
    return round(((ind_max-ind_min)**2-2*error_mean**2)**(1/2),2)

def extraer_por_temporada_ind(d_min,d_max,tabla):
    indices = []
    lista = tabla['Epoca'].values
    for i in range(len(lista)):
        if d_min <= lista[i] and lista[i] <= d_max:
            indices.append(i)
    df = tabla.iloc[indices]
    return df

def variabilidad_espectral(f1,f2):
    for j in range(len(si.names)):
        name = nombres_largos[si.names[j]]
        tab = si.tabla_color_mag(
            extraer_por_nombre(si.names[j], new_datos),
            f1,
            f2,
            f2
        )
        print(name)
        proof = tab.sort_values(by=['Epoca'])
        if len(proof['Epoca']) > 1:
            indices, err = si.indice_espectral(
                proof['Color_index'].values, 
                proof['Color_index_err'].values, 
                [frecuencias[f1],frecuencias[f2]]
            )
            if test_chi2(indices, err):
                resultado = amplitud_var_ind(indices, err)
                cadena = "& "+str(resultado)
            else:
                cadena = ''
            print(
                "&",
                math.floor(min(np.array(proof['Epoca']))),
                '-',
                math.ceil(max(np.array(proof['Epoca']))),
                "&",
                test_chi2(indices,err),
                cadena
            )
        else:
            continue

### Importar una tabla de valores críticos de chi2 para usar directamente en la prueba
chi2_critic_values = pd.read_csv("chi2.csv")