from func import *

# Recibe un dia para detectar cuasi simultaneidad
# en una lista de dias (mjd) con base en dia
def indices_observaciones_simultaneas(dia,mjd):
    indices_cuasi = []
    dia = int(dia)
    for i in mjd.index:
        if mjd[i] < dia+1 and mjd[i] >= dia:
            indices_cuasi.append(i)
    return indices_cuasi

# Regresa una lista, cada elemento de la lista es un DataFrame con
# observaciones ocurridas en la misma noche (MJD)
def detecta_simultaneas(tabla):
    tablas_simultaneas = []
    mjd = tabla['MJD']
    while len(mjd) > 0:
        indice = mjd.index[0]
        dia = mjd[indice]
        indices_cuasi = indices_observaciones_simultaneas(dia, mjd)
        new_df = tabla.loc[indices_cuasi]
        tablas_simultaneas.append(new_df)
        tabla = tabla.drop(indices_cuasi)
        mjd = tabla['MJD']
    return tablas_simultaneas

# Longitudes efectivas de los filtros
longitudes = {
    'B' : 4330E-10,
    'V' : 5750E-10,
    'R' : 6340E-10,
    'I' : 8040E-10
}
# Ancho de los filtros, en metros
anchos_banda_metros = {
    'B' : 950E-10,
    'V' : 1400E-10,
    'R' : 400E-10,
    'I' : 1660E-10
}
c = 299_792_458
filtros = ['B','V','R','I']
freq = []
anchos = []
for i in filtros:
    freq.append(c/longitudes[i])
    anchos.append(anchos_banda_metros[i] * c/longitudes[i]**2)
# Frecuencias efectivas en Hertz de cada filtro
frecuencias = dict(zip(filtros, freq))
# Anchos de banda en unidades de frecuencia
anchos_banda = dict(zip(filtros, anchos))

# Funcion para graficar las SED dada la tabla de una fuente fija
# Grafica todas las SED de las observaciones cuasi-simultaneas
# en la misma grafica
def graficar_SED(tabla):
    plt.figure(figsize=[10,6])
    tablas_simultaneas = detecta_simultaneas(tabla)
    for tabla in tablas_simultaneas:
        if len(tabla['Filter'].unique()) >= 4:
            tab = tablas_para_graficar(tabla)
            epocas = tab['Noche'].values
            eti = int(min(epocas))-59300
            plt.errorbar(
                tab['Freq_log'],
                tab['Flux_log'],
                yerr = tab['Flux_log_err'],
                label = eti,
                fmt = 'o--',
                ecolor = '#a3a3c2',
                markersize = 10,
                capsize = 5
            )
    size = 20
    plt.xlabel(r"$\log(\nu)$"+" [Hz]",fontsize=size)
    plt.ylabel(r"$\log(\nu F_\nu)$"+" [mJy Hz]",fontsize=size)
    plt.tick_params(axis='both',labelsize=size)
    plt.xticks([14.55,14.65,14.75,14.85])
    plt.ticklabel_format(useOffset=False, style='plain')
    plt.tight_layout()
    plt.legend(
        fontsize=size-5,title='MJD (+59300)',title_fontsize=size-5
        )
    name = tabla["Short_name"].values[0]
    plt.savefig(
        "D:/Documentos/Tesis/AnalisisDatos/SEDs/"+name+"_SED.png"
        )
    plt.show()

# Recibe una tabla de observaciones cuasi-simultaneas
# y crea un DataFrame con la informacion que se necesita
# para graficar la SED
def tablas_para_graficar(tabla):
    tab = pd.DataFrame(columns = [
            'Noche',
            'Freq',
            'Flux',
            'Flux_err',
            'Freq_log',
            'Flux_log',
            'Flux_log_err'
        ])
    for filtro in tabla['Filter'].unique():
        proof = extraer_por_filtro(filtro, tabla)
        flujos = proof['Flux'].values
        flujos_err = proof['Flux_err'].values
        if len(proof['MJD']) > 1:
            flux = sum(flujos)/len(flujos)
            flux_err = max(
                    [abs(flux-flujo) for flujo in flujos]+
                    flujos_err
                )
        else:
            flux = proof['Flux'][proof['Flux'].index[0]]
            flux_err = proof['Flux_err'][proof['Flux'].index[0]]
        tab = pd.concat([tab,pd.DataFrame({
            'Noche' : proof['MJD'][proof['Flux'].index[0]],
            'Freq' : frecuencias[filtro],
            'Flux' : flux*anchos_banda[filtro], # nuF_nu
            'Flux_err' : flux_err*anchos_banda[filtro],
            'Freq_log': np.log10(frecuencias[filtro]),
            'Flux_log': np.log10(flux),
            'Flux_log_err': (1/(flux*np.log(10)))*flux_err
        },index=[0])],ignore_index=True)
    tab = tab.sort_values(by=['Freq'])
    return tab