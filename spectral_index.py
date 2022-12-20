from func import *
from sed import *
from datos_completos import *
import scipy.stats as stats
import pingouin as pg

# Determinar índices espectrales

# tab es un DataFrame con las observaciones de una fuente
# fija
# filtro1 y filtro2 son cadenas que indican el indice
# de color a usar: filtro1 - filtro2
# filtro_x indica la magnitud a usar como
# variable independiente
def tabla_color_mag(tab,filtro1,filtro2,filtro_x):
    tabla_nueva = pd.DataFrame(columns = [
        'Noche',
        'x',
        'x_err',
        'Color_index',
        'Color_index_err'
    ])
    # Obtenemos las tablas con las observaciones
    # cuasi-simultaneas
    tablas_simultaneas = detecta_simultaneas(tab)
    for tabla in tablas_simultaneas:
        if len(tabla['Filter'].unique()) >= 2:
            epoca = int(tabla.MJD.values[0])
            proof1 = extraer_por_filtro(filtro1, tabla)
            proof2 = extraer_por_filtro(filtro2, tabla)
            proof3 = extraer_por_filtro(filtro_x, tabla)
            if (len(proof1['MJD'].values) == 0)\
                 or (len(proof2['MJD'].values) == 0)\
                 or (len(proof3['MJD'].values) == 0):
                continue
            # Tomar el pormedio de la magnitud en caso
            # de que haya mas de una por noche
            if len(proof1['Mag'].values) > 1:
                mags1 = proof1['Mag'].values
                mag1 = sum(mags1)/len(mags1)
                mag1_err = max(
                            [abs(mag1-mag) for mag in mags1]+\
                            [err for err in proof1['Mag_error'].values]
                        )
            else:
                mag1 = proof1['Mag'].values[0]
                mag1_err = proof1['Mag_error'].values[0]
            if len(proof2['Mag'].values) > 1:
                mags2 = proof2['Mag'].values
                mag2 = sum(mags2)/len(mags2)
                mag2_err = max(
                            [abs(mag2-mag) for mag in mags2]+
                            [err for err in proof2['Mag_error'].values]
                        )
            else:
                mag2 = proof2['Mag'].values[0]
                mag2_err = proof2['Mag_error'].values[0]
            if len(proof3['Mag'].values) > 1:
                Xs = proof3['Mag'].values
                x = sum(Xs)/len(Xs)
                x_err = max(
                            [abs(x-mag) for mag in Xs]+
                            proof3['Mag_error'].values
                        )
            else:
                x = proof3['Mag'].values[0]
                x_err = proof3['Mag_error'].values[0]
            color = mag1-mag2
            color_err = (mag1_err**2 + mag2_err**2)**0.5
            tabla_nueva = pd.concat([tabla_nueva,pd.DataFrame({
                'Noche' : epoca,
                'x' : x,
                'x_err' : x_err,
                'Color_index' : color,
                'Color_index_err' : color_err
            },index=[0])],ignore_index=True)
    return tabla_nueva

# Funcion generica para obtener un ajuste lineal, con sus
# incertidumbres de un conjunto de datos con incertidumbres
def min_cuad_con_incertidumbres(x,y,dx):
    n = len(x)
    x_arr = np.array(x)
    y_arr = np.array(y)
    dx_arr = np.array(dx)
    sumxy = (x_arr*y_arr).sum()
    sumx2 = (x_arr**2).sum()
    sumx = x_arr.sum()
    sumy = y_arr.sum()
    epsilon = max(dx_arr)
    
    a = n*sumxy-sumx*sumy
    d = n*sumx2-sumx**2
    f = sumx2*sumy-sumxy*sumx
    m = a/d
    b = f/d
    if n <= 2:
        deltam = 0
        deltab = 0
    else:
        deltam = epsilon*(n/d)**0.5
        deltab = epsilon*(1/n)**0.5
    return m,b,deltam,deltab

def ajuste_lineal_color_mag(tabla):
    mag_x = tabla['x'].values
    color = tabla['Color_index'].values
    mag_x_err = tabla['x_err'].values
    ajuste = min_cuad_con_incertidumbres(
        mag_x, 
        color,
        mag_x_err
    )
    return ajuste

def ajuste_lineal_esp_mag(tabla,frecuencias):
    mag_x = tabla['x'].values
    mag_x_err = tabla['x_err'].values
    color = tabla['Color_index'].values
    color_err = tabla['Color_index_err'].values
    indices, ind_err = indice_espectral(color, color_err, frecuencias)
    ajuste = min_cuad_con_incertidumbres(
        mag_x,
        indices,
        mag_x_err
    )
    return ajuste

names = [
    'oj287',
    's0954',
    'on231',
    '3c273',
    '3c279',
    'p1510',
    'da406',
    'mk501',
    'bllac'
    ]
# Funcion para graficar el diagrama color-magnitud
# para un indice de color dado de todas las fuentes
# en la muestra
def diagrama_color_mag(filtro1,filtro2,filtro_x,tabla=new_datos):
    fig, axs = plt.subplots(3,3)
    posiciones = [
        (0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)
        ]
    size = 18
    for i in range(len(names)):
        tab = tabla_color_mag(
            extraer_por_nombre(names[i], tabla), 
            filtro1, 
            filtro2,
            filtro_x
            )
        axs[posiciones[i]].errorbar(
            tab['x'],
            tab['Color_index'],
            xerr = tab['x_err'],
            yerr = tab['Color_index_err'],
            fmt = 'dk',
            ecolor = '#a3a3c2',
            barsabove = False,
            markersize = 10,
            capsize = 5,
        )
        axs[posiciones[i]].set_title(
                nombres_largos[names[i]],
                loc = 'center',
                y = 0.9,
                fontsize = size
                )
        x = np.linspace(min(tab['x']),max(tab['x']),10)
        m,b,dm,db = ajuste_lineal_color_mag(tab)
        axs[posiciones[i]].plot(
            x,
            m*x+b,
            color = 'r'
        )
    for ax in axs.flat:
        ax.set(
            xlabel = filtro_x+' (mag)',
            ylabel = filtro1+'-'+filtro2+' (mag)'
            )
        ax.xaxis.label.set_size(size)
        ax.yaxis.label.set_size(size)
        ax.tick_params(labelsize=size)
        xlims = ax.get_xlim()
        rango = xlims[1]-xlims[0]
        ax.xaxis.set_ticks(
            [round(rango*i/3+xlims[0],3) for i in range(0,4)]
            )
    fig.set_figheight(13)
    fig.set_figwidth(15)
    fig.tight_layout()
    plt.subplots_adjust(wspace=0.3)
    name = filtro1+'-'+filtro2
    plt.savefig(
        "D:/Documentos/Tesis/AnalisisDatos/DiagramasCM/"
        +name
        +".png"
        )
    plt.show()

# Funcion para realizar el analsis de correlacion
# y obtener parametros asociados a la grafica 
# color vs mag
# Recibe los filtros a usar para el diagrama
# y parametros para determinar la significancia
# de la correlacion
def correlacion_ind_color(f1,f2,f_x,p_cut=0.01,r_cut=0.2):
    redondeo = 2
    tabla = new_datos
    for i in range(len(names)):
        eti = 'pearson'
        nombre = nombres_largos[names[i]]
        proof = extraer_por_nombre(names[i], tabla)
        deltaT = int(max(proof['MJD']))-int(min(proof['MJD']))
        tab = tabla_color_mag(proof, f1, f2, f_x)
        if len(tab['Color_index'].values) <= 1:
            continue
        else:
            S,b,dS,db = ajuste_lineal_color_mag(tab)
            S,dS = round(S,redondeo),round(dS,redondeo)
            df = pg.corr(
                tab['x'],
                tab['Color_index'],
                method=eti
                )
            r = round(df['r'][eti],redondeo)
            obs = df['n'][eti]
            dr = 0.6745*(1-r**2)/(obs**0.5)
            dr = round(dr,redondeo)
            pvalue = df['p-val'][eti]
            ind_esp,dind = indice_espectral_promedio(
                tab['Color_index'], 
                tab['Color_index_err'], 
                [frecuencias[f1],frecuencias[f2]]
                )
            ind_esp = round(ind_esp,redondeo)
            dind = round(dind,redondeo)
            if (r > r_cut) & (pvalue < p_cut):
                comportamiento = 'BWB'
            elif (r < -r_cut) & (pvalue < p_cut):
                comportamiento = 'RWB'
            else:
                comportamiento = 'None'
            if pvalue < 1/100:
                pvalue = format(pvalue, '.1E')
            else:
                pvalue = round(pvalue,redondeo)
            print(
                nombre,'&',
                str(deltaT),'&',
                str(S),'('+str(dS)+')','&',
                str(r),'('+str(dr)+')','&',
                str(pvalue),'&',
                str(obs),'&',
                str(ind_esp),'('+str(dind)+')','&',
                comportamiento
            )

# colores es una lista con los valores de indice 
# de color para obtener el indice espectral
# frecuencias es una lista con las frecuencias de 
# los filtros usados para obtener el color
# el orden de las frecuencias sigue el orden 
# de los filtros para obtener el índice de color
def indice_espectral(colores,colores_err,frecuencias):
    indices = []
    indices_err = []
    for i in range(len(colores)):
        indice = 0.4 * colores[i]/np.log(frecuencias[0]/frecuencias[1])
        error = 0.4 * colores_err[i]/np.log(frecuencias[0]/frecuencias[1])
        indices.append(indice)
        indices_err.append(error)
    return indices, indices_err

def indice_espectral_promedio(colores,colores_err,frecuencias):
    color_ind_prom = np.mean(colores)
    indice_promedio = 0.4*color_ind_prom/\
            np.log(frecuencias[0]/frecuencias[1])
    err_ind = max(
        [abs(x-color_ind_prom) for x in colores]
        +list(colores_err)
        )
    err_ind = 0.4*err_ind/np.log(frecuencias[0]/frecuencias[1])
    return indice_promedio,err_ind

# Funcion para graficar los indices espectrales 
# vs MJD o mag para ver su variacion
def graficar_ind_esp(f1,f2,fx,tabla=new_datos,MJD=True):
    fig, axs = plt.subplots(3,3)
    posiciones = [
        (0,0),(0,1),(0,2),(1,0),(1,1),(1,2),(2,0),(2,1),(2,2)
        ]
    size = 18
    for i in range(len(names)):
        tab = tabla_color_mag(
            extraer_por_nombre(names[i], tabla),
            f1,
            f2,
            fx
            )
        indices, err = indice_espectral(
            tab['Color_index'].values, 
            tab['Color_index_err'].values, 
            [frecuencias[f1],frecuencias[f2]]
        )
        if MJD:
            x = np.array(tab['Noche'].values)-59300
            xerr = np.zeros(len(tab['Noche'].values))
            xlab = 'MJD (+59300)'
            eti = 'MJD'
        else:
            x = np.array(tab['x'].values)
            xerr = np.array(tab['x_err'].values)
            xlab = fx+' (mag)'
            eti = 'mag'
        axs[posiciones[i]].errorbar(
            x,
            indices,
            xerr = xerr,
            yerr = err,
            fmt = 'dk',
            ecolor = '#a3a3c2',
            barsabove = False,
            markersize = 10,
            capsize = 5
        )
        axs[posiciones[i]].set_title(
            nombres_largos[names[i]],
            loc = 'center',
            y = 0.9,
            fontsize = size
        )
        if not MJD:
            x = np.linspace(min(tab['x']),max(tab['x']),10)
            m,b,dm,db = ajuste_lineal_esp_mag(
                tab,
                [frecuencias[f1],frecuencias[f2]]
            )
            axs[posiciones[i]].plot(
                x,
                m*x+b,
                color = 'r'
            )
            if len(indices)<=2:
                ax = axs[posiciones[i]]
                xl = ax.get_xlim()
                rang = xl[1]-xl[0]
                ticks = [
                    round(rang*i/2+xl[0],3) for i in range(0,3)
                ]
                ax.xaxis.set_ticks(ticks)
    for ax in axs.flat:
        ax.set(
            xlabel = xlab,
            ylabel = r"$\alpha_{"+f1+f2+"}$"
        )
        ax.xaxis.label.set_size(size)
        ax.yaxis.label.set_size(size)
        ax.tick_params(labelsize = size)
    fig.set_figheight(13)
    fig.set_figwidth(15)
    plt.tight_layout()
    plt.subplots_adjust(wspace=0.3)
    name = tabla["Short_name"].values[0]
    plt.savefig(
        "D:/Documentos/Tesis/AnalisisDatos/SpecVar/"
        + f1
        + f2
        +'vs'
        + eti
        + ".png"
    )
    plt.show()

# El funcionamiento de esta funcion es equivalente a la
# correspondiente para el diagrama color-magnitud
def correlacion_ind_esp(f1,f2,fx,p_cut=0.01,r_cut=0.2):
    redondeo = 2
    tabla = new_datos
    for i in range(len(names)):
        eti = 'pearson'
        nombre = nombres_largos[names[i]]
        proof = extraer_por_nombre(names[i], tabla)
        deltaT = int(max(proof['MJD'])-min(proof['MJD']))
        tab = tabla_color_mag(proof, f1, f2, fx)
        if len(tab['Color_index'].values) <= 1:
            continue
        else:
            colores = tab['Color_index'].values
            colores_err = tab['Color_index_err'].values
            freqs = [frecuencias[f1],frecuencias[f2]]
            S,_,dS,_ = ajuste_lineal_esp_mag(
                tab, 
                freqs
            )
            S,dS = round(S,redondeo),round(dS,redondeo)
            indices,_ = indice_espectral(colores,colores_err,freqs)
            df = pg.corr(
                tab['x'],
                indices,
                method = eti
            )
            r = round(df['r'][eti],redondeo)
            obs = df['n'][eti]
            dr = 0.6745*(1-r**2)/(obs**0.5)
            dr = round(dr,redondeo)
            pvalue = df['p-val'][eti]            
            if (r > r_cut) & (pvalue < p_cut):
                comportamiento = 'FWB'
            elif (r < -r_cut) & (pvalue < p_cut):
                comportamiento = 'SWB'
            else:
                comportamiento = 'None'
            if pvalue < 1/100:
                pvalue = format(pvalue, '.1E')
            else:
                pvalue = round(pvalue,redondeo)
            print(
                nombre,'&',
                str(deltaT),'&',
                str(S),'('+str(dS)+')','&',
                str(r),'('+str(dr)+')','&',
                str(pvalue),'&',
                str(obs),'&',
                comportamiento
            )