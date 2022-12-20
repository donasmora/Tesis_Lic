import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

### Para separar las tablas y graficar datos

def extraer_por_nombre(name,tabla):
    indices = []
    lista = tabla['Short_name'].values
    for i in range(len(lista)):
        if lista[i] == name:
            indices.append(i)
    df = tabla.iloc[indices]
    return df

def extraer_por_filtro(filtro,tabla):
    indices = []
    lista = tabla['Filter'].values
    for i in range(len(lista)):
        if lista[i] == filtro:
            indices.append(i)
    df = tabla.iloc[indices]
    return df

def extraer_por_temporada2(d_min,d_max,tabla):
    indices = []
    lista = tabla['MJD'].values
    for i in range(len(lista)):
        if d_min <= lista[i] and lista[i] <= d_max:
            indices.append(i)
    df = tabla.iloc[indices]
    return df

def graficar_lightcurve_completa(name,filtros,tabla):
    fig, ax1 = plt.subplots(figsize=(10,4))
    for i in range(len(filtros)):
        if filtros[i] == 'B':
            eti = 'B'
            color = 'b'
        elif filtros[i] == 'V':
            color = 'g'
            eti = 'V'
        elif filtros[i] == 'R':
            color = 'r'
            eti = 'R'
        elif filtros[i] == 'I':
            color = 'orange'
            eti = 'I'
        sub_tabla = extraer_por_filtro(filtros[i],extraer_por_nombre(name,tabla))
        dias_modificados = np.array(sub_tabla['MJD'].values)-59300
        ax1.errorbar(dias_modificados,sub_tabla['Mag'],yerr=sub_tabla['Mag_error'],fmt='d',color=color,markersize=8,capsize=4)
    size = 15
    ejes = plt.axis()
    plt.axis(ymin=ejes[3],ymax=ejes[2])
    ax1.set_xlabel("MJD (+59300)",fontsize=size)
    ax1.set_ylabel("Mag",fontsize=size)
    plt.tick_params(axis='both', labelsize=size)
    
    if ejes[0] < 32:
        if ejes[1] < 65:
            ticks = [0,0.7]
            ticks_labels = ['Abril','Mayo']
        else:
            ticks = [0,0.2,0.8]
            ticks_labels = ['Abril','Mayo','Junio']
    elif ejes[0] < 65 and ejes[0] >= 32:
        ticks = [0,0.65]
        ticks_labels = ['Mayo','Junio']
    else:
        ticks = [0]
        ticks_labels = ['Junio']

    ax2 = ax1.twiny()
    plt.xticks(ticks,ticks_labels)
    plt.tick_params(axis='both', labelsize=size)

    plt.tight_layout()
    plt.savefig("D:/Documentos/Tesis/AnalisisDatos/LC-mag/"+name+'completa_graf.png')
    plt.show()

def graficar(name,filtros,tabla):
    plt.figure(figsize=(10,6))
    for i in range(len(filtros)):
        if filtros[i] == 'B':
            eti = 'B'
            color = 'b'
        elif filtros[i] == 'V':
            color = 'g'
            eti = 'V'
        elif filtros[i] == 'R':
            color = 'r'
            eti = 'R'
        elif filtros[i] == 'I':
            color = 'orange'
            eti = 'I'
        sub_tabla = extraer_por_filtro(filtros[i],extraer_por_nombre(name,tabla))
        dias_modificados = np.array(sub_tabla['MJD'].values)-59300
        plt.errorbar(dias_modificados,sub_tabla['Mag'],yerr=sub_tabla['Mag_error'],fmt='d',color=color,markersize=15,capsize=7)
    size = 30
    ejes = plt.axis()
    plt.axis(ymin=ejes[3],ymax=ejes[2])
    plt.xlabel("MJD (+59300)",fontsize=size)
    plt.ylabel("Mag",fontsize=size)
    plt.tick_params(axis='both', labelsize=size)
    plt.ticklabel_format(useOffset=False, style='plain')
    dist = ejes[1]-ejes[0]
    plt.xticks(np.linspace(ejes[0]+0.1*dist,ejes[1]-0.1*dist,4))
    plt.tight_layout()
    meses = tabla['Month'].values
    plt.savefig("D:/Documentos/Tesis/AnalisisDatos/LC-mag/"+name+meses[0]+'graf.png')
    plt.show()

### Para clalcular los flujos y sus incertidumbres

def mag_to_flux(mag,mag_err,filtro):
    if filtro == 'B':
        F0 = 4130
    elif filtro == 'V':
        F0 = 3781
    elif filtro == 'R':
        F0 = 2941
    elif filtro == 'I':
        F0 = 2635
    else:
        print('Error: Filter unknown')
    F = F0 * 1000*10**(-0.4*mag)
    F_err = 0.4*1000*np.log(10)*F0*10**(-0.4*mag)*mag_err
    return [round(F,2),round(F_err,2)]

### Graficar flujos

def graficar_flujo(name,filtros,tabla):
    plt.figure(figsize=(10,6))
    for i in range(len(filtros)):
        if filtros[i] == 'B':
            eti = 'B'
            color = 'b'
        elif filtros[i] == 'V':
            color = 'g'
            eti = 'V'
        elif filtros[i] == 'R':
            color = 'r'
            eti = 'R'
        elif filtros[i] == 'I':
            color = 'orange'
            eti = 'I'
        sub_tabla = extraer_por_filtro(filtros[i],extraer_por_nombre(name,tabla))
        dias_modificados = np.array(sub_tabla['MJD'].values)-59300
        plt.errorbar(dias_modificados,sub_tabla['Flux'],yerr=sub_tabla['Flux_err'],fmt='d',color=color,markersize=15,capsize=7)
    size = 30
    ejes = plt.axis()
    plt.xlabel("MJD (+59300)",fontsize=size)
    plt.ylabel(r'$F_\nu$'+" (mJy)",fontsize=size)
    plt.tick_params(axis='both', labelsize=size)
    dist = ejes[1]-ejes[0]
    plt.xticks(np.linspace(ejes[0]+0.1*dist,ejes[1]-0.1*dist,4))
    plt.ticklabel_format(useOffset=False, style='plain')
    plt.tight_layout()
    meses = tabla['Month'].values
    plt.savefig("D:/Documentos/Tesis/AnalisisDatos/LC-flux/"+name+meses[0]+'graf_flux.png')
    plt.show()

def graficar_flujo_completa(name,filtros,tabla):
    fig, ax1 = plt.subplots(figsize=(10,4))
    for i in range(len(filtros)):
        if filtros[i] == 'B':
            eti = 'B'
            color = 'b'
        elif filtros[i] == 'V':
            color = 'g'
            eti = 'V'
        elif filtros[i] == 'R':
            color = 'r'
            eti = 'R'
        elif filtros[i] == 'I':
            color = 'orange'
            eti = 'I'
        sub_tabla = extraer_por_filtro(filtros[i],extraer_por_nombre(name,tabla))
        dias_modificados = np.array(sub_tabla['MJD'].values)-59300
        ax1.errorbar(dias_modificados,sub_tabla['Flux'],yerr=sub_tabla['Flux_err'],fmt='d',color=color,markersize=8,capsize=4)
    size = 15
    ejes = plt.axis()
    ax1.set_xlabel("MJD (+59300)",fontsize=size)
    ax1.set_ylabel(r'$F_\nu$'+" (mJy)",fontsize=size)
    plt.tick_params(axis='both', labelsize=size)

    if ejes[0] < 32:
        if ejes[1] < 65:
            ticks = [0,0.7]
            ticks_labels = ['Abril','Mayo']
        else:
            ticks = [0,0.2,0.8]
            ticks_labels = ['Abril','Mayo','Junio']
    elif ejes[0] < 65 and ejes[0] >= 32:
        ticks = [0,0.65]
        ticks_labels = ['Mayo','Junio']
    else:
        ticks = [0]
        ticks_labels = ['Junio']
    
    ax2 = ax1.twiny()
    plt.xticks(ticks,ticks_labels)
    plt.tick_params(axis='both', labelsize=size)

    plt.tight_layout()
    plt.savefig("D:/Documentos/Tesis/AnalisisDatos/LC-flux/"+name+'completa_graf_flux.png')
    plt.show()

### Para identificar los outliers

def min_cuad(x,y):
    n = len(x)
    x_a = np.array(x)
    y_a = np.array(y)
    sumxy = sum(x_a*y_a)
    sumx = sum(x_a)
    sumy = sum(y_a)
    sumx2 = sum(x_a*x_a)
    a = n*sumxy-sumx*sumy
    b = n*sumx2-sumx**2
    c = sumx2*sumy-sumxy*sumx
    pend = a/b
    ordenada = c/b
    return [pend, ordenada]

def desv_stand_ajuste(x,y):
    n = len(x)
    x_a = np.array(x)
    y_a = np.array(y)
    estad = min_cuad(x_a,y_a)
    suma = sum((y_a-(estad[0]*x_a+estad[1]))**2)
    return (suma/n)**(1/2)

def outliers(x,y):
    estad = min_cuad(x,y)
    sigma = desv_stand_ajuste(x,y)
    indices = []
    for i in y.index:
        dist = abs(y[i]-estad[0]*x[i]-estad[1])
        if dist > 3*sigma:
            indices.append(i)
    return indices


### Ctes necesarias
nombres_largos = {
    'e1553' : '1ES 1553+113',
    'mk501' : 'Markarian 501',
    'on231' : 'ON 231',
    '3c273' : '3C 273',
    '3c279' : '3C 279',
    'oj287' : 'OJ 287',
    's0954' : 'S4 0954+650',
    'mk421' : 'Markarian 421',
    'oj248' : 'OJ 248',
    'p1222' : 'PKS 1222+216',
    'ot546' : 'OT 546',
    'p1510' : 'PKS 1510-089',
    'q4c38' : '4C 38.41',
    '3c371' : '3C 371',
    'da406' : 'DA 406',
    '3c345' : '3C 345',
    'bllac' : 'BL Lacertae',
    'e1959' : '1ES 1959+650',
    'q4c51' : '4C 51.37',
    '3c454' : '3C 454.3',
    'ct102' : 'CTA 102',
}

meses_comp = {
    'ABR' : 'Abril',
    'MAY' : 'Mayo',
    'JUN' : 'Junio'
}