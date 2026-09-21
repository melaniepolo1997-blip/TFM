#%%
#======================
#Cargar librerías 
#======================

import pandas as pd 
import requests
from io import StringIO
import re
import numpy as np
import time
import unicodedata
import os
from pathlib import Path

#%%

Base_DIR = Path(__file__).resolve().parent

Data_DIR = Base_DIR / "data"
Results_DIR = Base_DIR / "results"

Results_DIR.mkdir(exist_ok=True)
#%%
#===============================================
#Conectar a KEGG PATHWAY Database
#===============================================

url_kegg_pathway = "https://rest.kegg.jp/list/pathway"

respuesta = requests.get(url_kegg_pathway,
                         timeout=60)
respuesta.raise_for_status()
#Significado de resultados de respuesta:
#200 → solicitud correcta
#400 → solicitud incorrecta
#404 → recurso no encontrado
#500 → problema en el servidor

#==================================
#Descargar las rutas metabólicas
#==================================

datos = StringIO(respuesta.text)

rutas = pd.read_csv(
    datos,
    sep = "\t",
    header = None,
    names=["Pathway_ID", "Pathway"])

#==================================
#Búsqueda Hidrocarburos de interés
#==================================

#Primer filtro usando la lista de hidrocarburos del dataframe de clasificación
hidrocarburos = pd.read_csv(Data_DIR / "Clasificación tipo hidrocarburos.csv")
print(hidrocarburos["Ruta de degradacion"])
busqueda = "|".join(
    re.escape(hidrocarburo)
    for hidrocarburo in hidrocarburos["Hidrocarburo"]
)
print(busqueda)
ruta_hidrocarburos = rutas["Pathway"].str.contains(
    busqueda,
    case=False,
    na=False
)
rutas_filtro1 = rutas[ruta_hidrocarburos].copy()
print(rutas_filtro1)

#Segundo filtro: buscamos que exista la palabra degradation
filtro_degradation = rutas_filtro1["Pathway"].str.contains(
    "degradation",
    case = False,
    na = False)
rutas_filtro2 = rutas_filtro1[filtro_degradation].copy()
print(rutas_filtro2)

#Añadimos un nuevo filtro para todas las rutas
#Este filtro nos permitirá añadir rutas generales de degradacion de hidrocarburos
#Primero vemos la clasificación por categoría
print(hidrocarburos["Tipo hidrocarburo"].unique())
#Creamos una varialble con nuestros tipos de hidrocarburos
categorias_clasificacion = [
    "polycyclic aromatic hydrocarbon",
    "monocyclic aromatic hydrocarbon",
    "aliphatic hydrocarbon"]

busqueda_general = "|".join(
    re.escape(ruta)
    for ruta in categorias_clasificacion)
filtro_pathways_generales = rutas["Pathway"].str.contains(
    busqueda_general,
    case = False,
    na = False)
rutas_generales_kegg = rutas[filtro_pathways_generales].copy()
print("Número de rutas generales encontradas:")
print(rutas_generales_kegg.shape[0])

print("\nRutas encontradas:")
print(rutas_generales_kegg)

#Unimos las rutas que encontramos a partir de los hidrocarburos y las rutas generales
rutas_candidatas = pd.concat(
    [rutas_filtro2, rutas_generales_kegg]
    ).drop_duplicates(
        subset="Pathway_ID"
        ).reset_index(drop=True)
print(rutas_candidatas)
#%%
#=====================================================
#Obtención de información de Modules de cada Pathway
#=====================================================

#Definimos una funcion para poder obtener la info de todas las rutas
def obtener_info_pathways(pathway_id):
    url_2 = f"https://rest.kegg.jp/get/{pathway_id}"
    respuesta = requests.get(url_2)
    if respuesta.status_code != 200: 
        print("Error en pathway:", pathway_id)
        return None
    return respuesta.text

#Guardamos la informacion que obtenemos
informacion_pathways = {}

for pathway_id in rutas_candidatas["Pathway_ID"]:

    informacion_pathways[pathway_id] = obtener_info_pathways(pathway_id)

print("\nPathways:", len(informacion_pathways))

#=======================================
#Extraer la información de los MODULES
#=======================================
datos_modules = []

for pathway_id, informacion in informacion_pathways.items():
    lineas = informacion.split("\n")
    guardar = False
    for linea in lineas:
        if linea.startswith("MODULE"):
            guardar = True
            partes = linea.split()
            if len(partes) >= 3:
                module_id = partes[1]
                module_name = " ".join(partes[2:])
                datos_modules.append({
                    "Pathway_ID": pathway_id,
                    "Module_ID": module_id,
                    "Module": module_name
                })

        elif guardar and linea.strip() and not linea.startswith(" "):
            guardar = False
        elif guardar:
            partes = linea.split()
            if partes and re.match(r"^M\d+$", partes[0]):
                module_id = partes[0]
                module_name = " ".join(partes[1:])
                datos_modules.append({
                    "Pathway_ID": pathway_id,
                    "Module_ID": module_id,
                    "Module": module_name
                })
tabla_modules = pd.DataFrame(datos_modules)
print(tabla_modules)

#=======================
#Obtener MODULES únicos
#=======================
#Normalizar
tabla_modules_limpia = tabla_modules.copy()

tabla_modules_limpia["Module_ID"] = (
    tabla_modules_limpia["Module_ID"]
    .astype(str)
    .str.strip()
)
#Agrupamos
modules_unicos = (
    tabla_modules_limpia
    .groupby("Module_ID", as_index=False)
    .agg(
        Module=("Module", "first"),
        Pathways=("Pathway_ID", lambda x: "; ".join(sorted(set(x))))
    )
)
print(modules_unicos.head(20))

#==================================
#información de cada MODULE único
#==================================

def obtener_info_module(module_id):   
    url_modules = f"https://rest.kegg.jp/get/{module_id}" 
    respuesta = requests.get(url_modules)    
    if respuesta.status_code != 200:
        print("Error en:", module_id)
        return None    
    return respuesta.text

informacion_modules = {}

for module_id in modules_unicos["Module_ID"]:
    informacion_modules[module_id] = obtener_info_module(module_id)
#%%
#==================================
#Obtener los KO de cada MODULE 
#==================================
datos_ko = []

for module_id, informacion in informacion_modules.items():    
    lineas = informacion.split("\n")   
    guardar = False   
    for linea in lineas:       
        if linea.startswith("ORTHOLOGY"):
            guardar = True            
            contenido = linea[len("ORTHOLOGY"):].strip()        
        elif guardar and linea.strip() and not linea.startswith(" "):
            guardar = False
            continue        
        elif guardar:
            contenido = linea.strip()       
        else:
            continue        
        kos = re.findall(r"K\d{5}", contenido)        
        if not kos:
            continue        
# Eliminar KO de la línea para conservar con la descripción
        descripcion = re.sub(
            r"K\d{5}(?:\+K\d{5})*",
            "",
            contenido,
            count=1
        ).strip()
        
        for ko in kos:
            
            datos_ko.append({
                "Module_ID": module_id,
                "KO": ko,
                "Function": descripcion
            })
#Guardar información en un dataframe
tabla_ko = pd.DataFrame(datos_ko)
print("Módulos en modules_unicos:", modules_unicos["Module_ID"].nunique())
print("Módulos con KO:", tabla_ko["Module_ID"].nunique())
#Modules sin KO
modulos_sin_ko = modules_unicos[
    ~modules_unicos["Module_ID"].isin(tabla_ko["Module_ID"])
]

print(modulos_sin_ko)
#Unimos la información de MODULES y KO
tabla_modules_ko = modules_unicos.merge(
    tabla_ko,
    on="Module_ID",
    how="left"
)
print(tabla_modules_ko.head(20))
#%%
#======================================================================
#Buscar organismos relacionados a cada KO
#======================================================================
datos_organismos = []

for module_id, informacion in informacion_modules.items():    
    lineas = informacion.split("\n")   
    guardar = False    
    for linea in lineas:       

#Inicio de la sección COMPLETE
        if linea.startswith("COMPLETE"):
            guardar = True
            
            contenido = linea[len("COMPLETE"):].strip()
            
        #Fin de la sección COMPLETE
        elif guardar and linea.strip() and not linea.startswith(" "):
            guardar = False
            continue
        
        # Líneas continuadas de COMPLETE
        elif guardar:
            contenido = linea.strip()
        
        else:
            continue
        
        # Si no hay contenido, continuamos
        if not contenido:
            continue
        
        # Dividir código de organismo y nombre
        partes = contenido.split()
        
        if len(partes) >= 2:
            
            organism_code = partes[0]
            organism = " ".join(partes[1:])
            
            datos_organismos.append({
                "Module_ID": module_id,
                "Organism_code": organism_code,
                "Organism": organism
            })
tabla_organismos = pd.DataFrame(datos_organismos)

print(tabla_organismos)

#===============================
#Eliminar organismos duplicados
#===============================
organismos_modules = tabla_organismos[
    ["Organism_code", "Organism", "Module_ID"]
].drop_duplicates().reset_index(drop=True)

print(organismos_modules)
#====================================
#Resumen MODULES por organismo
#====================================
#Relacionamos los organismos con los Modules en los que aparece
organismos_resumen = (
    organismos_modules
    .groupby(["Organism_code", "Organism"], as_index=False)
    .agg(
        Modules=("Module_ID", lambda x: "; ".join(sorted(set(x))))
    )
)

print(organismos_resumen.head(20))
#Contamos cuantos modules comparte cada organismo 
organismos_resumen["Numero_modules"] = (
    organismos_resumen["Modules"]
    .str.split("; ")
    .str.len()
)
print(organismos_resumen)

#=========================================
# Caracterización de los perfiles
#=========================================
#Reducción del número de organismos candidatos
#Se tomará en cuenta
#1) Si hay organismos que repiten MODULES
#2) Diversidad taxonomica
#3) Diversidad de MODULES

organismos_detalle = organismos_modules.merge(
    modules_unicos[["Module_ID", "Pathways"]],
    on="Module_ID",
    how="left"
)

print(organismos_detalle.head())

perfiles = (
    organismos_detalle
    .groupby(["Organism_code", "Organism"])
    .agg(
        Modules=("Module_ID", lambda x: set(x)),
        Pathways=("Pathways", lambda x: set(
            pathway
            for grupo in x.dropna()
            for pathway in grupo.split("; ")
        ))
    )
    .reset_index()
)
perfiles["Numero_modules"] = perfiles["Modules"].apply(len)
perfiles["Numero_pathways"] = perfiles["Pathways"].apply(len)

perfiles["Genus"] = (
    perfiles["Organism"]
    .str.split()
    .str[0]
)

print("Organismos:", len(perfiles))
print("Géneros:", perfiles["Genus"].nunique())
print("Módulos:", modules_unicos["Module_ID"].nunique())
print("Pathways:", modules_unicos["Pathways"].nunique())

#============================
#Pathways representados
#============================

pathways_totales = set()

for pathways in modules_unicos["Pathways"].dropna():
    for pathway in pathways.split(";"):
        pathways_totales.add(pathway.strip())

print("Pathways representados:", len(pathways_totales))
print(sorted(pathways_totales))
#%%
#======================
#Pathways sin MODULES
#======================

pathways_sin_module = ["map00633", "map00642"]

for pathway_id in pathways_sin_module:
    
    informacion = informacion_pathways[pathway_id]
    
    print("\n" + "=" * 60)
    print(pathway_id)
    print("=" * 60)
    
    for linea in informacion.split("\n"):
        if linea.startswith("KO_PATHWAY"):
            print(linea)
#Descargar los KO de esas rutas 
informacion_ko_pathways = {}

for pathway_id in pathways_sin_module:
    
    ko_pathway_id = "ko" + pathway_id.replace("map", "")
    
    url = f"https://rest.kegg.jp/get/{ko_pathway_id}"
    
    respuesta = requests.get(url)
    
    print(
        pathway_id,
        "→",
        ko_pathway_id,
        "→ código:",
        respuesta.status_code
    )
    
    if respuesta.status_code == 200:
        informacion_ko_pathways[pathway_id] = respuesta.text
#Extraer KO
datos_ko_pathways = []

for pathway_id, informacion in informacion_ko_pathways.items():
    
    kos = sorted(set(
        re.findall(r"\bK\d{5}\b", informacion)
    ))
    
    for ko in kos:
        
        datos_ko_pathways.append({
            "Pathway_ID": pathway_id,
            "KO": ko
        })

tabla_ko_pathways = pd.DataFrame(datos_ko_pathways)

print(tabla_ko_pathways)

print(
    tabla_ko_pathways
    .groupby("Pathway_ID")["KO"]
    .nunique()
)

#================================================
#Unificar todos los KO de los MODULES de interés
#================================================

tabla_ko_completa = tabla_modules_ko[
    ["Module_ID", "Pathways", "KO", "Function"]
].copy()

tabla_ko_completa["Pathway_ID"] = (
    tabla_ko_completa["Pathways"]
    .str.split("; ")
)

tabla_ko_completa = tabla_ko_completa.explode(
    "Pathway_ID"
)

#KO de pathways sin MODULE
tabla_ko_sin_module = tabla_ko_pathways.copy()

tabla_ko_sin_module["Module_ID"] = pd.NA
tabla_ko_sin_module["Function"] = pd.NA

tabla_ko_sin_module = tabla_ko_sin_module[
    ["Module_ID", "Pathway_ID", "KO", "Function"]
]

#KO asociados a los MODULES
tabla_ko_modulos = tabla_ko_completa[
    ["Module_ID", "Pathway_ID", "KO", "Function"]
]

#Unir información
tabla_ko_total = pd.concat(
    [
        tabla_ko_modulos,
        tabla_ko_sin_module
    ],
    ignore_index=True
)

print(tabla_ko_total.head(20))

#====================
#Contar número de KO 
#====================

print(
    tabla_ko_total
    .groupby("Pathway_ID")["KO"]
    .nunique()
)
print(
    "\nNúmero total de KO únicos:",
    tabla_ko_total["KO"].nunique()
)

pathways_ko = sorted(
    tabla_ko_total["Pathway_ID"]
    .dropna()
    .unique()
)

print("Pathways con KO:")
for pathway in pathways_ko:
    print(pathway)
#%%
#===========================================
#Selección de candidatos según la diversidd
#===========================================
#Diversidad de: géneros, modules, pathways y perfiles funcional
perfiles["Perfil_modules"] = perfiles["Modules"].apply(
    lambda x: frozenset(x)
)

#Similitud entre perfiles funcionales
def similitud_jaccard(conjunto_a, conjunto_b):

    conjunto_a = set(conjunto_a)
    conjunto_b = set(conjunto_b)

    if not conjunto_a and not conjunto_b:
        return 1.0

    union = conjunto_a | conjunto_b

    if not union:
        return 0.0

    interseccion = conjunto_a & conjunto_b

    return len(interseccion) / len(union)

#Calcular Score de selección
def calcular_score(
    candidato,
    modulos_cubiertos,
    pathways_cubiertos,
    generos_seleccionados,
    perfiles_seleccionados
):
    
    # 1. Nuevos MODULES
    modulos_nuevos = (
        candidato["Modules"] - modulos_cubiertos
    )
    
    score_modulos = len(modulos_nuevos) * 10
    
    # 2. Nuevos PATHWAYS 
    pathways_nuevos = (
        candidato["Pathways"] - pathways_cubiertos
    )
    
    score_pathways = len(pathways_nuevos) * 6
    
    # 3. Número de MÓDULOS
    # Priorizar organismos con varios módulos,
    # pero sin permitir que este criterio domine.
    
    score_numero_modules = (
        min(candidato["Numero_modules"], 5) * 2
    )
    
    # 4. Diversidad taxonómica
    if candidato["Genus"] not in generos_seleccionados:
        score_genero = 8
    else:
        score_genero = 0
    
    # 5. Redundancia funcional
    
    penalizacion_redundancia = 0
    
    for perfil in perfiles_seleccionados:
        
        similitud = similitud_jaccard(
            candidato["Perfil_modules"],
            perfil
        )
        
        # Se penaliza cuando existe una similitud funcional
        
        if similitud >= 0.75:
            penalizacion_redundancia += (
                (similitud - 0.75) * 20
            )
    
    # SCORE FINAL
    score_total = (
        score_modulos
        + score_pathways
        + score_numero_modules
        + score_genero
        - penalizacion_redundancia
    )
    
    return score_total

#==================
#Selección inicial 
#==================

seleccionados = []

modulos_cubiertos = set()
pathways_cubiertos = set()
generos_seleccionados = set()
perfiles_seleccionados = []

for module_id in modules_unicos["Module_ID"]:
    
    # Organismos que poseen este módulo
    candidatos = perfiles[
        perfiles["Modules"].apply(
            lambda x: module_id in x
        )
    ].copy()
    
    # Eliminar organismos ya seleccionados
    candidatos = candidatos[
        ~candidatos["Organism_code"].isin(seleccionados)
    ]
    
    if candidatos.empty:
        continue
    
    # Ordenar según el número de MODULES y pathways
    candidatos = candidatos.sort_values(
        ["Numero_modules", "Numero_pathways"],
        ascending=False
    )
    
    candidato = candidatos.iloc[0]
    
    # Guardar organismo
    seleccionados.append(
        candidato["Organism_code"]
    )
    
    # Actualizar información
    modulos_cubiertos.update(
        candidato["Modules"]
    )
    
    pathways_cubiertos.update(
        candidato["Pathways"]
    )
    # Actualizar géneros  
    generos_seleccionados.add(
        candidato["Genus"]
    )
    # Guardar perfil
    perfiles_seleccionados.append(
        candidato["Perfil_modules"]
    )
    
print(
    "Organismos seleccionados inicialmente:",
    len(seleccionados)
)

print(
    "Módulos cubiertos:",
    len(modulos_cubiertos)
)

print(
    "Pathways cubiertos:",
    len(pathways_cubiertos)
)

print(
    "Géneros representados:",
    len(generos_seleccionados)
)

#=======================================
#Selección 150 organismos  candidatos
#=======================================

while len(seleccionados) < 150:
    
    candidatos = perfiles[
        ~perfiles["Organism_code"].isin(seleccionados)
    ].copy()
    
    if candidatos.empty:
        break
    
    
    scores = []
    for _, fila in candidatos.iterrows():
        score = calcular_score(
            fila,
            modulos_cubiertos,
            pathways_cubiertos,
            generos_seleccionados,
            perfiles_seleccionados)
        scores.append(score)

    candidatos["Score_seleccion"] = scores
    
    
#Elegir el mejor candidato según el score
    
    mejor = candidatos.sort_values(
        "Score_seleccion",
        ascending=False
    ).iloc[0]
    
#Guardar organismo
    
    seleccionados.append(
        mejor["Organism_code"]
    )
    
#Actualizar MODULES
    
    modulos_cubiertos.update(
        mejor["Modules"]
    )
    
#Actualizar pathways
    
    pathways_cubiertos.update(
        mejor["Pathways"]
    )
    
#Actualizar géneros
    
    generos_seleccionados.add(
        mejor["Genus"]
    )
    
#Actualizar perfiles
    
    perfiles_seleccionados.append(
        mejor["Perfil_modules"]
    )

#====================================
#Crear data frame con los candidatos
#====================================

candidatos_150 = perfiles[
    perfiles["Organism_code"].isin(seleccionados)
].copy()

candidatos_150 = candidatos_150.sort_values(
    ["Numero_modules", "Numero_pathways"],
    ascending=False
).reset_index(drop=True)

candidatos_150["Modules"] = candidatos_150["Modules"].apply(
    lambda x: "; ".join(sorted(x))
)

candidatos_150["Pathways"] = candidatos_150["Pathways"].apply(
    lambda x: "; ".join(sorted(x))
)

candidatos_150 = candidatos_150.drop(
    columns=["Perfil_modules"]
)

#Pathways objetivo
pathways_objetivo = [
    "map00361",
    "map00622",
    "map00623",
    "map00624",
    "map00626",
    "map00633",
    "map00642"
]

candidatos_150["Pathways_objetivo"] = "; ".join(
    pathways_objetivo
)

nombres_pathways = {
    "map00361": "Chlorocyclohexane and chlorobenzene degradation",
    "map00622": "Xylene degradation",
    "map00623": "Toluene degradation",
    "map00624": "Polycyclic aromatic hydrocarbon degradation",
    "map00626": "Naphthalene degradation",
    "map00633": "Nitrotoluene degradation",
    "map00642": "Ethylbenzene degradation"
}

candidatos_150["Pathways_objetivo"] = "; ".join(
    f"{pathway}: {nombres_pathways[pathway]}"
    for pathway in pathways_objetivo
)
#Comprobar los resultados
print("Número de candidatos:", len(candidatos_150))

print("Géneros representados:", candidatos_150["Genus"].nunique())

print("Módulos representados:",
    len(
        set().union(
            *[
                set(modules.split("; "))
                for modules in candidatos_150["Modules"]
            ]
        )
    )
)

print(
    "Pathways representados según Modules:",
    len(
        set().union(
            *[
                set(pathways.split("; "))
                for pathways in candidatos_150["Pathways"]
            ])))

print( "Pathways objetivo:", len(pathways_objetivo))

#Guardamos en un dataframe
candidatos_150.to_csv(
    Results_DIR /"candidatos_150.csv",
    index=False
)

candidatos_150.to_excel(
    Results_DIR /"candidatos_150.xlsx",
    index=False
)

#%%
#======================================================
#Obtener los genomas de los 150 organismos candidato
#======================================================

url_genomes = "https://rest.kegg.jp/list/genome"

respuesta_genomes = requests.get(url_genomes)

#Diccionario de genomas de KEGG
genome_dict = {}

for linea in respuesta_genomes.text.strip().splitlines():
    
    partes = linea.split("\t", maxsplit=1)
    
    if len(partes) == 2:
        
        genome_id = partes[0]
        informacion = partes[1]
        partes_organismo = informacion.split(";", maxsplit=1)
        
        if len(partes_organismo) == 2:
            
            organism_code = partes_organismo[0].strip()
            organism_name = partes_organismo[1].strip()
            
            genome_dict[organism_code] = {
                "Genome_ID": genome_id,
                "Organism": organism_name
            }

print("Número de organismos encontrados en KEGG:", len(genome_dict))

#Asociar el Genome_ID a los candidatos

datos_genomas = []

for _, fila in candidatos_150.iterrows():
    
    organism_code = fila["Organism_code"]
    organism = fila["Organism"]
    
    if organism_code in genome_dict:
        
        datos_genomas.append({
            "Organism_code": organism_code,
            "Organism": organism,
            "Genome_ID": genome_dict[organism_code]["Genome_ID"],
            "KEGG_Organism": genome_dict[organism_code]["Organism"]
        })
    
    else:
        
        datos_genomas.append({
            "Organism_code": organism_code,
            "Organism": organism,
            "Genome_ID": pd.NA,
            "KEGG_Organism": pd.NA
        })

tabla_genomas = pd.DataFrame(datos_genomas)

print("Candidatos seleccionados:", len(candidatos_150))

print("Genomas encontrados en KEGG:",
    tabla_genomas["Genome_ID"].notna().sum())

print("Genomas no encontrados:",
    tabla_genomas["Genome_ID"].isna().sum())

#Identificar organismos sin genoma 
genomas_faltantes = tabla_genomas[
    tabla_genomas["Genome_ID"].isna()
][
    ["Organism_code", "Organism"]
]

print(genomas_faltantes)

#Buscar directamente en KEGG estos 3 organismos 
organismos_faltantes = [
    "Pseudomonas putida",
    "Pseudomonas sp.",
    "Thauera aromatica"]

for nombre in organismos_faltantes:
    
    print(nombre)
    
    coincidencias = []
    
    for linea in respuesta_genomes.text.splitlines():
        
        if "\t" not in linea:
            continue
        
        genome_id, informacion = linea.split("\t", maxsplit=1)
        
        partes = informacion.split(";", maxsplit=1)
        
        if len(partes) != 2:
            continue
        
        organism_code = partes[0].strip()
        organism_name = partes[1].strip()
        
        if nombre.lower() in organism_name.lower():
            coincidencias.append(
                (
                    genome_id,
                    organism_code,
                    organism_name
                )
            )
    
    for coincidencia in coincidencias:
        print(coincidencia)
#Resolver los genomas de organismos faltantes
               
for module_id, informacion in informacion_modules.items():
    
    for linea in informacion.split("\n"):
        
        if linea.strip().startswith("COMPLETE"):
            if any(codigo in linea
                   for codigo in ["303", "306", "59405"]
            ):
                
                print(
                    module_id,
                    "→",
                    linea.strip()
                )
#Corregir el organismo Thauera aromatica
tabla_genomas.loc[
    tabla_genomas["Organism_code"] == "59405",
    "Genome_ID"
] = "T05381"

tabla_genomas.loc[
    tabla_genomas["Organism_code"] == "59405",
    "KEGG_Organism"
] = "Thauera aromatica"

tabla_genomas.loc[
    tabla_genomas["Organism_code"] == "59405",
    "Organism_code"
] = "tak"

#Eliminar organismos sin genoma disponible

codigos_eliminar = ["303", "306"]

candidatos_148 = candidatos_150[
    ~candidatos_150["Organism_code"].isin(codigos_eliminar)
].copy()

tabla_genomas_148 = tabla_genomas[
    ~tabla_genomas["Organism_code"].isin(codigos_eliminar)
].copy()

print("Candidatos:", len(candidatos_148))

print("Registros en la tabla de genomas:", 
      len(tabla_genomas_148))

print("Genomas encontrados:",
    tabla_genomas_148["Genome_ID"].notna().sum()
)

print("Genomas faltantes:",
    tabla_genomas_148["Genome_ID"].isna().sum()
)

#====================
#Guardar candidatos
#====================

candidatos_148.to_csv(Results_DIR /"candidatos_148.csv",
    index=False)

candidatos_148.to_excel(Results_DIR /"candidatos_148.xlsx",
    index=False)
#Comprobamos la información 
print(tabla_genomas_148["Genome_ID"].nunique())
#%%
#======================
#Creación de la matriz 
#======================
kos_objetivo = sorted(
    tabla_ko_total["KO"]
    .dropna()
    .unique()
)

print("Número de KO objetivo:", len(kos_objetivo))
print(kos_objetivo)

def obtener_kos_organismo(organism_code):
    
    url_matrix = f"https://rest.kegg.jp/link/ko/{organism_code}"
    
    respuesta = requests.get(url_matrix)
    
    if respuesta.status_code != 200:
        print("Error en:", organism_code)
        return set()    
    kos = set()
    
    for linea in respuesta.text.strip().splitlines():
        
        partes = linea.split("\t")
        
        if len(partes) == 2:
            
            ko = partes[1].replace("ko:", "").strip()
            
            kos.add(ko)
    
    return kos

kos_por_organismo = {}

for organism_code in tabla_genomas_148["Organism_code"]:
    
    print("Procesando:", organism_code)
    
    kos_por_organismo[organism_code] = (
        obtener_kos_organismo(organism_code)
    )
#  
set_kos_objetivo = set(kos_objetivo)

kos_objetivo_por_organismo = {}

for organism_code, kos in kos_por_organismo.items():
    
    kos_objetivo_por_organismo[organism_code] = (
        kos.intersection(set_kos_objetivo)
    )

matriz_features = pd.DataFrame(
    0,
    index=tabla_genomas_148["Organism_code"],
    columns=kos_objetivo
)

for organism_code, kos in kos_objetivo_por_organismo.items():
    
    for ko in kos:
        
        matriz_features.loc[
            organism_code,
            ko
        ] = 1
matriz_features = matriz_features.reset_index()

matriz_features = matriz_features.rename(
    columns={
        "index": "Organism_code"
    }
)
print(
    "Organismos:",
    matriz_features.shape[0]
)

print(
    "Features KO:",
    matriz_features.shape[1] - 1
)

print(
    matriz_features.head()
)

#Comprobamos
#Número de organismos que presentan cada KO
frecuencia_ko = matriz_features.drop(
    columns=["Organism_code"]
).sum(axis=0)

print(frecuencia_ko.sort_values())

print("KO presentes en 0 organismos:")
print(
    frecuencia_ko[frecuencia_ko == 0]
)
print("KO presentes en 1 organismo:")
print(
    frecuencia_ko[frecuencia_ko == 1]
)
print("KO presentes en más del 95% de los organismos:")
print(
    frecuencia_ko[
        frecuencia_ko >= 0.95 * len(matriz_features)
    ]
)
#Construir la matriz depurada
kos_eliminar = frecuencia_ko[
    frecuencia_ko == 0
].index.tolist()

print("KO que se eliminarán:", len(kos_eliminar))
print(kos_eliminar)

matriz_features_limpia = matriz_features.drop(
    columns=kos_eliminar
)

print(
    "Organismos:",
    matriz_features_limpia.shape[0]
)

print(
    "Features KO:",
    matriz_features_limpia.shape[1] - 1
)

#Guardamos la información de los ko eliminados
ko_eliminados = pd.DataFrame({
    "KO": kos_eliminar,
    "Motivo": "Ausente en los 148 organismos de referencia"
})

print(ko_eliminados)

#%%
# =======================================
# Matrices de relación MODULE Y PATHWAY
# =======================================

#Matriz de relacion Organismo x KO

#Crear una copia de la matriz
matriz_ko = matriz_features.copy()

# Asegurar que el índice sea Organism_code
matriz_ko = matriz_ko.set_index("Organism_code")

# Asegurar que todos los KO sean numéricos
matriz_ko = matriz_ko.astype(int)

print("Organismos:", matriz_ko.shape[0])

print("KO:", matriz_ko.shape[1])

#Definir los MODULES
modules = sorted(
    tabla_ko_total[
        "Module_ID"
    ]
    .dropna()
    .unique()
)

print("\nMODULES encontrados:", len(modules))

print(modules)

#Crear Matrices MODULE

matriz_module_presencia = pd.DataFrame(
    0,
    index=matriz_ko.index,
    columns=modules,
    dtype=int
)


matriz_module_cobertura = pd.DataFrame(
    0.0,
    index=matriz_ko.index,
    columns=modules,
    dtype=float
)

#Calcular MODULE por MODULE

for module in modules:

#KO definidos para este módulo
    kos_module = (
        tabla_ko_total.loc[
            tabla_ko_total["Module_ID"] == module,
            "KO"
        ]
        .dropna()
        .unique()
        .tolist()
    )
    if len(kos_module) == 0:
        continue


#KO que realmente existen como columnas en la matriz de organismos
    kos_validos = [
        ko
        for ko in kos_module
        if ko in matriz_ko.columns
    ]


    if len(kos_validos) == 0:

        continue


#Presencia

    matriz_module_presencia[
        module
    ] = (
        matriz_ko[
            kos_validos
        ]
        .sum(axis=1)
        > 0
    ).astype(int)


#Cobertura

    numero_ko_requeridos = len(
        kos_module
    )


    numero_ko_presentes = (
        matriz_ko[
            kos_validos
        ]
        .sum(axis=1)
    )


    matriz_module_cobertura[
        module
    ] = (
        numero_ko_presentes /
        numero_ko_requeridos
    )

#Definir los pathways

pathways = sorted(
    tabla_ko_total[
        "Pathway_ID"
    ]
    .dropna()
    .unique()
)


print("Pathways encontrados:", len(pathways))

print(pathways)

#Creación de matrices pathways

matriz_pathway_presencia = pd.DataFrame(
    0,
    index=matriz_ko.index,
    columns=pathways,
    dtype=int
)


matriz_pathway_cobertura = pd.DataFrame(
    0.0,
    index=matriz_ko.index,
    columns=pathways,
    dtype=float
)


#Calcular pathways 

for pathway in pathways:

#KO definidos para este pathway
    kos_pathway = (
        tabla_ko_total.loc[
            tabla_ko_total["Pathway_ID"] == pathway,
            "KO"
        ]
        .dropna()
        .unique()
        .tolist()
    )

    if len(kos_pathway) == 0:
        continue


#KO disponibles en la matriz
    kos_validos = [
        ko
        for ko in kos_pathway
        if ko in matriz_ko.columns
    ]


    if len(kos_validos) == 0:
        continue

#Presencia

    matriz_pathway_presencia[
        pathway
    ] = (
        matriz_ko[
            kos_validos
        ]
        .sum(axis=1)
        > 0
    ).astype(int)


#Cobertura

    numero_ko_requeridos = len(
        kos_pathway
    )


    numero_ko_presentes = (
        matriz_ko[
            kos_validos
        ]
        .sum(axis=1)
    )


    matriz_pathway_cobertura[
        pathway
    ] = (
        numero_ko_presentes /
        numero_ko_requeridos
    )

#Verificar las matrices pathways

print(matriz_pathway_presencia.shape)

print(matriz_pathway_presencia.head())

print(matriz_pathway_cobertura.shape)

print(matriz_pathway_cobertura.head())


#Resumen resultados cobertura

print("\nCobertura de Modules:", 
      matriz_module_cobertura.describe())

print("\nCobertura de Pathways:",
      matriz_pathway_cobertura.describe())

#Guardar los resultados

matriz_module_presencia.to_csv(
    Results_DIR / "matriz_organismo_Module_presencia.csv",
    encoding="utf-8-sig"
)

matriz_module_cobertura.to_csv(
    Results_DIR / "matriz_organismo_Module_cobertura.csv",
    encoding="utf-8-sig"
)

matriz_pathway_presencia.to_csv(
    Results_DIR / "matriz_organismo_Pathway_presencia.csv",
    encoding="utf-8-sig"
)

matriz_pathway_cobertura.to_csv(
    Results_DIR / "matriz_organismo_Pathway_cobertura.csv",
    encoding="utf-8-sig"
)

#%%
#Creación de una Tabla funcional entre organismos, pathways, MODULES y KO

matriz_ko = matriz_features.copy()

matriz_ko = (
    matriz_ko
    .set_index("Organism_code")
)

matriz_ko_larga = (
    matriz_ko
    .reset_index()
    .melt(
        id_vars="Organism_code",
        var_name="KO",
        value_name="Present"
    )
)


matriz_ko_larga["Present"] = (
    matriz_ko_larga["Present"]
    .astype(int)
)


print(
    "Registros organismo × KO:",
    len(matriz_ko_larga)
)

#Relación funcional pathway, MODULE y KO

relacion_funcional = (
    tabla_ko_total[
        [
            "Pathway_ID",
            "Module_ID",
            "KO"
        ]
    ]
    .drop_duplicates()
    .copy()
)


print(
    "Registros únicos:",
    len(relacion_funcional)
)

#Analizar los organismos y los KO de referencia

tabla_funcional = (
    relacion_funcional
    .merge(
        matriz_ko_larga,
        on="KO",
        how="left"
    )
)

#Agregar el nombre de los organismos

nombres_organismos = (
    candidatos_148[
        [
            "Organism_code",
            "Organism"
        ]
    ]
    .drop_duplicates(
        subset=["Organism_code"]
    )
    .copy()
)

nombres_organismos["Organism_code"] = (
    nombres_organismos["Organism_code"]
    .astype(str)
    .str.strip()
    .str.lower()
)


tabla_funcional["Organism_code"] = (
    tabla_funcional["Organism_code"]
    .astype(str)
    .str.strip()
    .str.lower()
)


tabla_funcional = (
    tabla_funcional
    .merge(
        nombres_organismos,
        on="Organism_code",
        how="left"
    )
)

#Ordenar

tabla_funcional = tabla_funcional[
    [
        "Organism_code",
        "Organism",
        "Pathway_ID",
        "Module_ID",
        "KO",
        "Present"
    ]
]

tabla_funcional = (
    tabla_funcional
    .sort_values(
        [
            "Organism_code",
            "Pathway_ID",
            "Module_ID",
            "KO"
        ]
    )
    .reset_index(drop=True)
)

#Resultados

print(tabla_funcional.head(20))

#Comprobar si existen nombres faltantes
nombres_faltantes = (
    tabla_funcional["Organism"]
    .isna()
    .sum()
)

print("Registros sin nombre:", nombres_faltantes)


if nombres_faltantes > 0:
    print(
        tabla_funcional.loc[
            tabla_funcional["Organism"].isna(),
            "Organism_code"
        ]
        .unique()
    )


#Almacenar

archivo_salida = (
    Results_DIR /
    "tabla_funcional_organismo_pathway_module_KO.csv"
)

tabla_funcional.to_csv(
    archivo_salida,
    index=False,
    encoding="utf-8-sig"
)

print(archivo_salida)

#%%
#===================================================
#Normalización de los identificadores de organsimos
#===================================================
#Crear tabla de equivalencias a partir de tabla_genomas_148

equivalencias_organismos = (
    tabla_genomas_148[
        [
            "Organism_code",
            "Organism",
            "Genome_ID",
            "KEGG_Organism"
        ]
    ]
    .drop_duplicates()
    .copy()
)

#Crear diccionario

organismo_a_kegg = dict(
    zip(
        equivalencias_organismos["Organism"],
        equivalencias_organismos["Organism_code"]
    )
)


#%%
# ============================================================
#Creación de matriz de realcion KO, RGANISMO, GENE_ID y GENE
# ============================================================

kos_objetivo = (
    tabla_ko_total["KO"]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
    .tolist()
)


#ORGANISMOS KEGG DE LOS 148 CANDIDATOS

organismos_148 = (
    tabla_genomas_148[
        [
            "Organism_code",
            "Organism",
            "Genome_ID"
        ]
    ]
    .drop_duplicates()
    .copy()
)


#Normalizar códigos

organismos_148["Organism_code"] = (
    organismos_148["Organism_code"]
    .astype(str)
    .str.strip()
    .str.lower()
)


codigos_148 = set(
    organismos_148["Organism_code"]
)

#Consultar los KO

resultados_genes = []


for i, ko in enumerate(kos_objetivo, start=1):

    url_ko = (
        f"https://rest.kegg.jp/get/{ko}"
    )
    try:

        respuesta = requests.get(
            url_ko,
            timeout=60
        )
    except requests.RequestException as error:

        print(
            "  Error de conexión:",
            error
        )
        continue
    if respuesta.status_code != 200:

        print(
            "  Error HTTP:",
            respuesta.status_code
        )
        continue

    lineas = (
        respuesta.text
        .splitlines()
    )
    leyendo_genes = False
    encontrados_ko = 0
    for linea in lineas:

#Inicio de GENES

        if linea.startswith("GENES"):
            leyendo_genes = True
            contenido = (
                linea[5:]
                .strip()
            )
        elif leyendo_genes:
            if linea.startswith(" "):
                contenido = (
                    linea.strip()
                )
            else:
                break
        else:
            continue

# Procesar línea

        if ":" not in contenido:
            continue
        codigo_organismo = (
            contenido
            .split(":")[0]
            .strip()
            .lower()
        )
        if codigo_organismo not in codigos_148:
            continue
        resto = (
            contenido
            .split(":", 1)[1]
            .strip()
        )
        if not resto:
            continue
        gene_id = (
            resto
            .split()[0]
        )
        resultados_genes.append(
            {
                "Organism_code":
                    codigo_organismo,
                "KO":
                    ko,
                "Gene_ID":
                    gene_id
            }
        )
        encontrados_ko += 1

#Crear el dataframe

tabla_organismo_KO_gene = pd.DataFrame(
    resultados_genes
)

#Eliminar los duplicados

tabla_organismo_KO_gene = (
    tabla_organismo_KO_gene
    .drop_duplicates()
    .reset_index(drop=True)
)

#Añadir la información del organismo

tabla_organismo_KO_gene = (
    tabla_organismo_KO_gene
    .merge(
        organismos_148,
        on="Organism_code",
        how="left"
    )
)

#===========================================
#Construcción de la matriz Organismo x Gene
#===========================================

#Realizar una copia de la tabla

tabla_genes = (
    tabla_organismo_KO_gene
    .copy()
)

#Normalizar texto

for columna in [
    "Organism_code",
    "KO",
    "Gene_ID"
]:

    tabla_genes[columna] = (
        tabla_genes[columna]
        .astype(str)
        .str.strip()
    )

#Eliminar duplicados

tabla_genes = (
    tabla_genes[
        [
            "Organism_code",
            "Organism",
            "Genome_ID",
            "KO",
            "Gene_ID"
        ]
    ]
    .drop_duplicates()
    .reset_index(drop=True)
)


#Crear matriz Organismo x GENE_ID

tabla_genes["Present"] = 1

matriz_organismo_gene = (
    tabla_genes
    .pivot_table(
        index="Organism_code",
        columns="Gene_ID",
        values="Present",
        aggfunc="max",
        fill_value=0
    )
)

matriz_organismo_gene = (
    matriz_organismo_gene
    .astype(int)
)

#Analizar los 148 organismos candidatos

codigos_148 = (
    tabla_genomas_148[
        "Organism_code"
    ]
    .astype(str)
    .str.strip()
    .str.lower()
    .unique()
)


matriz_organismo_gene = (
    matriz_organismo_gene
    .reindex(
        codigos_148,
        fill_value=0
    )
)

#Frecuencia de los genes

frecuencia_genes = (
    matriz_organismo_gene
    .sum(axis=0)
    .sort_values(
        ascending=False
    )
)

#Relación del GENE_ID y KO

tabla_gene_KO = (
    tabla_genes[
        [
            "Gene_ID",
            "KO"
        ]
    ]
    .drop_duplicates()
    .sort_values(
        [
            "KO",
            "Gene_ID"
        ]
    )
)

#Relación entre GENE, KO, MODULE y PATHWAY

tabla_gene_funcional = (
    tabla_genes[
        [
            "Organism_code",
            "Organism",
            "Genome_ID",
            "Gene_ID",
            "KO"
        ]
    ]
    .merge(
        tabla_ko_total[
            [
                "Module_ID",
                "Pathway_ID",
                "KO"
            ]
        ]
        .drop_duplicates(),
        on="KO",
        how="left"
    )
    .drop_duplicates()
    .sort_values(
        [
            "Organism_code",
            "Pathway_ID",
            "Module_ID",
            "KO",
            "Gene_ID"
        ]
    )
    .reset_index(drop=True)
)

print(tabla_gene_funcional.head(20))

#Guardar los resultados

archivo_matriz_gene = (
    Results_DIR / "matriz_organismo_x_gene.csv"
)

archivo_tabla_gene = (
    Results_DIR / "tabla_organismo_KO_Gene.csv"
)

archivo_gene_funcional = (
    Results_DIR / "tabla_gene_KO_module_pathway.csv"
)

matriz_organismo_gene.to_csv(
    archivo_matriz_gene,
    encoding="utf-8-sig"
)


tabla_genes.to_csv(
    archivo_tabla_gene,
    index=False,
    encoding="utf-8-sig"
)


tabla_gene_funcional.to_csv(
    archivo_gene_funcional,
    index=False,
    encoding="utf-8-sig"
)

# ==========================================
#Calcular el número de genes por organismo
# ==========================================

genes_por_organismo = (
    tabla_organismo_KO_gene
    .groupby(
        [
            "Organism_code",
            "Organism",
            "Genome_ID"
        ]
    )
    .agg(
        Genes_objetivo=(
            "Gene_ID",
            "nunique"
        ),

        KO_presentes=(
            "KO",
            "nunique"
        )
    )
    .reset_index()
)

organismos_base = (
    tabla_genomas_148[
        [
            "Organism_code",
            "Organism",
            "Genome_ID"
        ]
    ]
    .drop_duplicates()
)


genes_por_organismo = (
    organismos_base
    .merge(
        genes_por_organismo,
        on=[
            "Organism_code",
            "Organism",
            "Genome_ID"
        ],
        how="left"
    )
)


#Organismos sin ningún KO objetivo

genes_por_organismo[
    "Genes_objetivo"
] = (
    genes_por_organismo[
        "Genes_objetivo"
    ]
    .fillna(0)
    .astype(int)
)


genes_por_organismo[
    "KO_presentes"
] = (
    genes_por_organismo[
        "KO_presentes"
    ]
    .fillna(0)
    .astype(int)
)

#Ordenar

genes_por_organismo = (
    genes_por_organismo
    .sort_values(
        "Genes_objetivo",
        ascending=False
    )
    .reset_index(drop=True)
)
print(genes_por_organismo.head(20))


#Estadísticas

print(
    genes_por_organismo[
        "Genes_objetivo"
    ].describe()
)


print(
    genes_por_organismo[
        "KO_presentes"
    ].describe()
)


# 6. EXTREMOS
print(
    genes_por_organismo[
        [
            "Organism_code",
            "Organism",
            "Genes_objetivo",
            "KO_presentes"
        ]
    ]
    .head(10)
)

print(
    genes_por_organismo[
        [
            "Organism_code",
            "Organism",
            "Genes_objetivo",
            "KO_presentes"
        ]
    ]
    .sort_values(
        "Genes_objetivo"
    )
    .head(10)
)


#Candidatos sin los genes objetivo

sin_genes = (
    genes_por_organismo[
        "Genes_objetivo"
    ] == 0
).sum()


print(
    "\nOrganismos sin genes de los KO objetivo:",
    sin_genes
)


#Guardar

archivo = (Results_DIR / "resumen_genes_por_organismo.csv")

genes_por_organismo.to_csv(
    archivo,
    index=False,
    encoding="utf-8-sig"
)

print(archivo)

#%%
# ============================================================
#Filtrar los organismos según el número de genes objetivo
# ============================================================

umbral_genes = 8

organismos_filtrados = (
    genes_por_organismo[
        genes_por_organismo["Genes_objetivo"] >= umbral_genes
    ]
    .copy()
    .reset_index(drop=True)
)

print(
    genes_por_organismo[
        genes_por_organismo["Genes_objetivo"] < umbral_genes
    ][
        [
            "Organism_code",
            "Organism",
            "Genes_objetivo",
            "KO_presentes"
        ]
    ]
)
      
# ==============================
# Conjunto final de organismos
# ==============================

UMBRAL_GENES = 8

codigos_finales = set(
    organismos_filtrados[
        "Organism_code"
    ]
    .astype(str)
    .str.strip()
    .str.lower()
)

#======================================
# Tabla final de organismos candidatos
#======================================

tabla_organismos_final = (
    tabla_genomas_148[
        tabla_genomas_148[
            "Organism_code"
        ]
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(codigos_finales)
    ]
    .copy()
    .reset_index(drop=True)
)

print(
    tabla_organismos_final[
        "Organism_code"
    ]
    .duplicated()
    .sum()
)

#Guardar
archivo_final = (Results_DIR / "organismos_finales_141.csv")

tabla_organismos_final.to_csv(
    archivo_final,
    index=False,
    encoding="utf-8-sig"
)

print(archivo_final)      

#%%
# =================================================
#Matriz final y depurada de organismos candidatos
# =================================================

codigos_finales = set(
    tabla_organismos_final[
        "Organism_code"
    ]
    .astype(str)
    .str.strip()
    .str.lower()
)


# ========================================
#Filtrar tabla de organismos x KO x gene
# ========================================

tabla_genes_final = (
    tabla_organismo_KO_gene[
        tabla_organismo_KO_gene[
            "Organism_code"
        ]
        .astype(str)
        .str.strip()
        .str.lower()
        .isin(codigos_finales)
    ]
    .copy()
)


print(tabla_genes_final.head(20))

# ======================
#Matriz organismo × KO
# ======================

matriz_KO_final = (
    tabla_genes_final
    .assign(Present=1)
    .pivot_table(
        index="Organism_code",
        columns="KO",
        values="Present",
        aggfunc="max",
        fill_value=0
    )
)

# Asegurar los KO objetivo

kos_objetivo = (
    tabla_ko_total[
        "KO"
    ]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
)


matriz_KO_final = (
    matriz_KO_final
    .reindex(
        index=sorted(codigos_finales),
        columns=sorted(kos_objetivo),
        fill_value=0
    )
    .astype(int)
)

# ============================================================
#Tabla KO x MODULE x pathway
# ============================================================

tabla_funcional = (
    tabla_ko_total[
        [
            "Module_ID",
            "Pathway_ID",
            "KO"
        ]
    ]
    .drop_duplicates()
)


# ================================
#Crear matriz organismo x MODULE
# ================================

tabla_organismo_module = (
    matriz_KO_final
    .reset_index()
    .melt(
        id_vars="Organism_code",
        var_name="KO",
        value_name="Present"
    )
    .merge(
        tabla_funcional,
        on="KO",
        how="left"
    )
)


matriz_module_final = (
    tabla_organismo_module[
        tabla_organismo_module[
            "Present"
        ] == 1
    ]
    .assign(Value=1)
    .pivot_table(
        index="Organism_code",
        columns="Module_ID",
        values="Value",
        aggfunc="max",
        fill_value=0
    )
)


modules_objetivo = (
    tabla_ko_total[
        "Module_ID"
    ]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
)

matriz_module_final = (
    matriz_module_final
    .reindex(
        index=sorted(codigos_finales),
        columns=sorted(modules_objetivo),
        fill_value=0
    )
    .astype(int)
)

print(matriz_module_final.shape)

# ================================
#Analizar la cobertura del MODULE 
# ================================

#La cobertura corresponde a: KO presentes del MODULE y
#KO requeridos por el módulo

conteo_KO_module = (
    tabla_funcional
    .groupby(
        "Module_ID"
    )[
        "KO"
    ]
    .nunique()
)


tabla_module_cobertura = (
    tabla_organismo_module
    .groupby(
        [
            "Organism_code",
            "Module_ID"
        ]
    )
    .agg(
        KO_presentes=(
            "Present",
            "sum"
        )
    )
    .reset_index()
)


tabla_module_cobertura[
    "KO_totales"
] = (
    tabla_module_cobertura[
        "Module_ID"
    ]
    .map(
        conteo_KO_module
    )
)


tabla_module_cobertura[
    "Coverage"
] = (
    tabla_module_cobertura[
        "KO_presentes"
    ]
    /
    tabla_module_cobertura[
        "KO_totales"
    ]
)


matriz_module_cobertura = (
    tabla_module_cobertura
    .pivot(
        index="Organism_code",
        columns="Module_ID",
        values="Coverage"
    )
    .fillna(0)
)


matriz_module_cobertura = (
    matriz_module_cobertura
    .reindex(
        index=sorted(codigos_finales),
        columns=sorted(modules_objetivo),
        fill_value=0
    )
)

print(matriz_module_cobertura.shape)

# =================================
#Crear matriz organismo × pathway
# =================================

tabla_organismo_pathway = (
    matriz_KO_final
    .reset_index()
    .melt(
        id_vars="Organism_code",
        var_name="KO",
        value_name="Present"
    )
    .merge(
        tabla_funcional,
        on="KO",
        how="left"
    )
)

matriz_pathway_final = (
    tabla_organismo_pathway[
        tabla_organismo_pathway[
            "Present"
        ] == 1
    ]
    .assign(Value=1)
    .pivot_table(
        index="Organism_code",
        columns="Pathway_ID",
        values="Value",
        aggfunc="max",
        fill_value=0
    )
)


pathways_objetivo = (
    tabla_ko_total[
        "Pathway_ID"
    ]
    .dropna()
    .astype(str)
    .str.strip()
    .unique()
)


matriz_pathway_final = (
    matriz_pathway_final
    .reindex(
        index=sorted(codigos_finales),
        columns=sorted(pathways_objetivo),
        fill_value=0
    )
    .astype(int)
)

print(matriz_pathway_final.shape)


# ==============================
#Analizar cobertura del pathway
# ==============================

conteo_KO_pathway = (
    tabla_funcional
    .groupby(
        "Pathway_ID"
    )[
        "KO"
    ]
    .nunique()
)


tabla_pathway_cobertura = (
    tabla_organismo_pathway
    .groupby(
        [
            "Organism_code",
            "Pathway_ID"
        ]
    )
    .agg(
        KO_presentes=(
            "Present",
            "sum"
        )
    )
    .reset_index()
)


tabla_pathway_cobertura[
    "KO_totales"
] = (
    tabla_pathway_cobertura[
        "Pathway_ID"
    ]
    .map(
        conteo_KO_pathway
    )
)


tabla_pathway_cobertura[
    "Coverage"
] = (
    tabla_pathway_cobertura[
        "KO_presentes"
    ]
    /
    tabla_pathway_cobertura[
        "KO_totales"
    ]
)


matriz_pathway_cobertura = (
    tabla_pathway_cobertura
    .pivot(
        index="Organism_code",
        columns="Pathway_ID",
        values="Coverage"
    )
    .fillna(0)
)


matriz_pathway_cobertura = (
    matriz_pathway_cobertura
    .reindex(
        index=sorted(codigos_finales),
        columns=sorted(pathways_objetivo),
        fill_value=0
    )
)

print(matriz_pathway_cobertura.shape)
    

# =============================
#Crear matriz organimos x gene
# =============================

matriz_gene_final = (
    tabla_genes_final
    .assign(Present=1)
    .pivot_table(
        index="Organism_code",
        columns="Gene_ID",
        values="Present",
        aggfunc="max",
        fill_value=0
    )
    .astype(int)
)


matriz_gene_final = (
    matriz_gene_final
    .reindex(
        index=sorted(codigos_finales),
        fill_value=0
    )
)

print(matriz_gene_final.shape)

#Guardar

matriz_KO_final.to_csv(
    Results_DIR / "matriz_final_organismo_x_KO.csv",
    encoding="utf-8-sig"
)

matriz_module_final.to_csv(
    Results_DIR /
    "matriz_final_organismo_x_Module.csv",
    encoding="utf-8-sig"
)

matriz_module_cobertura.to_csv(
    Results_DIR /
    "matriz_final_organismo_x_Module_cobertura.csv",
    encoding="utf-8-sig"
)

matriz_pathway_final.to_csv(
    Results_DIR /
    "matriz_final_organismo_x_Pathway.csv",
    encoding="utf-8-sig"
)

matriz_pathway_cobertura.to_csv(
    Results_DIR /
    "matriz_final_organismo_x_Pathway_cobertura.csv",
    encoding="utf-8-sig"
)

matriz_gene_final.to_csv(
    Results_DIR /
    "matriz_final_organismo_x_Gene.csv",
    encoding="utf-8-sig"
)

#%%
#===============================================
#Relación funcional Gene, KO, MODULE y pathway 
#===============================================

tabla_gene_KO = (
    tabla_genes_final[
        [
            "Organism_code",
            "Organism",
            "Genome_ID",
            "Gene_ID",
            "KO"
        ]
    ]
    .drop_duplicates()
    .copy()
)

# =============================
#Relación KO, MODULE y pathway
# =============================

tabla_KO_Module_Pathway = (
    tabla_ko_total[
        [
            "KO",
            "Module_ID",
            "Pathway_ID"
        ]
    ]
    .drop_duplicates()
    .copy()
)


#Normalizar

for columna in [
    "KO",
    "Module_ID",
    "Pathway_ID"
]:

    tabla_KO_Module_Pathway[columna] = (
        tabla_KO_Module_Pathway[columna]
        .astype(str)
        .str.strip()
    )

print(
    "Registros:",
    len(tabla_KO_Module_Pathway)
)

print(
    "KO:",
    tabla_KO_Module_Pathway["KO"].nunique()
)

print(
    "Modules:",
    tabla_KO_Module_Pathway["Module_ID"].nunique()
)

print(
    "Pathways:",
    tabla_KO_Module_Pathway["Pathway_ID"].nunique()
)


#Unir las matrices GENE x KO, con KO x MODULE x pathway

tabla_cadena_funcional = (
    tabla_gene_KO
    .merge(
        tabla_KO_Module_Pathway,
        on="KO",
        how="left"
    )
    .drop_duplicates()
    .sort_values(
        [
            "Organism_code",
            "Pathway_ID",
            "Module_ID",
            "KO",
            "Gene_ID"
        ]
    )
    .reset_index(drop=True)
)

print(tabla_cadena_funcional.head(20))

#Comprobar si hay KO sin MODULE o pathway

sin_module = (
    tabla_cadena_funcional[
        "Module_ID"
    ]
    .isna()
    .sum()
)

sin_pathway = (
    tabla_cadena_funcional[
        "Pathway_ID"
    ]
    .isna()
    .sum()
)

resumen_funcional = (
    tabla_cadena_funcional
    .groupby(
        [
            "Organism_code",
            "Organism",
            "Genome_ID"
        ]
    )
    .agg(
        Genes_objetivo=(
            "Gene_ID",
            "nunique"
        ),

        KO_presentes=(
            "KO",
            "nunique"
        ),

        Modules_presentes=(
            "Module_ID",
            "nunique"
        ),

        Pathways_presentes=(
            "Pathway_ID",
            "nunique"
        )
    )
    .reset_index()
)

#Guardar

archivo_cadena = (
    Results_DIR /
    "cadena_gene_KO_module_pathway.csv"
)

archivo_resumen = (
    Results_DIR /
    "resumen_funcional_141.csv"
)

archivo_matriz_resumen = (
    Results_DIR /
    "matriz_resumen_funcional_141.csv"
)

tabla_cadena_funcional.to_csv(
    archivo_cadena,
    index=False,
    encoding="utf-8-sig"
)


#==================================
#Evaluación de features por nivel
#==================================

#Analizar la frecuencia de los genes

frecuencia_gene = (
    matriz_gene_final
    .sum(axis=0)
    .sort_values(ascending=False)
)

#Analizar la frecuencia de los KO

frecuencia_KO = (
    matriz_KO_final
    .sum(axis=0)
    .sort_values(ascending=False)
)

#Analizar la frecuencia de MODULE

frecuencia_module = (
    matriz_module_final
    .sum(axis=0)
    .sort_values(ascending=False)
)


#Analizar la frecuencia de los PATHWAY

frecuencia_pathway = (
    matriz_pathway_final
    .sum(axis=0)
    .sort_values(ascending=False)
)

#Resumen 

def resumen_features(
    frecuencia,
    numero_organismos,
    nombre
):

    print(
        "Número de features:",
        len(frecuencia)
    )

    print(
        "Features presentes en al menos un organismo:",
        (frecuencia > 0).sum()
    )

    print(
        "Features presentes en un solo organismo:",
        (frecuencia == 1).sum()
    )

    print(
        "Features presentes en ≤5%:",
        (
            frecuencia
            <= numero_organismos * 0.05
        ).sum()
    )

    print(
        "Features presentes en ≥95%:",
        (
            frecuencia
            >= numero_organismos * 0.95
        ).sum()
    )

    print(
        "Features constantes:",
        (
            (frecuencia == 0)
            |
            (frecuencia == numero_organismos)
        ).sum()
    )

    print(
        "\nDistribución:"
    )

    print(
        frecuencia.describe()
    )

N = len(codigos_finales)


resumen_features(
    frecuencia_gene,
    N,
    "GENES"
)


resumen_features(
    frecuencia_KO,
    N,
    "KO"
)


resumen_features(
    frecuencia_module,
    N,
    "MODULES"
)


resumen_features(
    frecuencia_pathway,
    N,
    "PATHWAYS"
)


#Features frecuentes
print("\nGENES:")
print(
    frecuencia_gene.head(20)
)


print("\nKO:")
print(
    frecuencia_KO.head(20)
)


print("\nMODULES:")
print(
    frecuencia_module
)


print("\nPATHWAYS:")
print(
    frecuencia_pathway
)

#Distribución de features por organismo 

resumen_por_organismo = pd.DataFrame({

    "Genes":
        matriz_gene_final.sum(axis=1),

    "KO":
        matriz_KO_final.sum(axis=1),

    "Modules":
        matriz_module_final.sum(axis=1),

    "Pathways":
        matriz_pathway_final.sum(axis=1)

})


print(resumen_por_organismo.describe())

#Guarar las frecuencias

tabla_frecuencia_gene = (
    frecuencia_gene
    .rename("Organismos_presentes")
    .reset_index()
    .rename(
        columns={
            "Gene_ID": "Feature"
        }
    )
)

tabla_frecuencia_KO = (
    frecuencia_KO
    .rename("Organismos_presentes")
    .reset_index()
    .rename(
        columns={
            "KO": "Feature"
        }
    )
)

tabla_frecuencia_module = (
    frecuencia_module
    .rename("Organismos_presentes")
    .reset_index()
    .rename(
        columns={
            "Module_ID": "Feature"
        }
    )
)

tabla_frecuencia_pathway = (
    frecuencia_pathway
    .rename("Organismos_presentes")
    .reset_index()
    .rename(
        columns={
            "Pathway_ID": "Feature"
        }
    )
)


tabla_frecuencia_gene.to_csv(
    Results_DIR /
    "frecuencia_features_gene.csv",
    index=False,
    encoding="utf-8-sig"
)

tabla_frecuencia_KO.to_csv(
    Results_DIR /
    "frecuencia_features_KO.csv",
    index=False,
    encoding="utf-8-sig"
)

tabla_frecuencia_module.to_csv(
    Results_DIR /
    "frecuencia_features_module.csv",
    index=False,
    encoding="utf-8-sig"
)

tabla_frecuencia_pathway.to_csv(
    Results_DIR /
    "frecuencia_features_pathway.csv",
    index=False,
    encoding="utf-8-sig"
)

resumen_por_organismo.to_csv(
    Results_DIR /
    "resumen_features_por_organismo.csv",
    encoding="utf-8-sig"
)

#%%

#===========================================
#Crear matrix definitiva de organismo x KO
#===========================================


matriz_KO_completa = matriz_KO_final.copy()

print("\nMatriz inicial:")
print(matriz_KO_completa.shape)


#Análisis de frecuencia de cada KO

frecuencia_KO = (
    matriz_KO_completa
    .sum(axis=0)
    .sort_values(ascending=False)
)

print(frecuencia_KO)


#Identificar los KO sin presencia 

KO_sin_presencia = (
    frecuencia_KO[
        frecuencia_KO == 0
    ]
    .index
    .tolist()
)


print("\n")
print("KO sin presencia en ningún organismo:")
print(KO_sin_presencia)

print(
    "\nNúmero de KO sin presencia:",
    len(KO_sin_presencia)
)


#Identificar a los KO más informativos 

KO_informativos = (
    frecuencia_KO[
        frecuencia_KO > 0
    ]
    .index
    .tolist()
)


print(
    "\nNúmero de KO informativos:",
    len(KO_informativos)
)


#Contruir la matriz KO definitiva

matriz_KO_definitiva = (
    matriz_KO_completa[
        KO_informativos
    ]
    .copy()
)


#Ordenar organismos

matriz_KO_definitiva = (
    matriz_KO_definitiva
    .reindex(
        sorted(codigos_finales)
    )
)

print(matriz_KO_definitiva.shape)


#Valores permitidos

valores = set(
    matriz_KO_definitiva
    .values
    .flatten()
)

print(
    "Valores presentes:",
    valores
)


#KO vacíos

KO_vacios_final = [
    ko
    for ko in matriz_KO_definitiva.columns
    if matriz_KO_definitiva[ko].sum() == 0
]

print(
    "KO sin organismos:",
    len(KO_vacios_final)
)


#Organismos vacíos

organismos_vacios = (
    matriz_KO_definitiva
    .sum(axis=1)
    .eq(0)
    .sum()
)

print(
    "Organismos sin ningún KO:",
    organismos_vacios
)


#Análisis de la frecuencia final

frecuencia_KO_definitiva = (
    matriz_KO_definitiva
    .sum(axis=0)
    .sort_values(ascending=False)
)

print(frecuencia_KO_definitiva)


#Análisis de los KO raros

KO_raros = (
    frecuencia_KO_definitiva[
        frecuencia_KO_definitiva <= 7
    ]
)

print("Número:", len(KO_raros))


#KO frecuentes

print(frecuencia_KO_definitiva.head(20))


#Resumen

KO_por_organismo = (
    matriz_KO_definitiva
    .sum(axis=1)
    .sort_values(ascending=False)
)

print(KO_por_organismo.describe())


print("\nOrganismos con más KO:", KO_por_organismo.head(10))


print("\nOrganismos con menos KO:", KO_por_organismo.tail(10))


#Guardar

archivo_matriz_KO = (
    Results_DIR /
    "MATRIZ_DEFINITIVA_141x107_ORGANISMO_KO.csv"
)


matriz_KO_definitiva.to_csv(
    archivo_matriz_KO,
    encoding="utf-8-sig"
)


#Guardar las frecuencias

tabla_frecuencia_KO_definitiva = (
    frecuencia_KO_definitiva
    .rename(
        "Organismos_presentes"
    )
    .reset_index()
)


tabla_frecuencia_KO_definitiva.columns = [
    "KO",
    "Organismos_presentes"
]


tabla_frecuencia_KO_definitiva[
    "Frecuencia"
] = (
    tabla_frecuencia_KO_definitiva[
        "Organismos_presentes"
    ]
    /
    len(codigos_finales)
)


tabla_frecuencia_KO_definitiva.to_csv(
    Results_DIR /
    "frecuencia_KO_definitiva.csv",
    index=False,
    encoding="utf-8-sig"
)


#Guardar la lista de los KO excluidos

tabla_KO_excluidos = pd.DataFrame({
    "KO": KO_sin_presencia,
    "Motivo": [
        "Sin presencia en los 141 organismos"
        for _ in KO_sin_presencia
    ]
})


tabla_KO_excluidos.to_csv(
    Results_DIR /
    "KO_excluidos_sin_presencia.csv",
    index=False,
    encoding="utf-8-sig"
)


#%%
#==============================
#Definir la variable objetivo
#==============================

#Definir el label 

LABEL_POSITIVO = 1
LABEL_REFERENCIA = 0


descripcion_labels = {
    1: "Potencial degradador",
    0: "Organismo de referencia"
}

print("Label 1:", descripcion_labels[1])

print("Label 0:", descripcion_labels[0])


# ===============================
#Crear matriz con los candidatos 
# ===============================

tabla_labels_141 = (
    tabla_organismos_final[
        [
            "Organism_code",
            "Organism",
            "Genome_ID",
            "KEGG_Organism"
        ]
    ]
    .copy()
)


tabla_labels_141[
    "Label"
] = pd.NA


tabla_labels_141[
    "Evidence"
] = pd.NA


tabla_labels_141[
    "Evidence_source"
] = pd.NA

tabla_labels_141[
    "Evidence_type"
] = pd.NA

tabla_labels_141[
    "Notes"
] = pd.NA

print(tabla_labels_141.head())

print(tabla_labels_141.columns.tolist())



#======================================
#Realcionar según el perfil funcional
#======================================

tabla_labels_141 = (
    tabla_labels_141
    .merge(
        resumen_funcional[
            [
                "Organism_code",
                "Genes_objetivo",
                "KO_presentes",
                "Modules_presentes",
                "Pathways_presentes"
            ]
        ],
        on="Organism_code",
        how="left"
    )
)

#Ordenar

tabla_labels_141 = (
    tabla_labels_141
    .sort_values(
        "Organism_code"
    )
    .reset_index(
        drop=True
    )
)

print(tabla_labels_141.head(20))


#Guardar

archivo_labels = (
    Results_DIR /
    "tabla_etiquetado_potencial_degradador_141.csv"
)


tabla_labels_141.to_csv(
    archivo_labels,
    index=False,
    encoding="utf-8-sig"
)

print(archivo_labels)

#%%

# ============================================================
# Organismos de referencia 
# ============================================================

BASE_KEGG = "https://rest.kegg.jp"

#KO objetivo

KO_objetivo = set(
    matriz_KO_definitiva.columns
)

print("\nKO objetivo:")
print(len(KO_objetivo))


#MODULES y pathways objetivo

modules_objetivo = sorted(
    tabla_ko_total[
        "Module_ID"
    ]
    .dropna()
    .astype(str)
    .unique()
)

pathways_objetivo = sorted(
    tabla_ko_total[
        "Pathway_ID"
    ]
    .dropna()
    .astype(str)
    .unique()
)


print(len(modules_objetivo))
print(modules_objetivo)


print(len(pathways_objetivo))
print(pathways_objetivo)


#Organismos candidatos

codigos_141 = set(
    tabla_organismos_final[
        "Organism_code"
    ]
    .astype(str)
)


print("\nOrganismos candidatos actuales:")
print(len(codigos_141))

#%%
# ============================================================
#Organsimos candidatos de referencia (positivos y negativos)
# ============================================================

organismos_referencia = {
    "eco": "Escherichia coli K-12 MG1655",
    "kpn": "Klebsiella pneumoniae",
    "sfl": "Shigella flexneri",
    "sty": "Salmonella enterica",
    "bsub": "Bacillus subtilis 168",
    "bcer": "Bacillus cereus ATCC 14579",
    "lmon": "Listeria monocytogenes",
    "lact": "Lactococcus lactis",
    "spn": "Streptococcus pneumoniae",
    "sae": "Staphylococcus aureus",
    "vch": "Vibrio cholerae",
    "paer": "Pseudomonas aeruginosa",
    "hin": "Haemophilus influenzae",
    "nmen": "Neisseria meningitidis"
}


#Excluir a los candidatos

organismos_referencia = {

    codigo: nombre

    for codigo, nombre
    in organismos_referencia.items()

    if codigo not in codigos_141
}

for codigo, nombre in (
    organismos_referencia.items()
):

    print(
        codigo,
        "->",
        nombre
    )


#Obtener KO

def obtener_KO_organismo(
    organism_code
):

    url = (
        f"{BASE_KEGG}/link/ko/"
        f"{organism_code}"
    )

    try:

        respuesta = requests.get(
            url,
            timeout=60
        )


        if respuesta.status_code != 200:

            return {
                "status": "ERROR",
                "status_code":
                    respuesta.status_code,
                "kos": set()
            }


        kos = set()


        for linea in (
            respuesta.text
            .splitlines()
        ):

            partes = linea.split(
                "\t"
            )


            if len(partes) < 2:
                continue


            ko = (
                partes[1]
                .strip()
            )


            ko = ko.replace(
                "ko:",
                ""
            )


            if ko.startswith("K"):

                kos.add(
                    ko
                )


        return {
            "status": "OK",
            "status_code": 200,
            "kos": kos
        }
    except Exception as e:

        return {
            "status": "ERROR",
            "status_code": None,
            "kos": set()
        }

#Obtener los pathways

def obtener_pathways_organismo(
    organism_code
):

    url = (
        f"{BASE_KEGG}/link/pathway/"
        f"{organism_code}"
    )

    try:

        respuesta = requests.get(
            url,
            timeout=60
        )


        if respuesta.status_code != 200:

            return set()


        pathways = set()


        for linea in (
            respuesta.text
            .splitlines()
        ):

            partes = linea.split(
                "\t"
            )


            if len(partes) < 2:
                continue


            pathway = (
                partes[1]
                .strip()
            )


            match = re.search(
                r"(\d{5})$",
                pathway
            )


            if match:

                pathways.add(
                    "map" +
                    match.group(1)
                )


        return pathways


    except Exception:

        return set()


#Cobertura MODULE

def calcular_cobertura_modules(
    kos_organismo
):

    resultados = []


    for module in modules_objetivo:

        datos_module = (
            tabla_ko_total[
                tabla_ko_total[
                    "Module_ID"
                ].astype(str)
                == module
            ]
        )


        kos_module = set(
            datos_module[
                "KO"
            ]
            .dropna()
            .astype(str)
        )


        if len(kos_module) == 0:

            continue


        presentes = (
            kos_module
            .intersection(
                kos_organismo
            )
        )


        cobertura = (
            len(presentes)
            /
            len(kos_module)
        )


        resultados.append({

            "Module_ID":
                module,

            "KO_module":
                len(kos_module),

            "KO_presentes":
                len(presentes),

            "Cobertura":
                cobertura

        })


    return resultados


#Procesar

resultados = []

for i, (
    codigo,
    nombre
) in enumerate(
    organismos_referencia.items(),
    start=1
):

    print(
        f"\n[{i}/"
        f"{len(organismos_referencia)}]"
        f" {codigo} | {nombre}"
    )


#KO

    resultado_KO = (
        obtener_KO_organismo(
            codigo
        )
    )


    kos = resultado_KO[
        "kos"
    ]


    status = resultado_KO[
        "status"
    ]


    print(
        "  Estado:",
        status
    )


    if status != "OK":

        print(
            "  HTTP:",
            resultado_KO[
                "status_code"
            ]
        )


        resultados.append({

            "Organism_code":
                codigo,

            "Organism":
                nombre,

            "KEGG_status":
                status,

            "KO_totales":
                0,

            "KO_objetivo":
                0,

            "KO_objetivo_lista":
                "",

            "Modules_con_KO":
                0,

            "Modules_completos":
                0,

            "Modules_lista":
                "",

            "Pathways_objetivo":
                0,

            "Pathways_lista":
                "",

            "Cobertura_module_promedio":
                0

        })


        time.sleep(
            0.5
        )

        continue


    print(
        "  KO totales:",
        len(kos)
    )


#Intersección con los KO objetivo

    kos_objetivo = (
        kos
        .intersection(
            KO_objetivo
        )
    )


    print(
        "  KO objetivo:",
        len(kos_objetivo)
    )


#Pathways

    pathways = (
        obtener_pathways_organismo(
            codigo
        )
    )


    pathways_presentes = (
        pathways
        .intersection(
            set(
                pathways_objetivo
            )
        )
    )


    print(
        "  Pathways objetivo:",
        len(
            pathways_presentes
        )
    )


#MODULES

    modules_info = (
        calcular_cobertura_modules(
            kos
        )
    )


    modules_con_KO = []

    modules_completos = []


    for registro in (
        modules_info
    ):

        if (
            registro[
                "KO_presentes"
            ] > 0
        ):

            modules_con_KO.append(
                registro[
                    "Module_ID"
                ]
            )


        if (
            registro[
                "Cobertura"
            ] >= 0.80
        ):

            modules_completos.append(
                registro[
                    "Module_ID"
                ]
            )


#Cobertura promedio

    if len(
        modules_info
    ) > 0:

        cobertura_promedio = (
            sum(
                registro[
                    "Cobertura"
                ]

                for registro
                in modules_info
            )
            /
            len(
                modules_info
            )
        )

    else:

        cobertura_promedio = 0


#Guardar

    resultados.append({
        "Organism_code":
            codigo,
        "Organism":
            nombre,
        "KEGG_status":
            status,
        "KO_totales":
            len(kos),
        "KO_objetivo":
            len(kos_objetivo),
        "KO_objetivo_lista":
            ";".join(
                sorted(
                    kos_objetivo
                )
            ),
        "Modules_con_KO":
            len(
                modules_con_KO
            ),
        "Modules_completos":
            len(
                modules_completos
            ),
        "Modules_lista":
            ";".join(
                sorted(
                    modules_con_KO
                )
            ),
        "Pathways_objetivo":
            len(
                pathways_presentes
            ),
        "Pathways_lista":
            ";".join(
                sorted(
                    pathways_presentes
                )
            ),
        "Cobertura_module_promedio":
            cobertura_promedio
    })
    time.sleep(
        0.5
    )


#Dataframe

tabla_referencia = pd.DataFrame(
    resultados
)


#Ordenar

tabla_referencia = (
    tabla_referencia
    .sort_values(
        [
            "KO_objetivo",
            "Modules_completos",
            "Pathways_objetivo"
        ],
        ascending=False
    )
    .reset_index(
        drop=True
    )
)

#%%
# ======================
#Clasificar candidatos
# ======================

def clasificar_referencia(
    fila
):

    if (
        fila[
            "KEGG_status"
        ] != "OK"
    ):

        return "ERROR_KEGG"


    ko = fila[
        "KO_objetivo"
    ]

    modules = fila[
        "Modules_completos"
    ]

    pathways = fila[
        "Pathways_objetivo"
    ]

    if (
        ko >= 5
        and modules <= 1
    ):

        return "A_negativo_dificil"

    if (
        2 <= ko <= 4
    ):

        return "B_intermedio"

    if ko <= 1:

        return "C_distante"

    if (
        ko >= 5
        and (
            modules > 1
            or pathways >= 3
        )
    ):

        return "D_revision_especial"


    return "E_revisar"


tabla_referencia[
    "Grupo_revision"
] = (
    tabla_referencia
    .apply(
        clasificar_referencia,
        axis=1
    )
)

tabla_referencia.loc[
    tabla_referencia[
        "Organism_code"
    ] == "paer",
    "Grupo_revision"
] = "D_revision_especial"


print(
    tabla_referencia[
        [
            "Organism_code",
            "Organism",
            "KO_objetivo",
            "Modules_con_KO",
            "Modules_completos",
            "Pathways_objetivo",
            "Cobertura_module_promedio",
            "Grupo_revision"
        ]
    ]
)

print(
    tabla_referencia[
        "Grupo_revision"
    ].value_counts(dropna=False)
)


negativos_dificiles = (
    tabla_referencia[
        tabla_referencia[
            "Grupo_revision"
        ]
        == "A_negativo_dificil"
    ]
    .copy()
)

print(
    negativos_dificiles[
        [
            "Organism_code",
            "Organism",
            "KO_objetivo",
            "Modules_completos",
            "Pathways_objetivo"
        ]
    ]
)

negativos_intermedios = (
    tabla_referencia[
        tabla_referencia[
            "Grupo_revision"
        ]
        == "B_intermedio"
    ]
    .copy()
)

print(
    negativos_intermedios[
        [
            "Organism_code",
            "Organism",
            "KO_objetivo",
            "Modules_completos",
            "Pathways_objetivo"
        ]
    ]
)

candidatos_distantes = (
    tabla_referencia[
        tabla_referencia[
            "Grupo_revision"
        ]
        == "C_distante"
    ]
    .copy()
)

print(
    candidatos_distantes[
        [
            "Organism_code",
            "Organism",
            "KO_objetivo",
            "Modules_completos",
            "Pathways_objetivo"
        ]
    ]
)


revision = (
    tabla_referencia[
        tabla_referencia[
            "Grupo_revision"
        ]
        == "D_revision_especial"
    ]
    .copy()
)


print(
    revision[
        [
            "Organism_code",
            "Organism",
            "KO_objetivo",
            "Modules_completos",
            "Pathways_objetivo"
        ]
    ]
)

#Guardar
archivo_referencia = (
    Results_DIR /
    "organismos_referencia_candidatos.csv"
)


tabla_referencia.to_csv(
    archivo_referencia,
    index=False,
    encoding="utf-8-sig"
)

negativos_dificiles.to_csv(
    Results_DIR /
    "referencia_A_negativos_dificiles.csv",
    index=False,
    encoding="utf-8-sig"
)

negativos_intermedios.to_csv(
    Results_DIR /
    "referencia_B_intermedios.csv",
    index=False,
    encoding="utf-8-sig"
)

candidatos_distantes.to_csv(
    Results_DIR /
    "referencia_C_distantes.csv",
    index=False,
    encoding="utf-8-sig"
)

revision.to_csv(
    Results_DIR /
    "referencia_D_revision_especial.csv",
    index=False,
    encoding="utf-8-sig"
)

#%%

# ==============================================
#Consolidar organismos candidatos de referencia
# ==============================================

clasificacion_referencia = {
    "eco": {
        "Label": 0,
        "Evidence": "No degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Organismo de referencia no degradador."
    },
    "bcer": {
        "Label": 1,
        "Evidence": "Degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Evidencia de capacidad degradadora."
    },
    "sfl": {
        "Label": 1,
        "Evidence": "Degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Evidencia de capacidad degradadora."
    },
    "bsub": {
        "Label": 1,
        "Evidence": "Degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Evidencia de capacidad degradadora."
    },
    "kpn": {
        "Label": 1,
        "Evidence": "Degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Evidencia de capacidad degradadora."
    },
    "paer": {
        "Label": 1,
        "Evidence": "Degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Evidencia de capacidad degradadora."
    },
    "vch": {
        "Label": 0,
        "Evidence": "No degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Organismo de referencia no degradador."
    },
    "lmon": {
        "Label": 0,
        "Evidence": "No degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Organismo de referencia no degradador."
    },
    "spn": {
        "Label": 0,
        "Evidence": "No degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Organismo de referencia no degradador."
    },
    "lact": {
        "Label": 0,
        "Evidence": "No degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Organismo de referencia no degradador."
    },
    "hin": {
        "Label": 0,
        "Evidence": "No degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Organismo de referencia no degradador."
    },
    "nmen": {
        "Label": 0,
        "Evidence": "No degradador",
        "Evidence_type": "Bibliográfica",
        "Notes": "Organismo de referencia no degradador."
    },
    "sty": {
        "Label": pd.NA,
        "Evidence": "Participación en comunidades degradadoras",
        "Evidence_type": "Indirecta",
        "Notes": (
            "La participación en una comunidad degradadora "
            "no demuestra capacidad degradadora autónoma."
        )
    },
    "sae": {
        "Label": pd.NA,
        "Evidence": "Capacidad descrita para cepas específicas",
        "Evidence_type": "Dependiente de cepa",
        "Notes": (
            "No se generaliza la capacidad degradadora "
            "a la cepa utilizada como referencia."
        )
    }
}
#Información de los organismos
tabla_info = tabla_referencia.copy()

tabla_info = tabla_info[
    [
        "Organism_code",
        "Organism",
        "KEGG_status",
        "KO_totales",
        "KO_objetivo",
        "Modules_con_KO",
        "Modules_completos",
        "Pathways_objetivo"
    ]
].copy()

#Construir tabla de etiquetado
registros = []


for codigo, info in (
    clasificacion_referencia.items()
):

    datos = tabla_info[
        tabla_info[
            "Organism_code"
        ]
        == codigo
    ]


    if len(datos) == 0:

        print(
            "\nADVERTENCIA:",
            codigo,
            "no está en tabla_referencia."
        )

        continue


    fila = datos.iloc[0]


    registros.append({

        "Organism_code":
            codigo,

        "Organism":
            fila["Organism"],

        "KEGG_status":
            fila["KEGG_status"],

        "KO_totales":
            fila["KO_totales"],

        "KO_objetivo":
            fila["KO_objetivo"],

        "Modules_con_KO":
            fila["Modules_con_KO"],

        "Modules_completos":
            fila["Modules_completos"],

        "Pathways_objetivo":
            fila["Pathways_objetivo"],

        "Label":
            info["Label"],

        "Evidence":
            info["Evidence"],

        "Evidence_type":
            info["Evidence_type"],

        "Evidence_source":
            pd.NA,

        "Notes":
            info["Notes"]

    })


tabla_referencia_etiquetada = pd.DataFrame(
    registros
)

tabla_referencia_etiquetada[
    "Label"
] = tabla_referencia_etiquetada[
    "Label"
].astype("Int64")

tabla_referencia_etiquetada = (
    tabla_referencia_etiquetada
    .sort_values(
        [
            "Label",
            "Organism"
        ],
        na_position="last"
    )
    .reset_index(
        drop=True
    )
)

print(
    tabla_referencia_etiquetada[
        [
            "Organism_code",
            "Organism",
            "Label",
            "Evidence",
            "Evidence_type",
            "KO_objetivo",
            "Modules_completos",
            "Pathways_objetivo",
            "Notes"
        ]
    ]
)

print(
    tabla_referencia_etiquetada[
        "Label"
    ]
    .value_counts(
        dropna=False
    )
)

referencia_positivos = (
    tabla_referencia_etiquetada[
        tabla_referencia_etiquetada[
            "Label"
        ] == 1
    ]
    .copy()
)

referencia_negativos = (
    tabla_referencia_etiquetada[
        tabla_referencia_etiquetada[
            "Label"
        ] == 0
    ]
    .copy()
)

referencia_NA = (
    tabla_referencia_etiquetada[
        tabla_referencia_etiquetada[
            "Label"
        ].isna()
    ]
    .copy()
)

print(
    referencia_positivos[
        [
            "Organism_code",
            "Organism",
            "Label",
            "KO_objetivo",
            "Modules_completos",
            "Pathways_objetivo"
        ]
    ]
)

print(
    referencia_negativos[
        [
            "Organism_code",
            "Organism",
            "Label",
            "KO_objetivo",
            "Modules_completos",
            "Pathways_objetivo"
        ]
    ]
)

print(
    referencia_NA[
        [
            "Organism_code",
            "Organism",
            "Label",
            "Evidence",
            "Evidence_type",
            "Notes"
        ]
    ]
)

#Comprobar duplicados

duplicados = (
    tabla_referencia_etiquetada[
        "Organism_code"
    ]
    .duplicated()
    .sum()
)


print(
    "Códigos duplicados:",
    duplicados
)

#Guardar el archivo maestro

archivo_maestro = (
    Results_DIR /
    "tabla_maestra_referencia_etiquetada.csv"
)


tabla_referencia_etiquetada.to_csv(
    archivo_maestro,
    index=False,
    encoding="utf-8-sig"
)

referencia_positivos.to_csv(
    Results_DIR /
    "referencia_positivos.csv",
    index=False,
    encoding="utf-8-sig"
)

referencia_negativos.to_csv(
    Results_DIR /
    "referencia_negativos.csv",
    index=False,
    encoding="utf-8-sig"
)

referencia_NA.to_csv(
    Results_DIR /
    "referencia_excluidos_NA.csv",
    index=False,
    encoding="utf-8-sig"
)

#%%

ARCHIVO_141 = (
    Results_DIR /
    "organismos_finales_141.csv"
)


organismos_141 = pd.read_csv(
    ARCHIVO_141,
    encoding="utf-8-sig"
)


organismos_141.columns = (
    organismos_141.columns
    .str.strip()
)

print(
    "Número de organismos:",
    len(organismos_141)
)


print(
    "Códigos únicos:",
    organismos_141["Organism_code"].nunique()
)


print(
    "\nPrimeros organismos:"
)


print(
    organismos_141[
        [
            "Organism_code",
            "Organism"
        ]
    ].head()
)

#Degradadores documentados e la literatura

degradadores_literatura = [
    {
        "Organism_literatura":
            "Arthrobacter chlorophenolicus A6",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Mycobacterium flavescens PYR-GCK",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Mycobacterium sp. JLS",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Mycobacterium sp. MCS",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Mycobacterium vanbaalenii PYR-1",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Polaromonas naphthalenivorans CJ2",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Rhodococcus erythropolis PR4",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Rhodococcus opacus B4",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Sphingomonas wittichii RW1",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Arthrobacter phenanthrenovorans Sphe3",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },

    {
        "Organism_literatura":
            "Cycloclasticus pugetii PS-1",
        "Hydrocarbon_type":
            "PAH",
        "Evidence":
            "Degradador de PAH"
    },
    {
        "Organism_literatura":
            "Alcanivorax borkumensis SK2",
        "Hydrocarbon_type":
            "Hidrocarburos/petróleo",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Marinobacter hydrocarbonoclasticus VT8",
        "Hydrocarbon_type":
            "Hidrocarburos/petróleo",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Arthrobacter sp. FB24",
        "Hydrocarbon_type":
            "Hidrocarburos/petróleo",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Geobacillus thermodenitrificans NG80-2",
        "Hydrocarbon_type":
            "Hidrocarburos/petróleo",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Alcanivorax sp. DG881",
        "Hydrocarbon_type":
            "Hidrocarburos/petróleo",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Marinobacter hydrocarbonoclasticus ATCC 49840",
        "Hydrocarbon_type":
            "Hidrocarburos/petróleo",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Marinobacter algicola DG893",
        "Hydrocarbon_type":
            "Hidrocarburos/petróleo",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Methylomicrobium album BG8",
        "Hydrocarbon_type":
            "Hidrocarburos",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Oceanicaulis alexandrii HTCC2633",
        "Hydrocarbon_type":
            "Hidrocarburos",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Oleispira antarctica RB-8",
        "Hydrocarbon_type":
            "Hidrocarburos/petróleo",
        "Evidence":
            "Degradador de hidrocarburos"
    },

    {
        "Organism_literatura":
            "Rhodococcus opacus B4 PD630",
        "Hydrocarbon_type":
            "Hidrocarburos",
        "Evidence":
            "Degradador de hidrocarburos"
    }
]


tabla_literatura = pd.DataFrame(
    degradadores_literatura
)


print(tabla_literatura)

#Normalizar nombres
def normalizar_nombre(nombre):

    if pd.isna(nombre):
        return ""

    nombre = str(nombre).lower().strip()

#Eliminar tildes
    nombre = unicodedata.normalize(
        "NFKD",
        nombre
    )

    nombre = "".join(
        c
        for c in nombre
        if not unicodedata.combining(c)
    )

#Eliminar caracteres especiales
    nombre = re.sub(
        r"[^a-z0-9\s]",
        " ",
        nombre
    )

#Espacios múltiples
    nombre = re.sub(
        r"\s+",
        " ",
        nombre
    ).strip()

    return nombre

#%%
# =======================================
#Buscar coincidencia entre las especies
# =======================================

coincidencias = []

for _, candidato in tabla_genomas_148.iterrows():

    nombre_candidato = normalizar_nombre(
        candidato["Organism"]
    )

    for _, literatura in tabla_literatura.iterrows():

        nombre_literatura = normalizar_nombre(
            literatura["Organism_literatura"]
        )

        if nombre_candidato == nombre_literatura:

            coincidencias.append({
                "Organism_code": candidato["Organism_code"],
                "Organism": candidato["Organism"],
                "Organism_literatura": literatura["Organism_literatura"],
                "Hydrocarbon_type": literatura["Hydrocarbon_type"],
                "Evidence": literatura["Evidence"]
            })


coincidencias_especie = pd.DataFrame(
    coincidencias
)


print("Coincidencias:",
    len(coincidencias_especie)
)


#Mostrar coincidencias

if len(coincidencias_especie) > 0:

    print("\nCoincidencias encontradas:")

    print(
        coincidencias_especie[
            [
                "Organism_code",
                "Organism",
                "Organism_literatura",
                "Hydrocarbon_type",
                "Evidence"
            ]
        ]
    )

else:

    print(
        "No se encontraron coincidencias"
    )

#%%
#===========================================
#Organismos de la literatura no encontrados
#===========================================

organismos_encontrados_literatura = set(
    coincidencias_especie[
        "Organism_literatura"
    ]
)


literatura_no_encontrada = (
    tabla_literatura[
        ~tabla_literatura[
            "Organism_literatura"
        ].isin(
            organismos_encontrados_literatura
        )
    ]
    .copy()
)

print(
    "Número:",
    len(literatura_no_encontrada)
)


print(
    literatura_no_encontrada[
        [
            "Organism_literatura",
            "Hydrocarbon_type",
            "Evidence"
        ]
    ]
)


#Candidatos sin coincidencias
codigos_con_evidencia = set(
    coincidencias_especie[
        "Organism_code"
    ]
)


candidatos_sin_evidencia = (
    organismos_141[
        ~organismos_141[
            "Organism_code"
        ].isin(
            codigos_con_evidencia
        )
    ]
    .copy()
)

print(
    "Número:",
    len(candidatos_sin_evidencia)
)


print(
    candidatos_sin_evidencia[
        [
            "Organism_code",
            "Organism"
        ]
    ]
    .head(30)
)


if len(organismos_141) > 0:

    porcentaje = (
        len(codigos_con_evidencia)
        /
        len(organismos_141)
        *
        100
    )

    print(
        "Porcentaje de candidatos con "
        "evidencia bibliográfica:",
        round(
            porcentaje,
            2
        ),
        "%"
    )


#Validación 


tabla_validacion_literatura = (
    coincidencias_especie[
        [
            "Organism_code",
            "Organism",
            "Organism_literatura",
            "Hydrocarbon_type",
            "Evidence"
        ]
    ]
    .copy()
)


tabla_validacion_literatura[
    "Literature_supported"
] = 1


#Guardar

ARCHIVO_VALIDACION = (
    Results_DIR /
    "validacion_141_degradadores_literatura.csv"
)


tabla_validacion_literatura.to_csv(
    ARCHIVO_VALIDACION,
    index=False,
    encoding="utf-8-sig"
)


ARCHIVO_NO_ENCONTRADOS = (
    Results_DIR /
    "degradadores_literatura_no_encontrados.csv"
)


literatura_no_encontrada.to_csv(
    ARCHIVO_NO_ENCONTRADOS,
    index=False,
    encoding="utf-8-sig"
)


ARCHIVO_SIN_EVIDENCIA = (
    Results_DIR /
    "candidatos_141_sin_evidencia_literatura.csv"
)


candidatos_sin_evidencia.to_csv(
    ARCHIVO_SIN_EVIDENCIA,
    index=False,
    encoding="utf-8-sig"
)

#%%
#=================================
#Normalizar degradadores externos
#=================================

ARCHIVO_LITERATURA = (
    Results_DIR /
    "degradadores_literatura_no_encontrados.csv"
)

degradadores_externos = pd.read_csv(
    ARCHIVO_LITERATURA,
    encoding="utf-8-sig"
)

degradadores_externos.columns = (
    degradadores_externos.columns
    .str.strip()
)

print(
    "Número de organismos:",
    len(degradadores_externos)
)

print(
    degradadores_externos[
        [
            "Organism_literatura",
            "Hydrocarbon_type",
            "Evidence"
        ]
    ].to_string(index=False)
)

#Búsqueda en KEGG

def buscar_kegg_genome(nombre):
    url_genome= (
        "https://rest.kegg.jp/find/genome/"
        + requests.utils.quote(nombre)
    )
    try:
        respuesta = requests.get(
            url_genome,
            timeout=60
        )
    except Exception as e:
        return {
            "status": "ERROR",
            "http_status": None,
            "results": [],
            "error": str(e)
        }
    if respuesta.status_code != 200:
        return {
            "status": "ERROR",
            "http_status": respuesta.status_code,
            "results": [],
            "error": respuesta.text[:500]
        }
    texto = respuesta.text.strip()
    if texto == "":
        return {
            "status": "NO_ENCONTRADO",
            "http_status": 200,
            "results": [],
            "error": None
        }
    resultados = []
    for linea in texto.splitlines():
        partes = linea.split("\t")
        if len(partes) >= 2:

            resultados.append({
                "KEGG_ID": partes[0],
                "KEGG_Name": partes[1]
            })
    return {
        "status": "OK",
        "http_status": 200,
        "results": resultados,
        "error": None
    }
resultados = []
for i, fila in degradadores_externos.iterrows():
    nombre = fila[
        "Organism_literatura"
    ]
    print(
        f"\n[{i + 1}/{len(degradadores_externos)}] "
        f"{nombre}"
    )
    resultado = buscar_kegg_genome(
        nombre
    )
    if resultado["status"] == "OK":
        print(
            "  Estado: OK"
        )
        print(
            "  Coincidencias:",
            len(
                resultado["results"]
            )
        )
        for candidato in resultado["results"]:
            resultados.append({
                "Organism_literatura":
                    nombre,
                "Hydrocarbon_type":
                    fila["Hydrocarbon_type"],
                "Evidence":
                    fila["Evidence"],
                "KEGG_ID":
                    candidato["KEGG_ID"],
                "KEGG_Name":
                    candidato["KEGG_Name"],
                "Search_status":
                    "OK"
            })
    elif resultado["status"] == "NO_ENCONTRADO":
        print(
            "  Estado: NO ENCONTRADO"
        )
        resultados.append({
            "Organism_literatura":
                nombre,
            "Hydrocarbon_type":
                fila["Hydrocarbon_type"],
            "Evidence":
                fila["Evidence"],
            "KEGG_ID":
                "",
            "KEGG_Name":
                "",
            "Search_status":
                "NO_ENCONTRADO"
        })
    else:
        print(
            "  Estado: ERROR"
        )
        print(
            "  HTTP:",
            resultado["http_status"]
        )

        resultados.append({
            "Organism_literatura":
                nombre,
            "Hydrocarbon_type":
                fila["Hydrocarbon_type"],
            "Evidence":
                fila["Evidence"],
            "KEGG_ID":
                "",
            "KEGG_Name":
                "",
            "Search_status":
                "ERROR_HTTP"
        })
    time.sleep(
        0.4
    )

#Guardar
tabla_kegg_externos = pd.DataFrame(
    resultados
)

print(
    tabla_kegg_externos.to_string(
        index=False
    )
)

#Resumen de la búsqueda

resumen_busqueda = (
    tabla_kegg_externos
    .groupby(
        "Organism_literatura",
        as_index=False
    )
    .agg(
        Numero_resultados=(
            "KEGG_ID",
            lambda x:
                x.replace(
                    "",
                    pd.NA
                ).dropna().nunique()
        ),
        Status=(
            "Search_status",
            lambda x:
                "OK"
                if "OK" in x.values
                else x.iloc[0]
        )
    )
)


print(resumen_busqueda.to_string(
        index=False
    )
)

#Organismos con una coincidencia
organismos_una_coincidencia = (
    resumen_busqueda[
        resumen_busqueda[
            "Numero_resultados"
        ] == 1
    ]
)

print(
    organismos_una_coincidencia.to_string(
        index=False
    )
)

#Organismos con varias coincidencias

organismos_multiples = (
    resumen_busqueda[
        resumen_busqueda[
            "Numero_resultados"
        ] > 1
    ]
)

print(
    organismos_multiples.to_string(
        index=False
    )
)

if len(organismos_multiples) > 0:

    print("\n")
    print(
        "Detalle coincidencias múltiples:"
    )

    nombres_multiples = set(
        organismos_multiples[
            "Organism_literatura"
        ]
    )

    print(
        tabla_kegg_externos[
            tabla_kegg_externos[
                "Organism_literatura"
            ].isin(
                nombres_multiples
            )
        ].to_string(
            index=False
        )
    )

#Organismos no encontrados

organismos_no_encontrados = (
    resumen_busqueda[
        resumen_busqueda[
            "Numero_resultados"
        ] == 0
    ]
)

print(
    organismos_no_encontrados.to_string(
        index=False
    )
)

tabla_revision_externos = (
    tabla_kegg_externos.copy()
)

tabla_revision_externos[
    "Seleccionado"
] = pd.NA

tabla_revision_externos[
    "Genome_ID"
] = pd.NA

tabla_revision_externos[
    "KEGG_Organism"
] = pd.NA

tabla_revision_externos[
    "Notes"
] = pd.NA

#Guardar
ARCHIVO_RESULTADOS = (
    Results_DIR /
    "normalizacion_degradadores_externos_KEGG.csv"
)


tabla_revision_externos.to_csv(
    ARCHIVO_RESULTADOS,
    index=False,
    encoding="utf-8-sig"
)

print(ARCHIVO_RESULTADOS)

#%%
#==========================================
#Recuperar KO de los degradadores externos
#==========================================

externos_confirmados = pd.DataFrame({
    "Organism_code": [
        "ach",
        "mgi",
        "mmc",
        "rer",
        "abo",
        "gtn",
        "mhc",
        "oai"
    ],
    "Organism": [
        "Arthrobacter chlorophenolicus A6",
        "Mycobacterium flavescens PYR-GCK",
        "Mycobacterium sp. MCS",
        "Rhodococcus erythropolis PR4",
        "Alcanivorax borkumensis SK2",
        "Geobacillus thermodenitrificans NG80-2",
        "Marinobacter hydrocarbonoclasticus ATCC 49840",
        "Oleispira antarctica RB-8"
    ],
    "Genome_ID": [
        "T00834",
        "T00498",
        "T00367",
        "T00881",
        "T00380",
        "T00496",
        "T02083",
        "T04028"
    ]

})

archivo_matriz_KO = (
    Results_DIR /
    "MATRIZ_DEFINITIVA_141x107_ORGANISMO_KO.csv"
)

matriz_KO_definitiva = pd.read_csv(
    archivo_matriz_KO,
    encoding="utf-8-sig"
)

#Normalizar los nombres
matriz_KO_definitiva.columns = [
    str(col).strip().upper()
    for col in matriz_KO_definitiva.columns
]

if "ORGANISM_CODE" in matriz_KO_definitiva.columns:
    matriz_KO_definitiva = (
        matriz_KO_definitiva
        .rename(
            columns={
                "ORGANISM_CODE":
                    "Organism_code"
            }
        )
    )

elif "ORGANISM_CODE" not in matriz_KO_definitiva.columns:
    if "UNNAMED: 0" in matriz_KO_definitiva.columns:
        matriz_KO_definitiva = (
            matriz_KO_definitiva
            .rename(
                columns={
                    "UNNAMED: 0":
                        "Organism_code"
                }
            )
        )

#Extraer los KO
KO_objetivo = []

for col in matriz_KO_definitiva.columns:
    match = re.fullmatch(
        r"K\d{5}",
        str(col).strip().upper()
    )
    if match:
        KO_objetivo.append(
            match.group(0)
        )

KO_objetivo = list(
    dict.fromkeys(
        KO_objetivo
    )
)

print(
    "Número de KO:",
    len(KO_objetivo)
)

print(
    "Primeros KO:",
    KO_objetivo[:20]
)
#=================
#Comprobación KO
#=================
print(
    "¿K01821 está en los KO objetivo?",
    "K01821" in KO_objetivo
)

print(
    "¿K00632 está en los KO objetivo?",
    "K00632" in KO_objetivo
)

#KEGG

def obtener_ko_genes_organismo(
    organism_code
):
    url = (
        "https://rest.kegg.jp/link/ko/"
        + organism_code
    )
    try:

        respuesta = requests.get(
            url,
            timeout=120
        )
    except Exception as e:

        print(
            "  Error de conexión:",
            e
        )
        return pd.DataFrame(
            columns=[
                "KO",
                "Gene_ID"
            ]
        )
    if respuesta.status_code != 200:
        print(
            "  Error HTTP:",
            respuesta.status_code
        )
        return pd.DataFrame(
            columns=[
                "KO",
                "Gene_ID"
            ]
        )
    registros = []
    for linea in respuesta.text.splitlines():

        partes = linea.strip().split(
            "\t"
        )
        if len(partes) != 2:
            continue
        gene_id = partes[0].strip()
        ko_id = partes[1].strip()

#Normalizar
        gene_id = re.sub(
            r"^[^:]+:",
            "",
            gene_id
        )
        ko_id = re.sub(
            r"^ko:",
            "",
            ko_id,
            flags=re.IGNORECASE
        )
        ko_id = ko_id.upper()
#Filtrar

        if ko_id in KO_objetivo:
            registros.append({
                "KO":
                    ko_id,
                "Gene_ID":
                    gene_id

            })
    return pd.DataFrame(
        registros
    )

todos_los_registros = []

for i, fila in externos_confirmados.iterrows():
    codigo = fila[
        "Organism_code"
    ]
    organismo = fila[
        "Organism"
    ]
    print("\n")
    print(
        f"[{i + 1}/8] "
        f"{codigo} | {organismo}"
    )
    tabla = obtener_ko_genes_organismo(
        codigo
    )
    numero_KO = (
        tabla["KO"]
        .nunique()
        if len(tabla) > 0
        else 0
    )
    numero_genes = (
        tabla["Gene_ID"]
        .nunique()
        if len(tabla) > 0
        else 0
    )
    print(
        "  KO objetivo encontrados:",
        numero_KO
    )
    print(
        "  Genes asociados:",
        numero_genes
    )
    if len(tabla) > 0:
        for _, registro in tabla.iterrows():
            todos_los_registros.append({
                "Organism_code":
                    codigo,
                "Organism":
                    organismo,
                "Genome_ID":
                    fila["Genome_ID"],
                "KO":
                    registro["KO"],
                "Gene_ID":
                    registro["Gene_ID"],
                "Present":
                    1,
                "Source":
                    "Literature",
                "Label":
                    1
            })
    time.sleep(
        0.5
    )

#Tabla KO x gene

tabla_externos_KO_gene = pd.DataFrame(
    todos_los_registros
)

print("Registros:",
    len(tabla_externos_KO_gene)
)

print("Organismos:",
    (
        tabla_externos_KO_gene[
            "Organism_code"
        ].nunique()
        if len(tabla_externos_KO_gene) > 0
        else 0
    )
)

print(
    "KO:",
    (
        tabla_externos_KO_gene[
            "KO"
        ].nunique()
        if len(tabla_externos_KO_gene) > 0
        else 0
    )
)

if len(tabla_externos_KO_gene) > 0:

    print(
        tabla_externos_KO_gene.head(20)
    )

matriz_externos_KO = pd.DataFrame(
    0,
    index=externos_confirmados[
        "Organism_code"
    ],
    columns=KO_objetivo
)

for _, fila in tabla_externos_KO_gene.iterrows():
    matriz_externos_KO.loc[
        fila["Organism_code"],
        fila["KO"]
    ] = 1

#Ko x organismo
KO_por_externo = (
    matriz_externos_KO
    .sum(axis=1)
    .sort_values(
        ascending=False
    )
)

print(KO_por_externo)

#Frecuencia KO 

frecuencia_KO_externos = (
    matriz_externos_KO
    .sum(axis=0)
    .sort_values(
        ascending=False
    )
)

print(frecuencia_KO_externos.head(20))

#KO ausentes

KO_sin_externos = (
    frecuencia_KO_externos[
        frecuencia_KO_externos == 0
    ]
    .index
    .tolist()
)

print("Número:",
    len(KO_sin_externos)
)

print(KO_sin_externos)

#Guardar
archivo_matriz_externos = (
    Results_DIR / "MATRIZ_8_DEGRADADORES_EXTERNOS_x_107_KO.csv"
)
matriz_externos_KO.to_csv(
    archivo_matriz_externos,
    encoding="utf-8-sig"
)

archivo_tabla_externos = (
    Results_DIR /
    "tabla_degradadores_externos_KO_gene.csv"
)


tabla_externos_KO_gene.to_csv(
    archivo_tabla_externos,
    index=False,
    encoding="utf-8-sig"
)

#%%
#=================
#Dataset para ML
#=================


archivo_matriz_candidatos = (
    Results_DIR /
    "MATRIZ_DEFINITIVA_141x107_ORGANISMO_KO.csv"
)

matriz_candidatos = pd.read_csv(
    archivo_matriz_candidatos,
    encoding="utf-8-sig"
)

# Normalizar nombres de Organism_code

if "Organism_code" not in matriz_candidatos.columns:
    posibles = [
        "ORGANISM_CODE",
        "Organism_Code",
        "organism_code",
        "Unnamed: 0"
    ]
    encontrado = None
    for columna in posibles:
        if columna in matriz_candidatos.columns:
            encontrado = columna
            break


    if encontrado is not None:
        matriz_candidatos = (
            matriz_candidatos
            .rename(
                columns={
                    encontrado:
                    "Organism_code"
                }
            )
        )

#Normalizar los KO

columnas_no_KO = [
    "Organism_code"
]


KO_columns = [
    col
    for col in matriz_candidatos.columns
    if str(col).strip().upper().startswith("K")
]


KO_columns = [
    str(col).strip().upper()
    for col in KO_columns
]

#Eliminar posibles duplicados
KO_columns = list(
    dict.fromkeys(
        KO_columns
    )
)

print(
    "Organismos:",
    len(matriz_candidatos)
)

print(
    "KO:",
    len(KO_columns)
)

#Información candidatos

archivo_organismos_finales = (
    Results_DIR /
    "organismos_finales_141.csv"
)

organismos_141 = pd.read_csv(
    archivo_organismos_finales,
    encoding="utf-8-sig"
)

# Seleccionar columnas básicas
columnas_basicas = [
    "Organism_code",
    "Organism"
]

for columna in [
    "Genome_ID",
    "KEGG_Organism"
]:

    if columna in organismos_141.columns:
        columnas_basicas.append(
            columna
        )

candidatos_info = (
    organismos_141[
        columnas_basicas
    ]
    .drop_duplicates(
        subset="Organism_code"
    )
)

#Crear Dataset

dataset_candidatos = (
    matriz_candidatos
    .merge(
        candidatos_info,
        on="Organism_code",
        how="left"
    )
)

# Label provisional

dataset_candidatos["Label"] = 1

dataset_candidatos[
    "Evidence_type"
] = "Genomic_potential"

dataset_candidatos[
    "Evidence_source"
] = "Candidate_selection"

dataset_candidatos[
    "Dataset_origin"
] = "Candidate_141"

dataset_candidatos[
    "Evidence_strength"
] = "Provisional"

dataset_candidatos[
    "Notes"
] = (
    "Candidato seleccionado por "
    "criterio genómico de potencial degradador."
)


#Cargar degradarores externos

archivo_externos = (
    Results_DIR /
    "MATRIZ_8_DEGRADADORES_EXTERNOS_x_107_KO.csv"
)

matriz_externos = pd.read_csv(
    archivo_externos,
    encoding="utf-8-sig"
)

if "Organism_code" not in matriz_externos.columns:

    primera_columna = (
        matriz_externos.columns[0]
    )

    matriz_externos = (
        matriz_externos
        .rename(
            columns={
                primera_columna:
                "Organism_code"
            }
        )
    )

#Información de los externos

externos_info = pd.DataFrame({

    "Organism_code": [
        "ach",
        "mgi",
        "mmc",
        "rer",
        "abo",
        "gtn",
        "mhc",
        "oai"
    ],
    "Organism": [
        "Arthrobacter chlorophenolicus A6",
        "Mycobacterium flavescens PYR-GCK",
        "Mycobacterium sp. MCS",
        "Rhodococcus erythropolis PR4",
        "Alcanivorax borkumensis SK2",
        "Geobacillus thermodenitrificans NG80-2",
        "Marinobacter hydrocarbonoclasticus ATCC 49840",
        "Oleispira antarctica RB-8"
    ],
    "Evidence_type": [
        "Literature",
        "Literature",
        "Literature",
        "Literature",
        "Literature",
        "Literature",
        "Literature",
        "Literature"
    ],
    "Evidence_source": [
        "Literature_catalog",
        "Literature_catalog",
        "Literature_catalog",
        "Literature_catalog",
        "Literature_catalog",
        "Literature_catalog",
        "Literature_catalog",
        "Literature_catalog"
    ]

})

externos_info = (
    externos_info
    .merge(
        matriz_externos,
        on="Organism_code",
        how="left"
    )
)
externos_info["Label"] = 1
externos_info[
    "Dataset_origin"
] = "External_literature"
externos_info[
    "Evidence_strength"
] = "Independent"
externos_info[
    "Notes"
] = (
    "Organismo con capacidad degradadora "
    "documentada independientemente."
)

#Cargar las referencias
archivo_referencias = (
    Results_DIR /
    "tabla_maestra_referencia_etiquetada.csv"
)


referencias = pd.read_csv(
    archivo_referencias,
    encoding="utf-8-sig"
)

print(
    referencias[
        [
            "Organism_code",
            "Organism",
            "Label"
        ]
    ]
)

#Matricz de referencias

def obtener_matriz_referencia(
    organism_code
):
    url = (
        "https://rest.kegg.jp/link/ko/"
        + organism_code
    )
    try:
        respuesta = requests.get(
            url,
            timeout=120
        )
    except:
        return {}
    if respuesta.status_code != 200:
        return {}
    resultado = {}
    for linea in respuesta.text.splitlines():

        partes = (
            linea.strip()
            .split("\t")
        )
        if len(partes) != 2:

            continue
        gene_id = partes[0]
        ko = (
            partes[1]
            .replace(
                "ko:",
                ""
            )
            .strip()
            .upper()
        )
        if ko in KO_columns:
            resultado[ko] = 1
    return resultado

archivo_matriz_referencias = (
    Results_DIR /
    "MATRIZ_REFERENCIAS_KO.csv"
)


if os.path.exists(
    archivo_matriz_referencias
):

    matriz_referencias = pd.read_csv(
        archivo_matriz_referencias,
        encoding="utf-8-sig"
    )
    if (
        "Organism_code"
        not in
        matriz_referencias.columns
    ):

        matriz_referencias = (
            matriz_referencias
            .rename(
                columns={
                    matriz_referencias.columns[0]:
                    "Organism_code"
                }
            )
        )
else:

    print("\n")
    print(
        "No existe matriz de referencias."
    )

    print(
        "Se construirá desde KEGG."
    )

    registros_ref = []
    for _, fila in referencias.iterrows():

        codigo = fila[
            "Organism_code"
        ]
        print(
            "Consultando referencia:",
            codigo
        )
        resultado = (
            obtener_matriz_referencia(
                codigo
            )
        )
        registro = {
            "Organism_code":
            codigo
        }
        for ko in KO_columns:

            registro[ko] = (
                resultado.get(
                    ko,
                    0
                )
            )
        registros_ref.append(
            registro
        )
        time.sleep(
            0.5
        )

    matriz_referencias = pd.DataFrame(
        registros_ref
    )


    matriz_referencias.to_csv(
        archivo_matriz_referencias,
        index=False,
        encoding="utf-8-sig"
    )

#Referencias final x label
referencias_final = (
    referencias[
        [
            "Organism_code",
            "Organism",
            "Label",
            "Evidence",
            "Evidence_source",
            "Evidence_type",
            "Notes"
        ]
    ]
    .merge(
        matriz_referencias,
        on="Organism_code",
        how="left"
    )
)


referencias_final[
    "Dataset_origin"
] = "Reference"


referencias_final[
    "Evidence_strength"
] = "Reference"

#Normaizar columnas
columnas_dataset = [
    "Organism_code",
    "Organism"
]

for columna in [
    "Genome_ID",
    "KEGG_Organism"
]:

    if columna in dataset_candidatos.columns:

        columnas_dataset.append(
            columna
        )


columnas_dataset += KO_columns


columnas_dataset += [
    "Label",
    "Evidence_type",
    "Evidence_source",
    "Dataset_origin",
    "Evidence_strength",
    "Notes"
]

candidatos_final = (
    dataset_candidatos
)

for columna in [
    "Genome_ID",
    "KEGG_Organism"
]:
    if columna not in externos_info.columns:
        externos_info[
            columna
        ] = pd.NA
for columna in [
    "Genome_ID",
    "KEGG_Organism"
]:
    if columna not in referencias_final.columns:
        referencias_final[
            columna
        ] = pd.NA

dataset_candidatos = (
    dataset_candidatos[
        columnas_dataset
    ]
)

externos_info = (
    externos_info[
        columnas_dataset
    ]
)

referencias_final = (
    referencias_final[
        columnas_dataset
    ]
)

dataset_maestro = pd.concat(
    [
        dataset_candidatos,
        externos_info,
        referencias_final
    ],
    ignore_index=True
)
#Verificar duplicados
duplicados = (
    dataset_maestro[
        dataset_maestro[
            "Organism_code"
        ]
        .duplicated(
            keep=False
        )
    ]
)

print(
    "Registros duplicados:",
    len(duplicados)
)
if len(duplicados) > 0:
    print(
        duplicados[
            [
                "Organism_code",
                "Organism",
                "Dataset_origin",
                "Label"
            ]
        ]
    )

#Resolver los duplicados
prioridad = {
    "Reference": 3,
    "External_literature": 2,
    "Candidate_141": 1
}


dataset_maestro[
    "_prioridad"
] = (
    dataset_maestro[
        "Dataset_origin"
    ]
    .map(prioridad)
    .fillna(0)
)


dataset_maestro = (
    dataset_maestro
    .sort_values(
        "_prioridad",
        ascending=False
    )
    .drop_duplicates(
        subset="Organism_code",
        keep="first"
    )
    .drop(
        columns="_prioridad"
    )
    .reset_index(
        drop=True
    )
)

print(
    "Organismos:",
    len(dataset_maestro)
)

print(
    "KO:",
    len(KO_columns)
)

print(
    "Dimensiones:",
    dataset_maestro.shape
)

print(
    dataset_maestro[
        "Label"
    ].value_counts(
        dropna=False
    )
)

print(
    dataset_maestro[
        "Dataset_origin"
    ].value_counts()
)

print(
    dataset_maestro[
        "Evidence_type"
    ].value_counts(
        dropna=False
    )
)

#Comprobar los valores de KO
matriz_X = dataset_maestro[
    KO_columns
].copy()

# Convertir a numérico
matriz_X = matriz_X.apply(
    pd.to_numeric,
    errors="coerce"
).fillna(0)

print(
    "Dimensiones X:",
    matriz_X.shape
)

print(
    "Valores:",
    sorted(
        matriz_X
        .stack()
        .unique()
        .tolist()
    )
)

y = dataset_maestro[
    "Label"
].copy()

print(
    y.value_counts(
        dropna=False
    )
)

sin_label = dataset_maestro[
    dataset_maestro[
        "Label"
    ].isna()
]

print(
    "Número:",
    len(sin_label)
)

if len(sin_label) > 0:
    print(
        sin_label[
            [
                "Organism_code",
                "Organism",
                "Dataset_origin"
            ]
        ]
    )

#Guardar dataset
archivo_dataset_maestro = (
    Results_DIR /
    "DATASET_MAESTRO_ML_107_KO.csv"
)


dataset_maestro.to_csv(
    archivo_dataset_maestro,
    index=False,
    encoding="utf-8-sig"
)


#Guardar matriz X
archivo_X = (
    Results_DIR /
    "X_ML_107_KO.csv"
)

matriz_X.to_csv(
    archivo_X,
    index=False,
    encoding="utf-8-sig"
)

#Guardar matriz Y

archivo_y = (
    Results_DIR /
    "y_ML_potencial_degradador.csv"
)

pd.DataFrame({
    "Organism_code":
        dataset_maestro[
            "Organism_code"
        ],
    "Label":
        dataset_maestro[
            "Label"
        ],
    "Dataset_origin":
        dataset_maestro[
            "Dataset_origin"
        ],
    "Evidence_type":
        dataset_maestro[
            "Evidence_type"
        ]
}).to_csv(
    archivo_y,
    index=False,
    encoding="utf-8-sig"
)

print(archivo_dataset_maestro)
print(archivo_X)
print(archivo_y)


#%%
#===================================
#Control de calidad del Dataset ML
#===================================

archivo_dataset = (
    Results_DIR /
    "DATASET_MAESTRO_ML_107_KO.csv"
)

dataset_maestro = pd.read_csv(
    archivo_dataset,
    encoding="utf-8-sig"
)

print("Dimensiones:", dataset_maestro.shape)
print("\nColumnas:")
print(dataset_maestro.columns.tolist())

#Identificar features
columnas_metadata = [
    "Organism_code",
    "Organism",
    "Genome_ID",
    "KEGG_Organism",
    "Label",
    "Evidence_type",
    "Evidence_source",
    "Dataset_origin",
    "Evidence_strength",
    "Notes"
]

#Un KO KEGG tiene formato K + 5 dígitos
patron_KO = re.compile(r"^K\d{5}$")

KO_features = [
    col
    for col in dataset_maestro.columns
    if patron_KO.match(str(col))
]

print("Número de KO:", len(KO_features))

print("\nPrimeros KO:")
print(KO_features[:20])

print("\nÚltimos KO:")
print(KO_features[-10:])

#Comprobación adicional
columnas_no_KO = [
    col
    for col in dataset_maestro.columns
    if col not in KO_features
    and col not in columnas_metadata
]

print("\nColumnas que no son KO ni metadata:")
print(columnas_no_KO)

if len(KO_features) != 107:
    raise RuntimeError(
        f"Se esperaban 107 KO, pero se encontraron "
        f"{len(KO_features)}."
    )

#Organismos únicos
print(
    "Registros:",
    len(dataset_maestro)
)

print(
    "Organismos únicos:",
    dataset_maestro["Organism_code"].nunique()
)

duplicados_codigo = (
    dataset_maestro["Organism_code"]
    .duplicated(keep=False)
)

print(
    "Códigos duplicados:",
    duplicados_codigo.sum()
)

if duplicados_codigo.sum() > 0:

    print("\nDuplicados:")
    print(
        dataset_maestro.loc[
            duplicados_codigo,
            [
                "Organism_code",
                "Organism",
                "Dataset_origin",
                "Label"
            ]
        ]
    )


#Duplicados por nombre
duplicados_nombre = (
    dataset_maestro["Organism"]
    .duplicated(keep=False)
)

print(
    "Nombres duplicados:",
    duplicados_nombre.sum()
)

if duplicados_nombre.sum() > 0:

    print(
        dataset_maestro.loc[
            duplicados_nombre,
            [
                "Organism_code",
                "Organism",
                "Dataset_origin",
                "Label"
            ]
        ]
    )

print(
    dataset_maestro["Label"]
    .value_counts(dropna=False)
)

print(
    dataset_maestro[
        "Dataset_origin"
    ].value_counts(dropna=False)
)


tabla_label_origen = pd.crosstab(
    dataset_maestro["Dataset_origin"],
    dataset_maestro["Label"],
    dropna=False
)

print(tabla_label_origen)

#Matrix X
X = dataset_maestro[
    KO_features
].copy()

print("Dimensiones:", X.shape)

print(
    "Valores únicos:",
    sorted(X.stack().dropna().unique().tolist())
)

#KO constantes
frecuencia_KO = X.sum(axis=0)

KO_constantes = [
    ko for ko in KO_features
    if X[ko].nunique() <= 1
]

print("Número:", len(KO_constantes))

if len(KO_constantes) > 0:
    print(KO_constantes)


#Frecuencia de cada KO
tabla_frecuencia = pd.DataFrame({
    "KO": KO_features,
    "Organismos_presentes": [
        frecuencia_KO[ko]
        for ko in KO_features
    ]
})

tabla_frecuencia["Frecuencia"] = (
    tabla_frecuencia["Organismos_presentes"]
    / len(dataset_maestro)
)

print(
    tabla_frecuencia
    .sort_values(
        "Organismos_presentes",
        ascending=False
    )
    .head(20)
)

#KO raros
KO_raros = tabla_frecuencia[
    tabla_frecuencia[
        "Organismos_presentes"
    ] <= 2
].copy()

print("Número:", len(KO_raros))

print(KO_raros)


#KO presentes en el label=1
dataset_label = dataset_maestro[
    dataset_maestro["Label"].notna()
].copy()

X_label = dataset_label[
    KO_features
]

y_label = dataset_label[
    "Label"
]

frecuencia_positivos = (
    X_label[y_label == 1]
    .sum(axis=0)
)

frecuencia_negativos = (
    X_label[y_label == 0]
    .sum(axis=0)
)

KO_exclusivos_positivos = [
    ko
    for ko in KO_features
    if frecuencia_positivos[ko] > 0
    and frecuencia_negativos[ko] == 0
]

print(KO_exclusivos_positivos)

#KO presentes en el label=0

KO_exclusivos_negativos = [
    ko
    for ko in KO_features
    if frecuencia_negativos[ko] > 0
    and frecuencia_positivos[ko] == 0
]

print( "Número:", len(KO_exclusivos_negativos))
print(KO_exclusivos_negativos)


#Número de KO x organismo
dataset_maestro[
    "KO_presentes"
] = X.sum(axis=1)

print(
    dataset_maestro[
        "KO_presentes"
    ].describe()
)

#KO x label
print(
    dataset_maestro[
        dataset_maestro["Label"].notna()
    ]
    .groupby("Label")[
        "KO_presentes"
    ]
    .describe()
)

#Organismos con pocos KO
organismos_pocos_KO = dataset_maestro[
    dataset_maestro["KO_presentes"] <= 5
][
    [
        "Organism_code",
        "Organism",
        "Dataset_origin",
        "Label",
        "KO_presentes"
    ]
].copy()

print("Número:", len(organismos_pocos_KO))
print(organismos_pocos_KO)

#Organismos con múltiples KO
organismos_muchos_KO = dataset_maestro[
    dataset_maestro["KO_presentes"] >= 40
][
    [
        "Organism_code",
        "Organism",
        "Dataset_origin",
        "Label",
        "KO_presentes"
    ]
].copy()

print("Número:", len(organismos_muchos_KO))

print(organismos_muchos_KO)


#Comparación de las frecuencias 
tabla_comparacion_KO = pd.DataFrame({
    "KO": KO_features,
    "Positivos": [
        frecuencia_positivos[ko]
        for ko in KO_features
    ],
    "Negativos": [
        frecuencia_negativos[ko]
        for ko in KO_features
    ]
})

tabla_comparacion_KO[
    "Diferencia"
] = (
    tabla_comparacion_KO["Positivos"]
    -
    tabla_comparacion_KO["Negativos"]
)

tabla_comparacion_KO = (
    tabla_comparacion_KO
    .sort_values(
        "Diferencia",
        ascending=False
    )
)

print(tabla_comparacion_KO.head(20))


#Comprobar si existen NA
organismos_NA = dataset_maestro[
    dataset_maestro["Label"].isna()
][
    [
        "Organism_code",
        "Organism",
        "Dataset_origin",
        "Evidence_type"
    ]
]

print(organismos_NA)

#Guardar
tabla_frecuencia.to_csv(
    Results_DIR /
    "QC_frecuencia_KO_dataset_ML.csv",
    index=False,
    encoding="utf-8-sig"
)

tabla_comparacion_KO.to_csv(
    Results_DIR /
    "QC_comparacion_KO_positivos_negativos.csv",
    index=False,
    encoding="utf-8-sig"
)

organismos_pocos_KO.to_csv(
    Results_DIR /
    "QC_organismos_pocos_KO.csv",
    index=False,
    encoding="utf-8-sig"
)

organismos_muchos_KO.to_csv(
    Results_DIR /
    "QC_organismos_muchos_KO.csv",
    index=False,
    encoding="utf-8-sig"
)