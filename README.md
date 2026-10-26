# TFM
Pipeline de Machine Learning para la predicción del potencial degradador de hidrocarburos en microorganismos
## Descripción 
Este repositorio contiene un pipeline bioinformático orientado a la identificación de microorganismos con potencial para la degradación de hidrocarburos. Para lo cual, se empleo análisis de rutas metabólica, anotación funcional y técnicas de Machine Learning (ML). 

El flujo de trabajo integra información de la base de datos Kyoto Encyclopedia of Genes and Genomes (KEGG), secuencias genómicas disponibles en NCBI y la anotación funcional usando KofamScan. A partir de la presencia de genes y grupos ortólogos de KEGG (KO) asociados a la degradación de hidrocarburos, se construyen matrices de características que permiten entrenar y evaluar los modelos de ML. 

El objetivo de este proyecto es desarrollar un procedimiento reproducible y que permita clasificar microorganismos según su asociación con la degradación de microorganismos, para así priorizar candidatos para una posible aplicación biológica. 

El presente proyecto representa una estimación basada en características genómicas y empleando organismos de referencia, Mas no constituye, por si sola, una demostración experimental. 

## Objetivos
- Recopilar y procesar información genómica y funcional de microorganismos relacionados con la degradación de hidrocarburos.
- Identificar genes y rutas metabólicas asociadas con la degradación de hidrocarburos.
- Construir un conjunto de datos que permita relacionar las características genómicas de los microorganismos con su potencial funcional.
- Evaluar diferentes modelos de machine learning para la identificación de microorganismos con potencial funcional relacionado con la degradación de hidrocarburos.
- Desarrollar un pipeline bioinformático reproducible para la clasificación de microorganismos con potencial degradador. 
- Aplicar el pipeline desarrollado a nuevos datos metagenómicos obtenidos a partir de ecosistemas de páramo para identificar y priorizar candidatos con potencial funcional para la degradación de hidrocarburos.

## Flujo de trabajo 

EL pipeline se encuentra organizado en las siguientes etapas: 

### Etapa 1. Definir las rutas metabólicas 

Se identificaron siete rutas metabólicas de KEGG relacionadas con la degradación de hidrocarburos:

| Código KEGG |	Ruta metabólica |
|:---|:---:|
| map00361 |	Degradación de clorociclohexano y clorobenceno |
| map00622 |	Degradación de xileno |
| map00623 |	Degradación de tolueno |
| map00624 |	Degradación de hidrocarburos aromáticos policíclicos |
| map00626 |	Degradación de naftaleno |
| map00633 |	Degradación de nitrotolueno |
| map00642 |	Degradación de etilbenceno |

A partir de estas rutas se definió los módulos y KO objetivo que serán usados en etapas posteriores. 

### Etapa 2: Recuperación y análisis de información funcional

Se recopiló información mediante la base de datos KEFF para identificar organismos asociados a las rutas metabólicas objetivo y así caracterizar su repertorio funcional. 

Se generaron tablas y matrices que representen la presencia/ausencia de genes, KO, módulos y rutas metabólicas. Con esta información se puede comparar los perfiles funcionales de los microrganismos y establecer criterios de selección de candidatos. 

### Etapa 3. Selección de microorganismos candidatos

La selección se realizó de manera progresiva, considerando la información funcional disponible y a disponibilidad de genomas. 

La selección se realizó de la siguiente manera: 

1.	Identificación de más de 1.000 organismos asociados con la búsqueda funcional inicial.
2.	Selección de organismos mediante criterios de diversidad y priorización.
3.	Exclusión de organismos sin genoma disponible para el análisis previsto.
4.	Aplicación del criterio mínimo de ocho genes objetivo, obteniendo un conjunto final de 141 organismos candidatos.

### Etapa 4. Selección de organismos de referencia y obtención de genomas 

Se definió un conjunto de organismos de referencia positivo y negativos a partir de información de artículos científicos sobre degradación de hidrocarburos. Estos organismos servirán para entrenar los modelos de clasificación. 

La obtención de genomas se realizando mediante la búsqueda en NCBI Assembly y se recopilaron sus identificadores y metadatos. El conjunto final utilizados para el entrenamiento consta de 53 organismos: 27 positivos y 26 negativos. 

### Etapa 5: Anotación funcional con KofanScan 

Los genomas de referencia se sometieron a anotación funcional utilizando la herramienta KofanScan, esta herramienta asigna KO utilizando permifes HMM y unbrales específicos para cada KO. 

La anotación se restringió a los KO objetivo (107 KO) definidos previamente. Los resultados se procesaron para construir una matriz binaria de presencia y ausencia: 

- 1: KO identificado mediante un resultado aceptado por KofamScan.
- 0: KO no identificado como presente en la anotación analizada.

### Etapa 6: Construcción de la matriz de entrenamiento 

A partir de la matriz de anotación funcional se eliminaron aquellos KO sin presencia en ninguno de los organismos de entrenamiento. Obteniendo un total de 78 KO informativos para la construcción de la matriz final de ML. 

La matriz contiene:
- 53 organismos.
- 78 variables KO binarias.
- Una variable de clase (Label).
- Identificadores del organismo y del genoma.

La etiqueta (label) distingue entre organismos positivos y negativos según la clasificación de referencia establecida para el entrenamiento. 

### Etapa 7: Entrenamiento y evaluación de modelos de ML 

A partir de información obtenida en bibliografía, se evaluaron dos algoritmos de clasificación supervisada: 

- Random Forest: clasificador basado en un conjunto de árboles de decisión.
-	Support Vector Machine (SVM): clasificador basado en máquinas de vectores de soporte.

El rendimiento se evaluó mediante las siguientes métricas: 
-	Accuracy.
-	Precision.
-	Recall.
-	F1-score.
-	ROC-AUC.

### Etapa 8: Predicción y priorización de candidatos

El modelos Random Forest (RF) se entrenó con los 53 organismos de referencia y se aplicó a la matriz de candidatos, utilizando exclusivamente los KO presentes en la matriz de entrenamiento. 

## Organización del repositorio

La estructura general del proyecto es la siguiente:

```text
TFM/
├── data/
│   ├── Clasificación tipo hidrocarburos.csv
│   ├── Organismos referencia positivos y negativos.xlsx
│   └── lista_53_genomas_KofamKOALA.csv
│
├── scripts/
│   ├── Annotation/
│   │   ├── KO_objetivo_107.txt
│   │   ├── build_matrix_kofamscan.py
│   │   ├── conteo_KO_53.csv
│   │   ├── resultados_53.zip
│   │   └── run_kofam_53_parallel.sh
│   │
│   ├── Assembly/
│   │   └── Flujo de Assembly.py
│   │
│   ├── entrenamiento ML/
│   │   └── Entrenamiento.py
│   │
│   └── preparación genomas/
│       ├── Preparación y descarga genomas.py
│       └── organismos de referencia.py
│
├── training/
│   └── MATRIZ_KOFAMSCAN_53x78.csv
│
├── Pipeline.py
├── LICENSE
└── README.md
```


Descripción de los componentes

| Componente | Función |
|:---|:---:|
| Pipeline.py | Script principal para integrar las etapas de análisis funcional y selección de candidatos |
| scripts/Assembly/	| Búsqueda, recuperación y revisión de información de ensamblajes genómicos |
| scripts/preparación genomas/ | Preparación de los conjuntos de organismos y descarga de genomas |
| scripts/Annotation/	| Configuración de los KO objetivo, ejecución de KofamScan y construcción de matrices de anotación |
| scripts/entrenamiento ML/	| Entrenamiento, evaluación y aplicación de los modelos de clasificación |
|training/ |	Matriz de características utilizada para el entrenamiento |
|data/ |	Archivos de entrada requeridos por las etapas del pipeline |

## Requisitos de software 

El flujo de trabajo requiere Python y bibliotecas utilizadas en los scripts (pandas, NumPy, Biopython, scikit-learn, Matplotlib, etc), además de Linux junto con Anaconda para el uso de herramientas específicas para la anotación funcional. 

## Ejecución del flujo de trabajo 

Orden de ejecución: 

1.	Preparar los archivos de entrada y la clasificación de organismos de referencia.
2.	Ejecutar el análisis de rutas metabólicas y la selección de candidatos.
3.	Obtener y revisar los ensamblajes genómicos de los organismos de referencia.
4.	Descargar o preparar las secuencias genómicas necesarias.
5.	Ejecutar la anotación funcional con KofamScan.
6.	Construir la matriz binaria de KO y seleccionar las variables utilizadas en el entrenamiento.
7.	Entrenar y evaluar Random Forest y SVM.
8.	Entrenar el modelo final seleccionado.
9.	Aplicar el modelo a los candidatos e integrar las predicciones con la información funcional.

## Resultados principales

El flujo de trabajo desarrollado contempla los siguientes productos:
-	Matriz funcional de 141 candidatos con 107 KO objetivo.
-	Matriz de anotación funcional de 53 organismos de referencia.
-	Matriz de entrenamiento con 78 KO informativos.
-	Métricas de validación cruzada de Random Forest y SVM.
-	Matriz de confusión e importancias de variables KO.
-	Predicciones y ranking de priorización de los 141 candidatos.
- Resumen de cobertura funcional mediante KO, módulos y rutas metabólicas.

## Limitaciones

La capacidad de generalización está condicionada por el número y la representatividad taxonómica de los organismos de referencia, la calidad de los ensamblajes y la precisión de la anotación funcional. La validación experimental y la evaluación en conjuntos de datos independientes constituyen pasos posteriores para comprobar las predicciones.

## Autoría
Proyecto desarrollado como parte del Trabajo de Fin de Máster en Bioinformática.
- **Autora:** Melanie Polo
- **Título del proyecto:** Desarrollo de un pipeline de Machine Learning para predecir funciones biológicas de microorganismos mediante análisis genómico y transcriptómico.
- **Institución:** Universidad Internacional de La Rioja (UNIR).

