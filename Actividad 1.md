# Actividad 1. Representación gráfica de datos y análisis de resultados

### Introducción

Se analiza a continuación un dataset que contiene información de la expresión de 46 genes en 65 pacientes, cada uno con distintos tipos de tratamiento y características tumorales. La expresión génica puede proporcionar información clave sobre cómo responden diferentes pacientes a tratamientos específicos, especialmente en el contexto del cáncer. Este informe analiza la expresión de genes claves en pacientes bajo dos tratamientos (A y B), con el objetivo de identificar patrones de expresión asociados a cada tratamiento.

En primer lugar, instalé y cargué los paquetes apropiados. En esta actividad utilicé tidyverse para asistir en el análisis de los datos, así como ggplot2 y pheatmap para realizar los gráficos pertinentes. Luego importé el dataset bajo el nombre “exp_genes”.

```{r, include=FALSE}
library(tidyverse)
library(ggplot2)
library(pheatmap)
exp_genes <- read_csv("~/Master en Bioinformática/Estadística y R para Ciencias de la Salud/Actividad 1/Dataset expresión genes.csv")
```

### Distribución de la expresión de genes comparados por tipo de tratamiento

En el dataset se observan, entre otros parámetros, los tratamientos que recibieron los pacientes (A y B) y la expresión de 46 genes. Realicé gráficos de cajas para 15 de estos genes, comparando en cada uno la expresión de cada gen en pacientes bajo el tratamiento A y B. Decidí utilizar un bucle "for" para automiatizar la generación de los gráficos. Como puede observarse en el siguiente código, añadí títulos que se actualizan en cada iteración para reflejar el nombre del gen, escribí etiquetas de ejes, modifiqué las etiquetas de los tratamientos y eliminé la leyenda del lado derecho de los gráficos para simplificar su visualización.

```{r boxplot, fig.width = 3, fig.height = 3, message=FALSE, warning=FALSE}
genes_expresados <- c("AQ_ALOX5", "AQ_CD274", "AQ_CHKA", "AQ_CSF2", "AQ_FOXO3", "AQ_IL6", "AQ_LDHA", "AQ_LIF", "AQ_MAPK1", "AQ_NOS2", "AQ_IFNG", "AQ_PDCD1", "AQ_PPARG", "AQ_TGFB1", "AQ_TNF")
for (gen in genes_expresados) {
  boxplot <- ggplot(exp_genes, aes_string(x = "trat", y = gen, fill = "trat")) +
    geom_boxplot() +
    ggtitle(paste("Expresión de", gen)) +
    labs(x = "Tratamiento", y = "Nivel de Expresión Génica") +
    scale_x_discrete(labels = c("tratA" = "A", "tratB" = "B")) +
    theme(legend.position = "none")
  print(boxplot)
}
```

Para la mayoría de los genes estudiados, puede observarse una distribución más amplia de los valores de expresión génica en los pacientes bajo el tratamiento A. En lo que respecta a la mayoría de los pacientes (esto es, aquellos incluidos en la caja, o sea dentro del rango intercuartílico), la expresión de estos genes suele ser mayor o igual en los pacientes que recibieron el tratamiento A. Aunque en algunos casos puntuales, la mediana de la expresión génica es muy similar para ambos tratamientos, para la mayoría de los genes, dicha mediana toma un valor superior para los pacientes bajo el tratamiento A.

Luego, y teniendo también en cuenta las funciones de los genes analizados, puede concluirse que en los pacientes bajo el tratamiento B se esperara una menor expresión de genes relacionados con el desarrollo y la diferenciación celular, la señalización celular, el ciclo celular, la respuesta inmune e inflamación, y la función y desarrollo de la célula inmunitaria, entre otros.

Dado que los genes estudiados son oncogenes o participan en procesos relacionados con la proliferación celular y la inflamación, una menor expresión podría estar asociada con una reducción en los mecanismos que favorecen el crecimiento tumoral. Esto sugiere que los pacientes bajo tratamiento B, con menor expresión de estos genes pro-tumorales, podrían experimentar un pronóstico más favorable en comparación con aquellos bajo tratamiento A. A su vez, la menor variabilidad en el tratamiento B sugiere una respuesta más homogénea en este grupo.

### Parámetros bioquímicos

El dataset también otorga datos sobre varios parámetros bioquímicos. Realicé histogramas para dichos parámetros con el fin de analizar su distribución en los pacientes estudiados. Nuevamente, utilicé un bucle "for" para automatizar la generación de los gráficos. El título y las etiquetas de los ejes se actualizan en cada iteración para reflejar el nombre del parámetro actual. Cada histograma incluye una línea de densidad roja superpuesta, generada con la función geom_density(), para facilitar la visualización de la distribución.

```{r histograma, fig.width = 3, fig.height = 2.5, message=FALSE, warning=FALSE}
bioquimica <- c("glucosa", "leucocitos", "linfocitos", "neutrofilos", "chol", "hdl", "hierro", "igA", "igE", "igG", "igN", "ldl", "pcr", "transferrina", "trigliceridos", "cpk")
for (parametro in bioquimica) {
  hist <- ggplot(exp_genes, aes_string(x = parametro)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", color = "black") +
    ggtitle(paste("Distribución de", parametro)) +
    xlab(parametro) +
    ylab("Densidad") +
    geom_density(adjust = 1.5, lwd = 1.5, linetype = 2, colour = "red")
  print(hist)
}
```

Puede observarse una distribución no normal para la gran mayoría de los parámetros, lo cual probablemente se deba a que la muestra de pacientes es heterogénea. No solo encontramos pacientes con dos tipos diferentes tratamientos, tres tipos de cáncer (colorrectal, pulmón y mama) y tres posibles extensiones tumorales (localizado, metastásico o regional), sino que también presentan diferentes comorbilidades: hipertensión, diabetes, cardiopatía, neuropatía, hepatopatía, hepatopatía, tabaquismo, etc. En resumen, la diversidad en los tipos de tumores, tratamientos, y comorbilidades contribuye a la variabilidad en los datos, lo que podría ser la causa de las distribuciones no normales de la mayoría de los parámetros bioquímicos.

Dentro de los parámetros estudiados, solo la transferrina y el colesterol siguen una distribución normal o simétrica, aunque los linfocitos se acercan bastante, con unos pocos valores anormalmente altos. De todas formas, no se forma una campana de Gauss perfecta, seguramente debido a la heterogeneidad en la población de pacientes, con variabilidad derivada de los diferentes tratamientos, tipos de tumores y comorbilidades. La presencia de valores atípicos contribuye a que no se logre una distribución normal ideal.

### Mapeo del perfil de expresión génica

Para finalizar el análisis, realicé un heatmap que mapea el nivel de expresión de cada gen por paciente.

Para lograr esto, primero seleccioné las columnas con niveles de expresión génica y escalé los datos para que sean comparables entre sí.

```{r, message=FALSE, warning=FALSE}
pax_genes_scaled <- exp_genes %>% 
  select(starts_with("AQ_")) %>% 
  scale()
```

A continuación, generé un heatmap (mapa de calor) de los datos de expresión génica. Añadí un título y modifiqué las etiquetas de la leyenda a la derecha para que reflejen si la expresión del gen es alta o baja, ya que, después de escalar el dataset, los valores originales no se mostraban en el heatmap. Utilicé las funciones cutree_rows y cutree_cols para agrupar genes y pacientes según patrones similares.

```{r Heatmap, fig.width = 8, fig.height = 5.5, message=FALSE, warning=FALSE}
set.seed(1995)
map <- pheatmap(pax_genes_scaled, 
         main = "Perfil de Expresión Génica",
         legend_breaks = c(1.5, 5, 7.9),
         legend_labels = c("Baja", "Alta", "Expresión\n"),
         legend = TRUE,
         cutree_rows = 2, 
         cutree_cols = 3,
         color = colorRampPalette(c("#2952e1", "#6699cc", "#d2466a", "#e12729"))(100))
print(map)
```
Como puede observarse en el gráfico, la mayoría de los pacientes presenta una disminución significativa en la expresión de la mayoría de los genes. Sin embargo, hay un pequeño grupo de seis pacientes (en la parte superior del gráfico) cuyo nivel de expresión es más alto en varios genes, sugiriendo que, en estos casos, el tratamiento no está logrando reducir la expresión génica de manera tan efectiva. Esto podría indicar que la terapia no está funcionando con la misma eficacia en estos pacientes.

En el resto de los pacientes (en la parte inferior del gráfico), puede notarse cómo en algunos casos la disminución de la expresión es más heterogénea que en otras; es decir, en algunos pacientes, la expresión de los genes se reduce de manera uniforme, mientras que en otros aún se mantienen niveles elevados de algunos genes, lo que sugiere una respuesta más variable al tratamiento. Esto podría reflejar diferencias individuales en la respuesta a la terapia, posiblemente relacionadas con características específicas de los pacientes, como el tipo de tumor o la presencia de comorbilidades.

Con lo que respecta a los genes, a la derecha del gráfico se observan cuatro genes con patrones de expresión diferenciados respecto al resto. En particular, los genes AQ_NOX5 y AQ_ADIPOQ presentan una expresión marcadamente baja en todos los pacientes menos uno. Por otro lado, los genes AQ_SLC2A4 y AQ_ARG1 muestran una expresión ligeramente inferior en comparación con los demás genes, aunque no tan baja como los anteriores. Esta diferencia en los niveles de expresión es evidente incluso en el grupo de pacientes que presenta una mayor expresión génica en general (es decir, los pacientes situados en la parte superior del gráfico). Los cuatro genes situados a la derecha del gráfico se destacan por su expresión más uniforme en todos los pacientes, en contraste con los demás genes, cuya expresión muestra una mayor variabilidad entre los pacientes.

### Conclusión

En resumen, este análisis ha permitido identificar patrones relevantes en la expresión génica de los pacientes bajo diferentes tratamientos. Estos hallazgos sobre una muestra tan heterogénea resaltan la importancia de considerar la variabilidad entre pacientes al evaluar la eficacia de los tratamientos y cómo la expresión génica puede servir como una herramienta poderosa para comprender las respuestas individuales a las terapias.
