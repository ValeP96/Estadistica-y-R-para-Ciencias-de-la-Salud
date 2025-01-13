# Actividad 2. Análisis descriptivos en R

## Introducción

Se analiza a continuación un dataset que contiene información de la expresión de 46 genes en 65 pacientes, cada uno con distintos tipos de tratamiento y características tumorales. Se pretende determinar la normalidad de la distribución de la expresión de los genes y calcular los descriptivos de la base de datos en función del tratamiento y del tipo de tumor, así como según la edad categorizando por la mediana.

En primer lugar, instalé y cargué los paquetes apropiados. En esta actividad utilicé tidyverse para asistir en el análisis de los datos, así como gtsummary y ggplot2 para realizar tablas y gráficos respectivamente, nortest para el test de Anderson-Darling y car para realizar la prueba de levene. Luego importé el dataset bajo el nombre “Exp_genes”.

```{r, include=FALSE}
library(tidyverse)
library(ggplot2)
library(gtsummary)
library(nortest)
library(car)

Exp_genes_original <- read_csv("~/Master en Bioinformática/Estadística y R para Ciencias de la Salud/Actividad 2/Dataset expresión genes.csv")
```

Antes de comenzar con el análisis, confirmé la ausencia de valores "NA" en el dataset.

```{r, include=FALSE}
colSums(is.na(Exp_genes_original))
```

## Evaluación de la normalidad de la distribución de datos de expresión génica

Al comenzar a analizar los datos correspondientes a la expresión génica, noté valores anormalmente altos en las columnas correspondientes a AQ_ADIPOQ y AQ_NOX5 para el paciente con id 14:

```{r atipicos}
atipicos <- Exp_genes_original %>%
  select(id, AQ_ADIPOQ, AQ_NOX5) %>%
  arrange(desc(AQ_ADIPOQ))

head(atipicos)
```

Luego, el paciente con ID 14 fue eliminado, reduciendo el tamaño de la muestra a 64 pacientes. Esto garantiza que los análisis posteriores no estén influenciados de manera desproporcionada por valores extremos.

```{r Exp_genes}
Exp_genes <- Exp_genes_original %>%
  filter(id != 14)
```

Después de confirmar la calidad de los datos, procedí a analizar la normalidad de la distribución de la expresión de cada gen.

De forma gráfica, puede observarse a través de histogramas que la distribución no es normal en ningún caso, observándose en la mayoría de los genes una asimetría positiva.

```{r histograma, fig.width = 3, fig.height = 2.5, message=FALSE, warning=FALSE}
genes <- names(Exp_genes)[startsWith(names(Exp_genes), "AQ")]

for (gen in genes) {
  hist <- ggplot(Exp_genes, aes_string(x = gen)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.5) +
    ggtitle(paste("Distribución de", gen))+
    xlab(gen) +
    ylab("Frecuencia")
 print(hist)
}
```

Confirmé la falta de normalidad utilizando el test de Anderson Darling, el cual se recomienda para n ≥ 50 y es robusto cuando hay valores atípicos y problemas en las colas de la distribución:

```{r AD}
length(genes) #46

pvalues_adtest <- data.frame(
  Gen = character(46),                     
  Anderson_Darling_pvalue = numeric(46)    
)

for (i in 1:46) {
  pvalues_adtest[i,1] <- genes[i]
  adtest_result <- ad.test(Exp_genes[[genes[i]]])
  pvalues_adtest[i,2] <- sprintf("%.5f", adtest_result$p.value)
}

print(arrange(pvalues_adtest, desc(Anderson_Darling_pvalue)))

```

Como puede observarse en la tabla, el p-valor para todas las expresiones génicas es menor a 0.05, por lo que se concluye que ningún gen mostró distribución normal según el test de Anderson-Darling. Por lo tanto, los siguientes análisis se basarán en medidas no paramétricas, utilizando la mediana y el rango intercuartílico como estadísticos descriptivos. Para futuros análisis, sería interesante aplicar transformaciones como logaritmos o raíces cuadradas para reducir la asimetría y mejorar la aproximación a la normalidad.

## Estadística descriptiva

Se resumen a continuación las características de los pacientes a través de la estadística descriptiva. En esta sección se presentan los resultados descriptivos de los parametros bioquímicos evaluados, los síntomas reportados, las características sociodemográficas y las comorbilidades de los pacientes, lo que nos permitirá entender mejor las características de la población estudiada.

#### Variables bioquímicas

Para realizar el análisis de las variables bioquímicas, en primer lugar, cree un vector que contiene dichos parámetros y determiné su tamaño:

```{r bq}
bioquimica <- c('glucosa', 'leucocitos', 'linfocitos', 'neutrofilos', 'chol', 'hdl', 'hierro', 'igA', 'igE', 'igG', 'igN', 'ldl', 'pcr', 'transferrina', 'trigliceridos', 'cpk')

length(bioquimica)
```

Cree luego una tabla para resumir los datos de interés. En primer lugar, determiné a través del test de Shapiro (ya que el n es pequeño) si la distribución de cada variable es o no paramétrica, ya que lo indicado para valores con distribución paramétrica es informar media y desviación estándar, mientras que para distribuciones no paramétricas utilizamos mediana y rango intercuartílico.

```{r bq2}
estadistica_bioquimica <- data.frame(
  "Parametro_bioquimico" = character(16), 
  "Media_o_mediana" = numeric(16),
  "SD" = numeric(16),
  "Q1" = numeric(16), 
  "Q3" = numeric(16),
  "Media_SD_o_P50_RIQ" = character(16),
  "Shapiro_pvalue" = numeric(16)
)

for (i in 1:16) {
  estadistica_bioquimica[i,1] <- bioquimica[i] #Completa la primer columna con los parametros bioquimicos
  shapiro_result_bq <- shapiro.test(Exp_genes[[bioquimica[i]]]) #Calcula el p-valor del test de Shapiro
  estadistica_bioquimica[i,7] <- sprintf("%.3f", shapiro_result_bq$p.value)#Guarda dicho valor en la última columna
}

#Creacion de una nueva columna que indica si la variable es parametrica en funcion del p valor
estadistica_bioquimica$Parametricos <- ifelse(as.numeric(estadistica_bioquimica$Shapiro_pvalue) > 0.05, "Sí", "No")

print(arrange(estadistica_bioquimica, desc(Shapiro_pvalue)))

```

Luego de analizar los resultados del test de Shapiro, asigné las variables con p\>0.05 al vector bioquimica_parametricos, mientras que las restantes se incluyeron en bioquimica_no_parametricos Finalmente, mediante un bucle y una estructura condicional, calculé y guardé los estadísticos descriptivos correspondientes.

```{r bq3}
bioquimica_parametricos <- c('chol', 'transferrina')
bioquimica_no_parametricos <- setdiff(bioquimica, bioquimica_parametricos)

for (i in 1:16) {
    if (bioquimica[i] %in% bioquimica_parametricos) {
    # Calcular media y desviación estándar para datos parametricos
    estadistica_bioquimica[i, 2] <- mean(Exp_genes[[bioquimica[i]]])
    estadistica_bioquimica[i, 3] <- sd(Exp_genes[[bioquimica[i]]])
    estadistica_bioquimica[i, 6] <- paste0( #La columna contiene media (SD)
      format(estadistica_bioquimica[i, 2], digits = 2), 
      " (", 
      format(estadistica_bioquimica[i, 3], digits = 2), 
      ")"
    )
    
  } else {
    # Calcular mediana y rango intercuartílico para datos no parametricos
    estadistica_bioquimica[i, 2] <- median(Exp_genes[[bioquimica[i]]])
    estadistica_bioquimica[i, 4] <- quantile(Exp_genes[[bioquimica[i]]], 0.25)
    estadistica_bioquimica[i, 5] <- quantile(Exp_genes[[bioquimica[i]]], 0.75)
    estadistica_bioquimica[i, 6] <- paste0( #La columna contiene mediana (Q1-Q3)
      format(estadistica_bioquimica[i, 2], digits = 2), 
      " (", 
      format(estadistica_bioquimica[i, 4], digits = 2), 
      "-", 
      format(estadistica_bioquimica[i, 5], digits = 2), 
      ")"
      )
  }
}

```

Finalmente, la siguiente tabla resume los valores de media y desviación estándar o mediana y rango intercuartílico según corresponda, para los parámetros bioquímicos evaluados en los pacientes del estudio:

```{r bq tabla}
print(select(estadistica_bioquimica, Parametro_bioquimico, Media_SD_o_P50_RIQ))
```

#### Sintomas

Para resumir la proporción de los síntomas presentados en los pacientes, utilicé `gtsummary`:

```{r sintomas}
resumen_sintomas <- Exp_genes %>%
  select(tos, disnea, expect, secrecion, dolor_garg, escalofrios, fiebre, cansancio, cefalea, mareo, 
         nauseas, vomitos, diarrea, dolor_hueso, dolor_abdo, perd_ape, anosmia, disgueusia) %>%
  tbl_summary() 

resumen_sintomas
```

#### Variables sociodemográficas

Ya que la edad es una variable sociodemográfica continua, primero evalué su distribución con el test de Shapiro-Wilk:

```{r edad}
shapiro.test(Exp_genes$edad)
```

El valor p\>0.05 sugiere que la edad de los pacientes sigue una distribución normal. Por lo tanto, calculé su media y desviación estándar. Tanto para esto como para resumir la proporción de sexo utilicé gtsummary:

```{r sd}
resumen_sociodemograficas <- Exp_genes %>%
  select (edad, sexo) %>%
  tbl_summary(statistic = list(all_continuous() ~ "{mean} ({sd})"))

resumen_sociodemograficas
```

#### Comorbilidades

Finalmente, para resumir la proporción de comorbilidades presentadas en los pacientes utilicé gtsummary:

```{r cmb}
resumen_comorbilidades <- Exp_genes %>%
  select(exfumador, hta, dm, alergia, cardiopatia, neumopatia, hepatopatia, colelitiasis, utolitiasis, ITU, renal, neuropatia)  %>%
  tbl_summary() 

resumen_comorbilidades
```

## Distribución de la expresión génica en función del tratamiento y tipo de tumor

El análisis de la expresión génica según el tratamiento y el tipo de tumor permite identificar patrones que podrían contribuir a comprender las diferencias biológicas entre los tipos de cáncer y su respuesta a los tratamientos.

En este estudio, se analizaron 46 genes en 64 pacientes clasificados según el tratamiento recibido (A o B) y el tipo de tumor (colorrectal, de pulmón y de mama). Habiendo estratificado a los pacientes de esta manera, y dado que la distribución de la expresión génica no es paramétrica, utilicé gtsummary para calcular la mediana y el rango intercuartílico de la expresión de cada gen. Para evaluar si la diferencia de expresión es significativa entre grupos, utilicé el test de Kruskal, prueba no paramétrica indicada para comparar la mediana de tres o más grupos independientes, cuando no se cumplen los supuestos de homogeneidad de varianzas.

```{r tabla2}
tabla2 <- Exp_genes %>%
  select(trat, tumor, all_of(genes)) %>%
  tbl_strata(strata=trat, #estratificamos 1ro por tto
              .tbl_fun = ~.x %>% 
                tbl_summary(by=tumor, #luego, estratificar por tumor
                          statistic = all_continuous() ~ "{median} ({p25} - {p75})",
                          digits = all_continuous() ~ function(x) format(x, digits = 2, scientific = TRUE)) %>%
                      add_p(test = list(all_continuous() ~ "kruskal.test"), 
                      pvalue_fun =label_style_pvalue(digits = 1)))

tabla2
```

Para representar estas diferencias de manera gráfica, generé diagramas de caja que reflejan la expresión génica en función del tratamiento y tipo de tumor:

```{r boxplot, fig.width = 4.5, fig.height = 3, message=FALSE, warning=FALSE}
for (gen in genes) {
  boxplot <- Exp_genes %>%
    ggplot(aes(x = trat, y = .data[[gen]], fill = tumor)) +
    geom_boxplot() +
    labs(title = gen, x = "Tratamiento", y = "Expresión Génica") +
    scale_x_discrete(labels = c("tratA" = "A", "tratB" = "B"))+
    theme_classic()
  print(boxplot)
}
```

En la tabla, los p-valores menores a 0.05 indican que se rechaza la hipótesis nula, lo que sugiere que al menos un grupo de pacientes, clasificado según el tratamiento y el tipo de tumor, tiene una mediana de expresión génica significativamente diferente a las otras. Por ejemplo, el gen AQ_IFNG muestra un valor p=0.001 para los pacientes bajo el tratamiento A, indicando que la expresión de este gen difiere significativamente en al menos uno de los tres tipos de tumores (colorrectal, de pulmón y de mama). Por otro lado, el valor p=0.071 para los pacientes bajo el tratamiento B no alcanza significancia estadística, lo que sugiere que no se detectó una diferencia significativa en la expresión génica según el tipo de tumor en este grupo.

En los pacientes bajo el tratamiento A, 27 genes (AQ_ADIPOQ, AQ_CCR5, AQ_FOXO3, AQ_IRS1, AQ_STAT3, AQ_CCL5, AQ_IFNG, AQ_CCL2, AQ_G6PD, AQ_NFE2L2, AQ_CD36, AQ_MAPK1, AQ_PPARG, AQ_FOXP3, AQ_IL1B, AQ_FASN, AQ_GPD2, AQ_JAK3, AQ_LDHA, AQ_IL10, AQ_BMP2, AQ_CPT1A, AQ_LIF, AQ_SLC2A4, AQ_JAK1, AQ_TLR3 y AQ_TLR4) mostraron diferencias significativas en su expresión según el tipo de tumor. En los pacientes tratados con B, 20 genes (AQ_ADIPOQ, AQ_CHKA, AQ_NOX5, AQ_PTGS2, AQ_ALOX5, AQ_PTAFR, AQ_IL1B, AQ_NFE2L2, AQ_NFKB1, AQ_TLR4, AQ_SOD1, AQ_CPT1A, AQ_JAK1, AQ_SREBF1, AQ_CD274, AQ_MAPK1, AQ_FOXO3, AQ_CCL5, AQ_GPX1 y AQ_PDCD1) presentaron diferencias similares.

Esto sugiere que, dependiendo del tratamiento y tipo de cáncer, ciertos oncogenes y genes relacionados con procesos inflamatorios y metabólicos podrían expresarse de manera diferencial.

Sería valioso analizar la expresión génica antes y después del tratamiento para evaluar su efecto diferencial en distintos tipos de tumores, y poder determinar así la efectividad de los tratamientos en diferentes tipos de cáncer.

## Distribución de la expresión génica en función de la edad

Para proceder con en análisis de la expresión génica en función de la edad, primero dividí las edades en dos categorías según la mediana, o sea, según la edad sea menor o mayor igual a la mediana

```{r cat}
percentil_50 <- median(Exp_genes$edad)

Exp_genes$edad_categoria <- ifelse(Exp_genes$edad < percentil_50, "Edad_cat < percentil 50", "Edad_cat ≥ percentil 50")
Exp_genes$edad_categoria <- as.factor(Exp_genes$edad_categoria)

```

Luego utilicé el test de Levene para evaluar la homogeneidad de varianzas para las expresiones génicas entre ambas categorías de edad. Primero cree el dataframe vacío y luego utilicé un bucle para ejecuta el test de Levene para cada gen y guardar el resultado en la columna pvalue del dataframe. Finalmente añadí una nueva columna que indica si el gen cumple con el supuesto de homogeneidad de varianzas:

```{r levene}
resultados_levene <- data.frame(
  Gen = character(46),                     
  pvalue = numeric(46)    
)

for (i in 1:46) {
  resultados_levene[i,1] <- genes[i]
  levene <- leveneTest(Exp_genes[[genes[i]]], Exp_genes$edad_categoria)
  resultados_levene[i,2] <- levene$`Pr(>F)`[1]
}

resultados_levene$Homogeneidad <- ifelse(resultados_levene$pvalue > 0.05, "Sí", "No")

print(arrange(resultados_levene, pvalue)) 

```

Como puede observarse en la tabla, solo el gen AQ_IL6 no cumple con el supuesto de homogeneidad de varianzas. Considerando además que las muestras son independientes y más de 30, debe utilizarse el test de Welch para determinar si hay una diferencia en la expresión de AQ_IL6 según la categoría de edad. Por el otro lado, para el resto de los genes, ya que hay homogeneidad de varianzas, podemos realizar el test de Student.

Ya que al utilizar t.test el método utiliza por defecto var.equal = FALSE, o sea el test de Welch, y add_p en gtsummary no reconoce el valor de x al intentar modificar los argumentos del test, opté por realizar la tabla de forma manual.

En primer lugar cree dos vectores según se cumpla o no el supuesto de homogeneidad de varianzas y se deba utilizar el test de Student o Welch:

```{r var}
welch_var <- c("AQ_IL6")
ttest_var <- setdiff(genes, welch_var)
```

Luego cree dos dataframes que se distinguen según la categoría de edad y contienen la información sobre la expresión génica:

```{r tablas Edad_cat}
Edad_cat_menor <- Exp_genes %>%
  select(edad_categoria, all_of(genes)) %>%
  filter(edad_categoria == "Edad_cat < percentil 50")

Edad_cat_mayor <- Exp_genes %>%
  select(edad_categoria, all_of(genes)) %>%
  filter(edad_categoria == "Edad_cat ≥ percentil 50")

```

A continuación, cree un nuevo dataframe que contenga la mediana y rango intercuartílico de la expresión de cada gen en pacientes menores a la mediana de edad:

```{r tabla3menor}
tabla3menor <- data.frame(
  "Gen" = character(46),   
  "Mediana menor" = numeric(46),
  "Q1 menor" = numeric(46),
  "Q3 menor" = numeric(46),
  "Edad_cat menor" = character(46)
)

for (i in 1:46) {
  tabla3menor[i,1] <- genes[i]
  tabla3menor[i,2] <- median(Edad_cat_menor[[genes[i]]]) #Calculo de mediana
  tabla3menor[i,3] <- quantile(Edad_cat_menor[[genes[i]]], 0.25) #Calculo de Q1
  tabla3menor[i,4] <- quantile(Edad_cat_menor[[genes[i]]], 0.75) #Calculo de Q3
  tabla3menor[i, 5] <- paste0( #La columna contiene mediana (Q1-Q3)
    format(tabla3menor[i, 2], scientific = TRUE, digits = 2),
    " (", 
    format(tabla3menor[i, 3], scientific = TRUE, digits = 2), 
    " - ", 
    format(tabla3menor[i, 4], scientific = TRUE, digits = 2), 
    ")"
  )
}
```

Repetí el proceso para pacientes mayores a la mediana:

```{r tabla3mayor}
tabla3mayor <- data.frame(
  "Gen" = character(46),   
  "Mediana mayor" = numeric(46),
  "Q1 mayor" = numeric(46),
  "Q3 mayor" = numeric(46),
  "Edad_cat mayor" = character(46)
)

for (i in 1:46) {
  tabla3mayor[i,1] <- genes[i]
  tabla3mayor[i,2] <- median(Edad_cat_mayor[[genes[i]]])
  tabla3mayor[i,3] <- quantile(Edad_cat_mayor[[genes[i]]], 0.25)
  tabla3mayor[i,4] <- quantile(Edad_cat_mayor[[genes[i]]], 0.75)
  tabla3mayor[i, 5] <- paste0( 
    format(tabla3mayor[i, 2], scientific = TRUE, digits = 2), 
    " (", 
    format(tabla3mayor[i, 3], scientific = TRUE, digits = 2), 
    " - ", 
    format(tabla3mayor[i, 4], scientific = TRUE, digits = 2), 
    ")"
  )
}
```

A continuación, uní las columnas conteniendo la lista de genes y la mediana y el rango intercuartílico de ambas categorías de edad en una nueva tabla:

```{r tabla3-1}
tabla3 <- tabla3menor %>%
  select("Gen", "Edad_cat.menor") %>%
  left_join(tabla3mayor %>% select("Gen", "Edad_cat.mayor"), by = "Gen")
```

Creé una nueva columna para guardar los p-valores y utilicé un bucle y una estructura condicional para calcular el test de Student o Welch según se cumpla o no el supuesto de homogeneidad de varianzas:

```{r tabla3-2}
tabla3$pvalue <- NA 

for (i in 1:46) {
  gen <- tabla3$Gen[i]  
  if (gen %in% ttest_var) { # Genes que complen con el supuesto de homogeneidad de varianzas
    test_result <- t.test(Exp_genes[[gen]] ~ Exp_genes$edad_categoria, var.equal = TRUE) # test de Student
  } else if (gen %in% welch_var) { #AQ_IL6
    test_result <- t.test(Exp_genes[[gen]] ~ Exp_genes$edad_categoria, var.equal = FALSE) # Welch
  }
  tabla3$pvalue[i] <- sprintf("%.3f", test_result$p.value) #Se guarda el valor en la nueva columna
}

```

Así, obtuve la siguiente tabla resumiendo la distribución de la expresión génica en función de la edad (como categoría):

```{r tabla3-3}
print(arrange(tabla3, pvalue))
```

En el análisis, los resultados obtenidos muestran que solo el gen AQ_IL6 presenta una diferencia estadísticamente significativa en su expresión según la categoría de edad (p\<0.05). Esto sugiere que la expresión de este gen, asociado a procesos inflamatorios y respuesta inmune, es significativamente mayor en pacientes menores a la mediana de edad.

Por otro lado, para el resto de los genes estudiados, no se detectaron diferencias significativas en la expresión según la edad. Esto podría indicar que, en este estudio, la mayoría de los genes analizados no están influenciados por la edad de los pacientes, o bien que el tamaño muestral limita la capacidad de detectar relaciones más sutiles.

Estos resultados pueden visualizarse también con gráficos de caja:

```{r boxplot_3, fig.width = 3, fig.height = 2.5, message=FALSE, warning=FALSE}
for (gen in genes) {
  boxplot_3 <- Exp_genes %>%
    ggplot(aes(x = edad_categoria, y = .data[[gen]], fill = edad_categoria)) +
    geom_boxplot() +
    labs(title = gen, x = "Categoría de Edad", y = "Expresión Génica") +
    theme_classic()+
    scale_x_discrete(labels = c("Edad_cat < percentil 50" = "< percentil 50", "Edad_cat ≥ percentil 50" = "≥ percentil 50")) +
    theme(legend.position = "none")
  print(boxplot_3)
}
```

En el caso del gen AQ_IL6, para los pacientes de la categoría de edad menor a la mediana, se observa una mayor variabilidad y valores más altos de expresión génica.

Sin embargo, de forma genérica no se observa una correlación estadísticamente significativa entre la edad y la expresión génica en los pacientes de este estudio.

## Conclusión

He utilizado de forma exitosa las herramientas de R para:

-   Evaluar la distribución de la expresión génica y determinar que no sigue una distribución normal, lo que llevó a utilizar el rango intercuartílico como medida central para los análisis posteriores

-   Calcular la mediana y el rango intercuartílico de la expresión génica en función del tratamiento y del tipo de tumor, encontrando diferencias significativas en la expresión génica de ciertos genes entre los pacientes estratificados por estas variables.

-   Examinar los descriptivos estadísticos de la expresión génica en función de la edad, categorizada por la mediana, y determinar que la mayoría de los genes no presentan una correlación significativa con la edad, con la excepción del gen AQ_IL6.

A través de este análisis estadístico, he identificado que la expresión de ciertos genes está correlacionada significativamente con el tipo de tratamiento y el tipo de tumor de los pacientes, sugiriendo que estos factores pueden influir en la variabilidad en la expresión génica. Sin embargo, la mayoría de los genes no muestran una relación significativa con la edad de los pacientes, salvo en el caso de AQ_IL6, donde se observó una mayor expresión en los pacientes menores a la mediana de edad.
