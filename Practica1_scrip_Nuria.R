#############################################################################
#
# PRACTICA 1                 NURIA LONGARES IBÁÑEZ
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 01 OCTUBRE 23:59
## Se requiere la entrega de este script completado con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #dimensiones de los datos
head (data) #primeras filas
tail(data) #últimas filas

# Hacemos un primer histograma para explorar los datos
hist(data) ###Archivo Histograma_de_data
hist(data, col= "purple", main="GSE5583 - Histogram") #histograma de color morado y nombre GSE5583 - Histogram  ###Archivo Histograma_de_data_color

# Transformamos los datos con un logaritmo
data_log=log2(data) #tranformación de datos a log en base 2
hist(data_log)  ###Archivo Histograma_de_data_log
hist(data_log, col= "purple", main="GSE5583 - Histogram_log")  ###Archivo Histograma_de_data_log_color

# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
# Si transformamos los datos a logarítmos estos forman una campana de Gauss y siguen una distribución normal en vez de tan sesgada, además las gráficas son más bonitas y más entendibles para publicar

# Hacemos un boxplot con los datos transformados
boxplot(data_log, col=c("blue","blue","blue","orange","orange","orange"), main="GSE5583-boxplots", las=2) ###Archivo Boxplots_GSE5583

# ¿Qué significan los parámetros que hemos empleado?
# boxplot para hacer un gráfico boxplot
# col para cambiar el color
# main para poner el nombre
# las para poner los ejes en vertical

# ¿Qué es un boxplot?
# Un boxplot es un diagrama de cajas y bigotes. Este es un método estandarizado para representar gráficamente una serie de datos numéricos a través de sus cuartiles. 
# De esta manera, se muestran a simple vista la mediana y los cuartiles de los datos, ​ y también pueden representarse sus valores atípicos

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación de los valores de expresión
hc = hclust(as.dist(1-cor(data_log))) #creamos variable hc con el clustering guardado
plot (hc, main= "Hierarchical_Clustering") #representamos hc como plot  ###Archivo Hierarchical_clustering

# ¿Es correcta la separación?
# Si, separa los 3 WT de los 3 KO

#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
# WT columnas de la 1 a la 3
# KO columnas de la 4 a la 6

wt <- data[,1:3]
ko <- data[,4:6]
class(wt) #nos dice el tipo de dato: "matrix" "array"
head(wt)
# Los datos que he generado son de tipo matriz

# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply(wt,1,mean) #calcular la media en todos wt
wt.mean
head(wt.mean)
ko.mean = apply(ko,1,mean) #calcular la media en todos ko
ko.mean
head(ko.mean)

# ¿Cuál es la media más alta?
max(wt.mean) #calcular media más alta de wt
max(ko.mean) #calcular media más alta de ko
# La media más alta es la de ko

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean, xlab="WT", ylab="KO", main="GSE5583-Scatter")  ###Archivo Scatter

# Añadir una línea diagonal con abline
abline(0,1,col="red")  #añade diagonal en color rojo en este caso
abline(h=2,col="blue") #añade línea horizontal en color azul en este caso
abline(v=5, col="green") #añade línea vertical en color verde en este caso

# ¿Eres capaz de añadirle un grid?
grid() #añade cuadrícula   ###Archivo Scatter_color

# Calculamos la diferencia entre las medias de las condiciones
diff.mean= wt.mean-ko.mean
diff.mean  

# Hacemos un histograma de las diferencias de medias
hist(diff.mean) #hacemos un histograma ###Archivo Histograma_diff.mean
hist(diff.mean, col= "pink") #hacemos un histograma de color rosa en este caso  ##Archivo Histograma_diff.mean_color

# Calculamos la significancia estadística con un t-test
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# ¿Cuántas valores tiene cada condición?
pvalue= NULL #lista vacía para guardar los pvalues
tstat= NULL #lista vacía para guardar las estadísticas del test
for(i in 1:nrow(data)) { #para cada gen
	x= wt[i,] #gene wt número 1
	y= ko[i,] #gene ko número 1

	#hacemos el test 
	t = t.test(x,y)

	#añadimos el p-value a la lista
	pvalue[i] = t$p.value
	#añadimos las estadísticas a la lista
	tstat[i] = t$statistic
}

head(pvalue) 

# Usamos los datos sin transformar para que sean datos fiables porque no sufren ninguna modificación
# Cada condición tiene 3 valores, tenemos 2 condiciones(wt y ko), por lo tanto tenemos 6 muestras

# Ahora comprobamos que hemos hecho TODOS los cálculos
length(pvalue)

# Hacemos un histograma de los p-values.
hist(pvalue) ###Archivo Histograma_pvalue

# ¿Qué pasa si le ponemos con una transformación de -log10?
hist(-log10(pvalue),col="green") ###Archivo Histograma_log10_pvalue
# Lo que pasa al cambiar a -log10 es que cambia nuestra distribución

# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "Diff.mean_Volcano") ###Archivo Volcanoplot_diff.mean

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico? Si
diff.mean_cutoff = 2 #limite diferencia de medias
pvalue_cutoff = 0.01 #límite de p-value
abline(v = diff.mean_cutoff, col= "blue", lwd=3) #lwd es para el ancho de la línea
#abline(v = -diff.mean_cutoff, col= "red", lwd=3) #no lo hacemos porque se solapa
abline(h = -log10(pvalue_cutoff), col= "green", lwd=3) ###Archivo Volcanoplot_diff.mean_significancia

# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff  #abs es para los valores absolutos
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data[filter_by_pvalue, ])

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data[filter_combined,]
dim(filtered)
head(filtered)
#426 genes cumplen los ccriterios, los mismos que filter_by_pvalue

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "Volcanoplot_log10_pvalue")
points(diff.mean[filter_combined], -log10(pvalue[filter_combined]), col="red") ###Archivo Volcanoplot_log10_pvalue

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
plot(diff.mean, -log10(pvalue), main = "Volcanoplot_log10_pvalue_separados")
points (diff.mean[filter_combined & diff.mean < 0], -log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0], -log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")  ###Archivo Volcanoplot_log10_pvalue_separados
#diff.mean = wt.mean - ko.mean

# Parece que están al revés porque hemos calculado la diferencia de medias con wt - ko entonces como los sobreexpresados son los ko, los valores de ko son 
# mayores que wt entonces nos sale diferencia negativa y por eso están en la izquierda de la gráfica, para que saliesen al revés podríamos multiplicar por -1 todo

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE) ###Archivo Heatmap_ordenado
 
heatmap(filtered)  ###Archivo Heatmap_filtered

#Los parámetros significan:
#filtered son los datos que vamos a representar 
#cexCol es el tamaño del eje x (columnas)
#Colv y Rowv son los dendogramas (columnas y filas)
#labRow es para quitar el nombre

#Para cambiar los colores del heatmap:
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7, col=hcl.colors(50))  ###Archivo Heatmap_ordenado_diferentecolor

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

# Hacemos nuestro heatmap
heatmap.2(filtered,Rowv=rowv, Colv=colv, cexCol=0.7, col = rev(redblue(256)), scale = "row")  ###Archivo GSE5583_DE_Heatmap


# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered,Rowv=rowv, Colv=colv, cexCol=0.7, col = rev(redblue(256)), scale = "row", labRow=FALSE)
dev.off()  ###Archivo GSE5583_DE_Heatmap_pdf
heatmap.2(filtered,Rowv=rowv, Colv=colv, cexCol=0.7, col = redgreen(75), scale = "row", labRow=FALSE) #cambio los colores pero
#no guardo en pdf  ###Archivo GSE5583_DE_Heatmap_redgreen


# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep = "\t", quote = FALSE) ###Archivo GSE5583_DE

