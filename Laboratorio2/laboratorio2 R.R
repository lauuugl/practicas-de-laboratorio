#laboratorio 2
# cargar libreria tidyverse y data set starwars

library (tidyverse)
data(starwars)

#ejemplos usando select
starwars %>% select (-name)
#Seleccionar sólo las columnas que tienen subraya (_)
starwars %>% select(contains("_"))
#Seleccionar sólo las columnas que empiezan con "s"
starwars %>% select(starts_with("s"))

#Crear un data frame con los nombres y planeta de origen (homeworld)
homeworld <- starwars %>% select(name, homeworld)

#filtrar datos 
# digamos que queremos filtrar por especie(solo humanos)
human <- starwars %>% filter(species == "Human")
#CREACION de DATA FRAMES
#Crear un nuevo datframe con todas las especies menos los Droides
#PUEDO GENERAR DATA FRAMES EXCLUYENDO DATOS 
starwars_nodroids <- starwars %>% filter(species != "Droid")
#PARA FILTRAR
#<	Menor que
#>	Mayor que
#==	Igual que
#<=	Menor o igual que
#>=	Mayor o igual que
#!=	Diferente que
#%in%	Pertenece al conjunto
#is.na	Es NA
#!is.na	No es NA

# Tally  cuenta las observaciones en cada grupo 
# group_by agrupar por una o mas variables 
#Usamos group_by y tally
starwars %>% group_by(species) %>% tally()
# me agrupa las variables de especie y me dice cuanto hay de cada una

#puedo añadir otra variable
starwars %>% group_by(species, gender) %>% tally()
#Si lo quieres guardar en el environment recuerda asignarle un nombre. Hago una tabla con los datos.
table_gender <- starwars %>% group_by(species, gender) %>% tally()

#na.rm=T quiere decir que elima los NA (valores No Asignados o sin datos)
# summarise resume cada grupo en una fila y se usa para calcular aprametros estadisticos con funciones basicas como mean, median, sd, min, max y IQR
starwars %>% group_by(species) %>% summarise(mean_height = mean(height, na.rm = T),mean_mass = mean(mass,na.rm = T))

# COMO CALCULARIAS LA DESVIACION ESTANDAR DE ESOS PARAMETROS?
starwars %>% group_by(species) %>% summarise(Standard_height = sd(height,na.rm =T), standard_mass = sd(mass, na.rm= T))

#Hacer un gráfico de la altura vs. la masa de los personajes
ggplot(starwars, aes(height, mass)) + geom_point()
#ggplot (data, mapping = aes (variables a plotear x,y)). si no le pongo el color al geom_point me la pone negra
#Puedes modificar el color 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red")

#Modificando el color y el punto
ggplot(starwars, aes(height, mass)) + geom_point(colour = "purple", pch = 3)
#agregar el pch cambiar el rombo por un +

#Modificando el color y el fondo 
ggplot(starwars, aes(height, mass)) + geom_point(colour = "red") + theme_light()
#theme light el fondo de la grafica

#indentificar cual es el que pesa mas
table_mass <- starwars %>% group_by(species) %>% summarise(mean_mass = mean(mass,na.rm = T))
# Como ya se cual pesa mas hago un data frame en que no se incluya
table_graphmassheight <- starwars %>% group_by(species) %>% filter (species != "Hutt")
# a partir de este data frame hago mi plot
ggplot(table_graphmassheight, aes(height,mass)) + geom_point(colour = "purple", pch= 3) + theme_light()

# Ejercicio 
#Descarga el dataset toy.csv cargalo en R studio usando la función read_csv de la libreria tidyverse. Tienes que poner la dirección donde has guardado el archivo descargado. En el ejemplo, el archivo está en la carpeta "Descargas"
toy <- read_csv("toy.csv") 
#Como estoy en mi directorio practicasR paso el doc al directorio y ahora si uso el comando de arriba 

#Inspecciona el dataset, haz un resumen de la media (mean) de las variables (Peso, Altura,IMC, IAS, CCintura). Agrupando por sexo.
toy %>% group_by(Sex) %>% summarise(mean_weightkg = mean(Weight_Kg, na.rm = T),mean_height = mean(Height_cm,na.rm = T), meanIMC = mean(IMC, na.nr= T), meanIAS = mean(IAS, na.nr= T), meanCcintura= mean(Ccintura, na.nr= T))

#Haz una tabla sólo con los pacientes femeninos ¿Cuántos registros cumplen las condiciones? ¿De estos cuantos tienen Sobrepeso (Overweight)?

tablewomen <- toy %>% filter(Sex != "Men")
tablewomenoverweight <- toy %>% filter(Sex != "Men") %>% filter(IMC_clas == "Overweight")
#Haz un gráfico usando ggplot relacionando el IMC (Indice de masa corporal) con el peso (Weight_Kg) de todos los pacientes.
ggplot(toy, aes(IMC, Weight_Kg)) + geom_point(colour = "red") + theme_light()

#Repítelo filtrando sólo los pacientes categorizados como "Overweight" y "Obesity".
#la | es como un Y esto
tableobesityoverweight <- toy %>% filter(IMC_clas == "Overweight" | IMC_clas == "Obesity")
ggplot(tableobesityoverweight, aes(IMC, Weight_Kg)) + geom_point(colour = "pink") + theme_light()

install.packages("ape")
install.packages("phangorn")
install.packages("phytools")
