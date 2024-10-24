library(ape)
library(phangorn)
library(phytools)

# crear matriz de distancia 
#Para empezar convertimos el archivo en un formato compatible con la librería phangorn, indicando el formato: “FASTA” y el tipo de secuencias: “AA” (aminoácidos)
fraxatin <- read.phyDat(file = "fraxatin_aligned.fasta.txt", 
                        format = "FASTA", type = "AA")
fraxatin
#11 sequences with 310 character and 180 different site patterns.
#The states are A R N D C Q E G H I L K M F P S T W Y V 


#Ahora lo que se hace es crear una matriz de distancia para poder crear árboles de distancia o parsimonia a través de ella. 
#En este caso la clase debe ser AAbin (Amino Acid Sequences), por eso transformamos el objeto fraxatin en este tipo de clase. 
#La función dist.aa (método o función de objeto o clase) calcula una matriz de distancias por pares de las secuencias de aminoácidos partir de una objecto de clase AAbin utilizando un modelo de evolución de aminoácidos (por ejemplo Dayhoff).

matrizdist <- as.AAbin(fraxatin)
matrizdist <- dist.aa(matrizdist)
matrizdist
#los arboles de distancia no representa inferencia evolutiva como el uso de UPGMA

#metodo UPGMA
#Con la matriz de distancia perdemos los caracteres a favor de las diferencias en caracteres entre especies. 
#Cuando aparece un valor de 0, significa que no hay diferencia al nivel de los caracteres (aminoácidos) entre las sequencias de dos especies. Ahora creamos un árbol con el método de grupo de pares no ponderados con media aritmética (UPGMA) usando la matriz de distancia que acabamos de calcular.
#en el metodo upgma la longitud de las ramas no significa nada todas son iguales
arbolUPGMA <- upgma(matrizdist)
plot(arbolUPGMA)
#Para personalizar los árboles podemos agregar argumentos a parámetros como cex, para el tamaño de la letra, edge.color, para el grosos de las ramas, etc. También se puede escoger entre diferentes visualizaciones de árbol como filograma, cladograma, radial y demás.
plot(arbolUPGMA, type= "p", cex=0.8, edge.width=2, edge.color="red", font=3)
plot(arbolUPGMA, type= "p", label.offset=0.0005, edge.lty=1, node.pos=2, cex=0.8, edge.width=2, edge.color="black", font=3)
plot(arbolUPGMA, type= "c", cex=0.8, edge.width=2, edge.color="blue", font=3)

#metodo NJ
arbolNJ <- nj(matrizdist)
plot(arbolNJ)
# longitud de las ramas significa algo: Si la longitud de dos ramas es indéntica, significa que las secuancias también son indénticas y en la matriz de distancia la diferencia es de 0.

# GRAFICAR ARBOLES
#Además de plot podemos graficar árboles con el método plotTree del paquete phytools, el cual es compatible con ape y con phangorn.
plotTree(arbolNJ)
# PONERLE MODIFICACIONES
plotTree(arbolNJ, ftype="b", fsize=0.8, offset=1, color="red", lwd=2)

#SIN CAMBIAR TOPOLOGIA.CAMBIAR EL ORDEN EN QUE LOS GRUPOS SON VISUALIZADOS.
#COMO ORDENAR DE FORMA ALFABETICA O CON LOS GRUPOS MAS DERIVADOS HACIA UNO DE LOS LADOS DEL ARBOL. PARA ESCALERIZAR HACIA LA DERECHA.

plotTree(ladderize(arbolNJ))
#GUARDAR EL ARBOL
write.tree(arbolNJ, file = "file_name.nex")
#PARA VERLO 
read.tree(file = "file_name.nex")

#PARA ENRAIZAR EL ARBOL
#Usamos la funcion root del paquete ape
#para poner por raíz las secuencias de fraxatina de Ornitorinco pasamos por argumento el nombre de la secuencia que corresponde a esta secuencia en el parámetro outgroup (Ornitorrinco).

arbolNJraiz <-root(arbolNJ, outgroup = "Ornitorrinco", r = TRUE)
plot(arbolNJraiz)

#para enraizar el arbol de UPGMA
arbolUPGMAraiz <-root(arbolUPGMA, outgroup = "Ornitorrinco", r=TRUE)
plot(arbolUPGMAraiz)

#VISUALIZACION DE DOS ARBOLES A LA VEZ

layout(matrix(c(1,2)), height=c(10,10))
par(mar=c(1,1,1,1))
plot(arbolUPGMAraiz, label.offset=0.0005, main="ARBOL UPGMA", cex=0.4)
plot(arbolNJraiz, label.offset=0.0005, main="ARBOL NJ", cex=0.4)

#ARBOLES PARSIMONIA 
#PARSIMONIA busca disminuir el número de pasos que explican un árbol evolutivo contando el número de cambios de cada uno de los caracteres.
#En los métodos de distancia (UPGMA y NJ) sólo se llega a un único árbol, mientras que en parsimonia se evalúan múltiples árboles.
#El árbol con menor número de cambios es el ideal.
#Hay que tener en cuenta que la parsimonia no usa todos los caracteres (aminoácidos de la secuecia).
#Son eliminados aquellos que son constantes en todos los taxones, o aquellos que son variables pero no informativos. Un caracter es informativo si tiene al menos dos estados de caracter y por lo menos dos de estos estados ocurren con una frecuencia mínima de dos.

#Para estimar árboles de máxima parsimonia existen varias posbilidades, la más sencilla es partir de árboles de distancia. 
#Se utiliza un árbol de inicio obtenido por distancia y se cuenta su número de pasos. Podemos estimar, por ejemplo, el número de pasos del árbol arbolUPGMAraiz.
parsimony(arbolUPGMAraiz, fraxatin)
#El árbol arbolUPGMAraiz tiene 313 pasos
parsimony(arbolUPGMA, fraxatin)
#algo importante es que aunque esté con raíz o no el número de pasos debe ser el mismo.

# generar arboles con maxima parsimonia a partir de arboles de distancias 
#MAXIMA PARSIMONIA UPGMA
mejorUPGMA <- optim.parsimony(arbolUPGMAraiz, fraxatin)
#MAXIMA PARSIMONIA NJ
mejorNJ <- optim.parsimony(arbolNJraiz, fraxatin)

#METODO HEURISTICO 
#Los métodos heurísticos buscan algunos árboles seleccionando un árbol inicial y cortando ramas y poniéndolas en otro lugar del árbol (branch swapping).
#FORMAS DE HACER BUSQUEDA HEURISTICA
#Intercambio de vecino más cercano (Nearest Neighbor Interchange, NNI): Consiste en intercambir ramas adyacentes hacia una rama interna.
#Podado de subárbol y reinjerto (Subtree pruning and regrafting, SPR): Se corta un pedazo de rama y se inserta en otra rama de manera azarosa.
#Bisección de árbol y reconexión (Tree bisection and reconnection, TBR): Se corta el árbol en la mitad y una porción del árbol que queda de coloca en alguna parte que queda de las ramas del otro árbol.

#estrategia alterna para Otra estrategia para hacer el proceso de búsqueda de árbol con mayor parsimonia es con el algoritmo de búsqueda pratchet. 
#El cual tiene los siguientes pasos:
#1. Generar un árbol de inicio con algún nivel de intercambio de ramas
#2. Seleccionar al azar un subconjunto de caracteres (aminoácidos) y darles más peso, esto quiere decir que se usan con mayor frecuencia dichos caracteres para hacer los análisis. 
#3. La cantidad de caracteres que se seleccionan es establecida por el usuario, típicamente entre 5 y 25% de los caracteres. 
#4.Realizar intercambio de ramas (branch swapping). Iterar pasos 1 a 3 entre 50 y 200 veces. 

#probamos usando la funcion de pratchet
fraxatin_parsimonia <- pratchet(fraxatin, all = TRUE)
#Parsimony score of initial tree: 307 
#Iteration: 100. Best parsimony score so far: 307
#Se genera un árbol, se pesan caracteres por defecto y se itera un número determiando de veces. El algoritmo ha encontrado múltiples ocasiones que el árbol más corto es de 307 pasos pero sólo 4 árboles con igual longitud y número de pasos aunque presenten diferente topología.

fraxatin_parsimonia
#4 phylogenetic trees con igual longitud y numero de pasos pero diferente topologia
# ENRAIZAR LOS 4 PARA COMPARAR
fraxatin_parsimoniaR <- root(phy = fraxatin_parsimonia, outgroup = "Ornitorrinco")
plot(fraxatin_parsimoniaR, cex = 0.6)

#para escoger solo un arbol con igual parsimonia hay que conseguir un árbol de consenso.

# escoger consenso metodo estricto: sólo los grupos monofiléticos presentes en todos los árboles son incluidos; si un grupo monofilétco no está presente en todos los árboles colapsa formando una politomía en el árbol consenso. 
# Para hacer un árbol de consenso estricto podemos usar el método ape con parámetro p de 1, que corresponde a un 100% de consenso entre ramas.
estrictode100 <- consensus(fraxatin_parsimoniaR, p = 1)
plot(estrictode100, cex = .6)

#un arbol estricto pero menos cambiamos el parametro p
# un resultado es un arbol con menos poliotomias mas resuelto
estrictode30 <- consensus(fraxatin_parsimoniaR, p = 0.3)
plot(estrictode30, cex = .6)

#BOOTSTRAP
#Para dar soporte a los árboles se puede hacer una serie de seudoréplicas con remplazamiento (bootstrapping) de una matriz.
#Los cambios entre seudoréplicas consisten en el uso diferencial de los caracteres. Se crean árboles en los que un mismo caracter se repite u otro no se usa.
#Con cada réplica se hace un árbol consenso; la veces que un grupo se repita en el conjunto de réplica es el valor de soporte del nodo.
#UN NODO CON UN SOPORTE MAYOR AL 80% ES UN BUEN SOPORTE 
arbolesbootstrap <- bootstrap.phyDat(fraxatin, FUN = pratchet, bs = 10)
#La rutina anterior genera entonces 10 árboles pseudoréplicas.
plot(arbolesbootstrap, cex = .6)

#generamos un arbol consenso; en este caso con un consenso al 60%:(en este vere los nodos que aparezcan en un 60% de las replicas creadas)
estricto60 <- consensus(arbolesbootstrap, p = 0.6)
plot(estricto60, cex = .6)

#MODELOS PROBABILISTICOS 
#Los caracteres en este caso no son sólo caracteres si no que también representan dinámicas evolutivas
# Todos los caracteres son útiles, incluso los que son constantes entre todos los taxones. Además estos cumplen con un modelo evolutivo.
#Árboles de máxima verosimilitud
#Se calcula la verosimilitud (probabilidad de obtener un dato según un modelo) de un árbol de acuerdo a un alineamiento de secuencias usando un modelo de sustitución de aminoácidos.
#1.Se genera un árbol enraizado de inicio de cualquier tipo, puede ser al azar.
#2.Se calcula la verosimilitud de cada sitio (de aminoácidos) usando un árbol.
#3.Se calcula la verosimilitud total del árbol por nodo, se consideran todos los escenario posibles que pudieron dar origen a los estados de caracter observados para dicho sitio. La verosimilitud de un sitio es, entonces, la suma de verosimilitudes de la reconstrución de este sitio (en forma de árbol) dado un modelo de evolución o sustitución.

#Para el ejemplo usamos de nuevo nuestro objeto fraxatin de clase phy y creamos un árbol al azar de 11 ramas (porque tenemos 11 secuencias) con rtree como punto de partida.
arbolazar <- rtree(n = 11, tip.label = names(fraxatin))
plot(arbolazar, cex = .5)

#enraizamos el arbol por las secuencias de Ornitorinco para poderlo visualizar mejor. 
#los «escalerizamos» hacia la derecha y le agregamos escala;
#aquí la longitud de la rama sí es significativa, indica cantidad de cambio en cuanto a sustituciones de aminoácidos.

arbolazarR <- root(phy = arbolazar, outgroup = "Ornitorrinco")
plot(ladderize(arbolazarR), cex = .5); add.scale.bar()

#A partir del arbol de arriba se puede iniciar la búsqueda del mejor árbol por máxima verosimilitud. 
#Lo primero que se hace es calcular la verosimilitud del árbol dadas las secuencias. Con pml (Phylogenetic maximum likelihood), podemos computar tal verosimilitud.

ajustado <- pml(arbolazarR, fraxatin)
ajustado

#INFORMACION QUE SE OBTIENE
#model: Mk (reporta un modelo de substitución general, el cual tal vez no se ajuste bien a los datos)
#loglikelihood: -4626.347 (Nos reporta la verosimilitud del árbol al azar que habíamos creado, que es -4626.347)
#unconstrained loglikelihood: -1479.871 
#Rate matrix: 

#Lo que hay que hacer es encontrar un árbol que optimice la verosimilitud usando un modelo de sustitución; 
#para esto vamos a usar el método optim.pml del paquete phangorn,
# el cual computa la verosimilitud de un árbol filogenético dado un alineamiento múltiple de secuencias y un modelo de evolución de AA. 
#Toma como argumentos un objeto de clase pml, el tipo de modelo que se quiere usar así como el tiempo de rearreglo para los árboles.

ajustadoconDay <- optim.pml(object = ajustado, model = "Dayhoff", rearrangement = "ratchet")

#ver arbol oculto
ajustadoconDay$tree
# la informacion que nos da: 
#Phylogenetic tree with 11 tips and 9 internal nodes.
#Tip labels:
  #Perro, Rhesus, Ornitorrinco, Raton, Vaca, Caballo, ...
#Node labels:
  #1, 1, 0.407, 1, 0.529, 0.469, ...
#Unrooted; includes branch lengths.

#enraizamos el arbol

ajustadoconDayraíz <- root(ajustadoconDay$tree, outgroup = "Ornitorrinco")
plot(ladderize(ajustadoconDayraíz), cex = .5); add.scale.bar()

# en vez de generar el arbol con el modelo de sustitucion de dayhoff 
# usamos otro modelo como Blosum62
ajustadoconBlo <- optim.pml(object = ajustado, model = "Blosum62", rearrangement = "ratchet")

#usar otro modelo como JTT
ajustadoconJTT <- optim.pml(object = ajustado, model = "JTT", rearrangement = "ratchet")

#AIC criterio de informacion de AKAIKE
AIC(ajustadoconDay, ajustadoconBlo, ajustadoconJTT)
#               df      AIC
#ajustadoconDay 19 5188.890
#ajustadoconBlo 19 5197.264
#ajustadoconJTT 19 5061.388

#La primera columna corresponde a los grados libertad. 
#el mejor modelo que se ajusta con los datos (con el AIC más bajo) es JTT modelo de Jones-Taylor-Thornton el que menos informacion pierde al explicar los datos observados.

# le pedimos que nos muestre el mejor arbol usando dos modelos tanto el de dayhoff como el de JTT
mejorarbol <- optim.pml(
  object = ajustadoconDay, 
  model = "JTT", 
  rearrangement = "ratchet")
mejorarbol
#model: JTT 
#loglikelihood: -2511.694 
#unconstrained loglikelihood: -1479.871 
#Rate matrix: JTT 

#Enraizamos el arbol 
mejorarbolR <- root(mejorarbol$tree, outgroup = "Ornitorrinco")
plot(ladderize(mejorarbolR), cex = 0.5); add.scale.bar()

# DE ESTA MANERA OBTENEMOS EL MEJOR ARBOL CON EL MODELO DE MAXIMA VEROSIMILITUD.
