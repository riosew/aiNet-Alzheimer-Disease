import random
import numpy
import math
import blosum
import time
import os

NFE =0
blosum = blosum.BLOSUM(62,default=0)


def leer_archivo_fasta(archivo):
    nombre = None
    secuencia = ""

    with open(archivo, 'r') as file:
        for linea in file:
            linea = linea.strip()  # Eliminar caracteres de espacio en blanco al principio y al final

            if linea.startswith('>'):
                # Si la línea comienza con '>', se considera como la línea de encabezado (nombre)
                if nombre is not None:
                    # Si ya hemos encontrado un nombre previamente, significa que hemos terminado con la secuencia anterior
                    break
                nombre = linea[1:]  # Elimina el '>' del encabezado y almacena el nombre
            else:
                # De lo contrario, la línea se considera parte de la secuencia
                secuencia += linea

    return nombre, secuencia



def eliminar_columnas_con_guiones(matriz):
    """
    Elimina las columnas de la matriz que contienen únicamente guiones.

    :param matriz: Lista de cadenas de texto que representa la matriz.
    :return: Nueva matriz con las columnas de solo guiones eliminadas.
    """
#    imprimir_matriz(matriz)
    columnas_validas = []
    for j in range(len(matriz[0])):
        columna_contiene_caracter = any(fila[j] != '-' for fila in matriz)
        if columna_contiene_caracter:
            columnas_validas.append(j)

    matriz_sin_guiones = []
    for fila in matriz:
        nueva_fila = ''.join(fila[j] for j in columnas_validas)
        matriz_sin_guiones.append(nueva_fila)

    return matriz_sin_guiones

def aplicar_mutacion_incremento(matriz, timer):
    """
    Aplica una mutación a la matriz, introduciendo guiones en posiciones aleatorias.

    :param matriz: Lista de cadenas de texto que representa la matriz original.
    :param porcentaje: Porcentaje de caracteres a mutar en cada fila.
    :return: Tupla que contiene la nueva matriz mutada y la lista de puntuaciones.
    """

    tiempo_inicial = time.time()
    tiempo_final = tiempo_inicial + timer
    while time.time() < tiempo_final:
        matriz_mutada = []
        for fila in matriz:
            nueva_fila = fila
            #cantidad_agregar = int(len(fila) * timer)
            cantidad_agregar = 1
            #cantidad_agregar = porcentaje # numero fijo
            for _ in range(cantidad_agregar):
                indice = random.randint(0, len(nueva_fila) - 1)
                if random.randint(0,100)  < 100:                                         #cada fila puede mutar con prob
                    nueva_fila = nueva_fila[:indice] + '-' + nueva_fila[indice:]
                    
            matriz_mutada.append(nueva_fila)
        matriz = matriz_mutada


    matriz = completar_filas(matriz)
    matriz = eliminar_columnas_con_guiones(matriz)
      
 
#    imprimir_matriz(matriz_mutada)
    puntuacionArray, puntos =  calculaFitness(matriz)
    return matriz, puntuacionArray



def completar_filas(matriz):
    # Encuentra la longitud de la fila más larga
    newmatriz = []
    longmax = max(len(fila) for fila in matriz)
    for fila in matriz:
        newfila = fila
        while len(newfila) < longmax:
            if random.randint(0,1)==1:       #inserta en cominezo o final aleatoriamente
                newfila = newfila + '-'
            else: newfila = '-' + newfila 
        newmatriz.append(newfila)
#    imprimir_matriz(newmatriz)
    return newmatriz


def calculaFitness(matriz_mutada):
    num_filas = len(matriz_mutada)
    num_columnas = len(matriz_mutada[0]) if num_filas > 0 else 0
    puntos_columna = 0
    puntuacion = []
    puntos = 0
    global NFE
    for columna in range(num_columnas):
        columna = [fila[columna] for fila in matriz_mutada]
#         print(f'Columna {columna}: {columna}')
        puntos_columna = calcular_puntuacion_blosum62(columna)
        puntuacion.append(puntos_columna)
        puntos = puntos + puntos_columna
        NFE = NFE + 1
    
    return puntuacion, puntos  


def calcular_puntuacion_blosum62(columna):
    puntos = 0
    for i in range(len(columna) - 1):
        aminoacido1 = columna[i]
        aminoacido2 = columna[i + 1]
#       print("smion3: ",aminoacido1, " : ",aminoacido2 )
        puntos += blosum[aminoacido1][aminoacido2]
        if aminoacido1=="-" or aminoacido2 == "-":
            puntos = puntos -2
            
    
#    print(puntos)
    return puntos




def leer_palabras():
    """
    Lee un conjunto de Caracteres desde la entrada estándar.

    :param num_filas: Número de filas a leer.
    :return: Cadena de Genoma.
    """
    palabras = []
    names = []

    name, seq = leer_archivo_fasta("C:\\ainet\\ADAM9.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\ADAM17.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\APP.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\MAPK3.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\MAPK8.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\MAPK13.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\PKN2.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\PRKACA.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\PRKCI.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\PRKCQ.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\PSEN1.FASTA")
    palabras.append(seq)
    names.append(name)
    name, seq = leer_archivo_fasta("C:\\ainet\\PSEN2.FASTA")
    palabras.append(seq)
    names.append(name)
    

#     palabras.append("MURCIELAGO")
#     palabras.append("CATAMARAN")
#     palabras.append("COMPATRIOTA")


    return names, palabras

def crear_matriz(palabras, max_longitud):
    """
    Crea una matriz a partir de una lista de caracteres, rellenando con guiones según sea necesario.

    :param caracteres: Lista de caracteres.
    :param max_longitud: Longitud máxima de una palabra en la lista, para estandarizar la longitud de las filas.
    :return: Matriz creada a partir de las palabras proporcionadas.
    """
    matriz = []
    for palabra in palabras:
        matriz.append(palabra.ljust(max_longitud, '-'))
    return matriz

def imprimir_matriz(matrix,names, opt):
    matriz, puntos = matrix    
    if opt=="matriz+fitness" or opt=="matriz":
        for fila in matriz:
            print(fila)
    
    if opt=="fitness" or opt=="matriz+fitness":
        puntuacion, puntos = calculaFitness(matriz)
        print("Best fitness last Gen: ",puntos)
        
    if opt=="file":
        with open("aiNet_out.FASTA", "w") as archivo:
            for name, fila in zip(names, matriz):
                archivo.write(">{}\n{}\n".format(name, fila))

#---------------------------------------------------------------
        
  
        
def mutation_rate(beta, normalized_cost):
  return (1.0/beta) * numpy.exp(-normalized_cost)        
        
def calculate_normalized_cost(fitness, rango):
#   print ("rango ", rango)
    normcost = 1.0-(fitness/rango)
#    print("normCost ", normcost)
    if abs(normcost) >= 709.78:
        print("------------------------------------------salvado")
        normcost = 709.77
    return normcost
  

def gradoMutacion(linfosito, rango):
    arrayPuntos, fitness = calculaFitness(linfosito)
    alpha = mutation_rate(100, calculate_normalized_cost(fitness, rango))  #beta = 100
    grado = alpha * numpy.random.randn()  #random gaussian
    #grado =   (len(linfosito)*len(linfosito[0])) * grado / 1000

#    print(grado)
    return abs(grado)   

#---------------------------------------------------------------


def primer_digito_no_cero(numero):
    numero_str = str(numero)  # Convierte el número a una cadena
    for caracter in numero_str:
        if caracter.isdigit() and caracter != '0':
            return int(caracter)
    return None  # Devuelve None si no se encuentra un dígito no cero



def calculaRango(poblacion):
    seqs, puntosI = poblacion[0]   #mejor y peor valor
    seqs, puntosF = poblacion[-1]
    
    worst = sum(map(lambda x: x, puntosI))
    best  = sum(map(lambda x: x, puntosF))
    rango = abs(best - worst)
    if rango == 0:
        rango=1
    return rango


def insertaRandomClones(matriz_original,poblacion,percent):
    N = int(len(poblacion) * percent)
    poblacion = poblacion + [aplicar_mutacion_incremento(matriz_original.copy(), 0.001) for _ in range(N)]
    return poblacion

def extincion(poblacion_ordenada):
    for linfo in poblacion_ordenada[1:]:
        if random.randint(0, 100) < 90:
            poblacion_ordenada.remove(linfo)
    return poblacion_ordenada
    
def reduccion(poblacion_ordenada, val):
    while len(poblacion_ordenada) > val:
        for linfo in poblacion_ordenada[50:]:
            if random.randint(0, 100) < 50:
                poblacion_ordenada.remove(linfo)
    return poblacion_ordenada    
    
    
def generaReporte(nCorrida, generacion,  popSize, puntos, fitnessProm, NFE, tfinal):
    archivo_existente = os.path.isfile("reporteAiNet.csv")

    with open("reporteAiNet.csv", "a") as archivo:
        if not archivo_existente:
            # Si el archivo no existe, escribe la cabecera
            archivo.write("nCorrida,Generacion,popSize,bestFitness,Fitness_Prom,NFE,time\n")

        archivo.write("{},{},{},{},{},{},{}\n".format(
            nCorrida,
            generacion,
            popSize,
            puntos,
            fitnessProm,
            NFE,
            tfinal
        ))
  

def main():
    
    for nCorrida in range(1, 31):   
        #--------------------------------------------------------------------------
            pobSize = 5
            clonRate = 0.2
            threshold = -16000
            cycles = 30
            nClones = 10
            rClonInsert = 0.3
            maxPob = 80
        #--------------------------------------------------------------------------
            names, palabras = leer_palabras()
            
            # Encontrar la longitud de la palabra más larga para estandarizar la matriz
            max_longitud = max(len(palabra) for palabra in palabras)
            tinicial = time.time()

            matriz_original = crear_matriz(palabras, max_longitud)

            N = pobSize
            poblacion = [aplicar_mutacion_incremento(matriz_original.copy(), clonRate) for _ in range(N)]


            # Definir el número de generaciones
            num_generaciones = cycles
            for generacion in range(1, num_generaciones+1):
                print(f"\nGeneración {generacion }\n{'='*15}")

                # Ordenar la población por puntuación
                poblacion_ordenada = sorted(poblacion, key=lambda x: sum(x[1]), reverse=True)
                print ("poblacion size: ", len(poblacion_ordenada))

             
           

                rango = calculaRango(poblacion_ordenada)

                clones = []
                

                for linfo in poblacion_ordenada:                        #elimina por valor de fitness tresshold
                    fitness = sum(linfo[1])
                    if fitness < threshold:
                        poblacion_ordenada.remove(linfo)
          

        #         if random.randint(0,100) < 25:
        #             print("------------------------------extincion")
        #             poblacion_ordenada = extincion(poblacion_ordenada)
        
        
        
                if len(poblacion_ordenada) > maxPob:
                    poblacion_ordenada = reduccion(poblacion_ordenada, maxPob)
                    print("nueva pov", len(poblacion_ordenada))        
                        
                        

                c=0
                fitnessProm=0
                for linfo in poblacion_ordenada:                              #calcula promedio
                    c=c+1
                    fitnessProm = fitnessProm + sum(linfo[1])
                    
                fitnessProm = int(fitnessProm/c)
                print("fitness Promedio: ",fitnessProm)
        #         
        #         for linfosito in poblacion_ordenada:                        #elimina por promedio
        #             fitness = sum(linfo[1])
        #             if fitness < fitnessProm:
        #                 poblacion_ordenada.remove(linfosito)

                    
                for _ in range(nClones):     
                    for clon in poblacion_ordenada.copy():                   #clona el resto
                        grado = gradoMutacion(clon,rango)
                        matriz, puntuacion = clon
                        clon = aplicar_mutacion_incremento(matriz, grado)
                        clones.append(clon)                                  #integra en poblacion

                
                poblacion_ordenada = poblacion_ordenada + clones
#                print("con Clones", len(poblacion_ordenada))   
                clones.clear()
                

                
                poblacion = poblacion_ordenada
                poblacion = insertaRandomClones(matriz_original, poblacion,rClonInsert)  #inserta random clones
                tfinal = time.time() - tinicial
                print("time", tfinal)
                #imprimir_matriz(poblacion[0], names,"fitness")                          #best fitness in Generacion
                matriz, puntos = poblacion[0]  
                puntuacion, punts = calculaFitness(matriz)
                
                global NFE
                generaReporte(nCorrida, generacion,  len(poblacion_ordenada), punts, fitnessProm, NFE , tfinal)



            # Obtener e imprimir la mejor matriz de la última generación}
            print("\nMejor matriz (última generación):")
            imprimir_matriz( poblacion[0],names, "matriz+fitness")
            imprimir_matriz( poblacion[0],names, "file")
            nCorrida = nCorrida +1
            NFE = 0                                                        #reset de NFE

if __name__ == "__main__":
    main()
