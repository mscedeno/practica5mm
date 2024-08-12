import math
import numpy as np
from matplotlib import pyplot as plt

def eclineal2(A,B):
    A = np.array(A)
    B = np.array(B)
    x = np.linalg.solve(A,B)
    return x

def nr2i(theta3k,theta4k,theta1,theta2,L1,L2,L3,L4,tol=0.01):
    a = -L1*math.cos(theta1)+L2*math.cos(theta2)
    b = -L1*math.sin(theta1)+L2*math.sin(theta2)
    error = tol+1
    while error>tol:
        F1 = a+L3*math.cos(theta3k)-L4*math.cos(theta4k)
        F2 = b+L3*math.sin(theta3k)-L4*math.sin(theta4k)
        J11 = -L3*math.sin(theta3k)
        J12 = +L4*math.sin(theta4k)
        J21 = L3*math.cos(theta3k)
        J22 = -L4*math.cos(theta4k)
        d = eclineal2([[J11,J12],[J21,J22]],[-F1,-F2])
        theta3k += d[0]
        theta4k += d[1]
        error = math.sqrt((d[0]**2)+(d[1]**2))
    return theta3k,theta4k

N = 360
L1 = 15.00 #in
L2 = 4.00 #pulg
L3 = 12.00 #pulg
L4 = 8.00 #pulg

g=386.01#in/s2

m2 = 0.002 #blobs=slug
m3 = 0.02
m4 = 0.10


I2 = 0.10 #blobs-pulg2
I3 = 0.20
I4 = 0.50

Rcg2 = 2 #pulg
Thetarcg2 = 0 #rad
Rcg4 = 4 #pulg radio a centro de gravedad 4
Rc_g4 = 4.95 #pulg radio del centro de gravedad 4 a punto C
thetarcg4 = 30 #grados
thetarcg4 *= math.pi/180 #rad
theta4rc_g4 = 75.52 #grados
theta4rc_g4 *= math.pi/180 #rad
thetarg4=thetarcg4
Rcg3 = 5 #pulg 
thetarcg3 = 0 #grados
Rc_g3 = L3-Rcg3

theta1 = 0 #grados
theta1 *= math.pi/180 #rad
theta2 = np.linspace(2*math.pi,0,N+1) #rad
theta3 = np.zeros(N+1)
theta4 = np.zeros(N+1)
thetap4 = 30 #grados
thetap4 *= math.pi/180 #rad
omega2 = 190 #rpm
omega2 *= math.pi/30 #rad/s2
tRev = 2*math.pi/abs(omega2) #tiempo de un revolucion [s]
t = np.linspace(0,tRev,N+1) #s
omega3 = np.zeros(N+1)
omega4 = np.zeros(N+1)
alpha3 = np.zeros(N+1)
alpha4 = np.zeros(N+1)
alpha2 =0
alpha2 *= math.pi/180 #rad
Vcg2 = np.zeros((N+1,3))
Vcg3 = np.zeros((N+1,3))
Vl2 = np.zeros((N+1,3))
Vcg4 = np.zeros((N+1,3))

Acg2 = np.zeros((N+1,3))
Acg3 = np.zeros((N+1,3))
Acg4 = np.zeros((N+1,3))

Vcg2mod = np.zeros((N+1,1))
Acg2mod = np.zeros((N+1,1))
Vcg3mod = np.zeros((N+1,1))
Acg3mod = np.zeros((N+1,1))
Vcg4mod = np.zeros((N+1,1))
Acg4mod = np.zeros((N+1,1))

variables = ['Tin','Ax','Ay','Bx','By','Cx','Cy','Dx','Dy']
n_variables = len(variables)
Mf = np.zeros((N+1,n_variables),dtype=np.float32)
var_index = {var:idx for idx,var in enumerate(variables)}

for i in range(0,N+1):
    #Calculo de theta 3 y theta 4
    theta3[i], theta4[i] = nr2i(5.6,5,theta1,theta2[i],L1,L2,L3,L4)
    #Calculo de omega 3 y omega 4
    a11 = -L3*math.sin(theta3[i])
    a12 = +L4*math.sin(theta4[i])
    a21 = +L3*math.cos(theta3[i])
    a22 = -L4*math.cos(theta4[i])

    b1 = L2*math.sin(theta2[i])*omega2
    b2 = -L2*math.cos(theta2[i])*omega2

    omega3[i],omega4[i] = eclineal2([[a11,a12],[a21,a22]],[b1,b2])
    #Calculo de alpha 3 y alpha 4
    a11 = -L3*math.sin(theta3[i])
    a12 = L4*math.sin(theta4[i])
    a21 = L3*math.cos(theta3[i])
    a22 = -L4*math.cos(theta4[i])

    b1 = L2*math.cos(theta2[i])*(omega2**2)+L3*math.cos(theta3[i])*(omega3[i]**2)-L4*math.cos(theta4[i])*(omega4[i]**2)
    b2 = L2*math.sin(theta2[i])*(omega2**2)+L3*math.sin(theta3[i])*(omega3[i]**2)-L4*math.sin(theta4[i])*(omega4[i]**2)

    alpha3[i],alpha4[i] = eclineal2([[a11,a12],[a21,a22]],[b1,b2])
    #Calculo de la velocidad lineal en el punto cg2
    omega2Vect = np.array([0,0,omega2])
    R2 = np.array([Rcg2*math.cos(theta2[i]),Rcg2*math.sin(theta2[i]),0])
    Vcg2 = np.cross(omega2Vect,R2)

    #Calculo de la velocidad lineal en el punto cg4
    omega4Vect = np.array([0,0,omega4[i]])
    theta5 = theta4[i] + thetap4
    R4 =  np.array([Rcg4*math.cos(theta5),Rcg4*math.sin(theta5),0])
    Vcg4[i] = np.cross(omega4Vect,R4)

    #Calculo de la velocidad lineal en el punto Vl2
    omega2Vect = np.array([0,0,omega2])
    R2l = np.array([L2*math.cos(theta2[i]),L2*math.sin(theta2[i]),0])
    Vl2 = np.cross(omega2Vect,Vcg2)

    #Calculo de la velocidad lineal en el punto Vcg3
    omega3Vect = np.array([0,0,omega3[i]])
    R3 =  np.array([Rcg3*math.cos(theta3[i]),Rcg3*math.sin(theta3[i]),0])
    Vcg3[i] = Vl2 + np.cross(omega3Vect,R3)
    Vcg3mod[i] = math.sqrt((Vcg3[i,0]**2)+(Vcg3[i,1]**2))

    #Calculo de la aceleracion lineal en el punto Acg2
    Acg2 = np.cross(omega2Vect,Vcg2)

    #Calculo de la aceleracion lineal en el punto Acg4
    alpha4Vect = np.array([0,0,alpha4[i]])

    Acg4[i] = np.cross(alpha4Vect,R4) + np.cross(omega4Vect,Vcg4[i])
    Acg4mod[i] = math.sqrt((Acg4[i,0]**2)+(Acg4[i,1]**2))

    #Calculo de la aceleracion lineal en el punto Al2
    Al2 = np.cross(omega2Vect,Vl2)

    #Calculo de la aceleracion lineal en el punto Acg3
    alpha3Vect = np.array([0,0,alpha3[i]])

    Acg3[i] = np.cross(alpha3Vect,R3) + np.cross(omega3Vect,Vcg3[i])
    Acg3mod[i] = math.sqrt((Acg3[i,0]**2)+(Acg3[i,1]**2))
    
    #Matriz de coeficientes
    matriz_coeficientes = np.zeros((n_variables,n_variables), dtype=np.float32)
    j=0
    matriz_coeficientes[j,var_index['Ax']]=1
    matriz_coeficientes[j,var_index['Bx']]=-1
    j=1
    matriz_coeficientes[j,var_index['Ay']]=1
    matriz_coeficientes[j,var_index['By']]=-1
    j=2
    matriz_coeficientes[j,var_index['Tin']]=1
    matriz_coeficientes[j,var_index['Ax']]=Rcg2*math.sin(theta2[i]+math.pi)
    matriz_coeficientes[j,var_index['Ay']]=-Rcg2*math.cos(theta2[i]+math.pi)
    matriz_coeficientes[j,var_index['Bx']]=Rcg2*math.sin(theta2[i])
    matriz_coeficientes[j,var_index['By']]=-Rcg2*math.cos(theta2[i])
    j=3
    matriz_coeficientes[j,var_index['Bx']]=1
    matriz_coeficientes[j,var_index['Cx']]=-1
    j=4
    matriz_coeficientes[j,var_index['By']]=1
    matriz_coeficientes[j,var_index['Cy']]=-1
    j=5
    matriz_coeficientes[j,var_index['Bx']]=Rcg3*math.sin(theta3[i]+math.pi)
    matriz_coeficientes[j,var_index['By']]=-Rcg3*math.cos(theta3[i]+math.pi)
    matriz_coeficientes[j,var_index['Cx']]=Rc_g3*math.sin(theta3[i])
    matriz_coeficientes[j,var_index['Cy']]=-Rc_g3*math.cos(theta3[i])
    j=6
    matriz_coeficientes[j,var_index['Cx']]=1
    matriz_coeficientes[j,var_index['Dx']]=-1
    j=7
    matriz_coeficientes[j,var_index['Cy']]=1
    matriz_coeficientes[j,var_index['Dy']]=-1
    j=8
    matriz_coeficientes[j,var_index['Cx']]=Rc_g4*math.sin(theta4[i]+theta4rc_g4)
    matriz_coeficientes[j,var_index['Cy']]=-Rc_g4*math.cos(theta4[i]+theta4rc_g4)
    matriz_coeficientes[j,var_index['Dx']]=Rcg4*math.sin(theta4[i]+thetarg4)
    matriz_coeficientes[j,var_index['Dy']]=-Rcg4*math.cos(theta4[i]+thetarg4)
    #Vector de soluciones
    vector_soluciones = np.zeros(n_variables, dtype=np.float32)
    j=0
    vector_soluciones[j]=m2*Acg2[0]
    j=1
    vector_soluciones[j]=m2*Acg2[1]+m2*g
    j=2
    vector_soluciones[j]=I2*alpha2
    j=3
    vector_soluciones[j]=m3*Acg3[i,0]
    j=4
    vector_soluciones[j]=m3*Acg3[i,1]+m3*g
    j=5
    vector_soluciones[j]=I3*alpha3[0]
    j=6
    vector_soluciones[j]=m4*Acg4[i,0]
    j=7
    vector_soluciones[j]=m4*Acg4[i,1]+m4*g
    j=8
    vector_soluciones[j]=I4*alpha4[0]
    #print(vector_soluciones[8])
    #Resolucion 
    vector_incognitas = np.linalg.solve(matriz_coeficientes,vector_soluciones)
    Mf[i] = vector_incognitas
#Grafica
#print('hola')
#print(Mf[:,var_index['Tin']])
plt.plot(t,Mf[:,var_index['Ax']])
plt.plot(t,Mf[:,var_index['Ay']])
plt.title('Fuerza union A')
plt.xlabel('tiempo [s]')
plt.ylabel('Fuerza en A [blob*in/s2]')
plt.legend(['Fuerza en x','Fuerza en y'])
plt.grid(True)
plt.show()

plt.plot(t,Mf[:,var_index['Bx']])
plt.plot(t,Mf[:,var_index['By']])
plt.title('Fuerza union B')
plt.xlabel('tiempo [s]')
plt.ylabel('Fuerza en B [blob*in/s2]')
plt.legend(['Fuerza en x','Fuerza en y'])
plt.grid(True)
plt.show()

plt.plot(t,Mf[:,var_index['Cx']])
plt.plot(t,Mf[:,var_index['Cy']])
plt.title('Fuerza union C')
plt.xlabel('tiempo [s]')
plt.ylabel('Fuerza en A [blob*in/s2]')
plt.legend(['Fuerza en x','Fuerza en y'])
plt.grid(True)
plt.show()

plt.plot(t,Mf[:,var_index['Dx']])
plt.plot(t,Mf[:,var_index['Dy']])
plt.title('Fuerza union D')
plt.xlabel('tiempo [s]')
plt.ylabel('Fuerza en D [blob*in/s2]')
plt.legend(['Fuerza en x','Fuerza en y'])
plt.grid(True)
plt.show()

plt.plot(t,Mf[:,var_index['Tin']])
plt.title('Torque')
plt.xlabel('tiempo [s]')
plt.ylabel('Torque[blob*in2/s2]')
plt.grid(True)
plt.show()
