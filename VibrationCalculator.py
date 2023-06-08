"""
Project 3 kod dosyasıdır
numerik yöntemler ile titreşim analizi için yazılmıştır.
"""
import math
def Solver(x1_0,v1_0,a1_0,x2_0,v2_0,a2_0,mag1,mag2,forceFreq1,forceFreq2,m1,m2,c1,c2,k1,k2,tfin):
    t=[0]
    dt = float(1/(((k1/m1)**0.5)/(2*math.pi))/120)
    print(f"dt : {dt}")
    num_of_steps = round((tfin/dt)+1)
    print(f"number of steps: {num_of_steps}")
    X1 = [0] *(num_of_steps+1)
    V1 = [0] *(num_of_steps+1)
    A1 = [0] *(num_of_steps+1)
    G1 = [0] *(num_of_steps+1)

    X2 = [0] *(num_of_steps+1)
    V2 = [0] *(num_of_steps+1)
    A2 = [0] *(num_of_steps+1)
    G2 = [0] *(num_of_steps+1)
    
    X1[0] = x1_0
    V1[0] = v1_0
    A1[0] = a1_0

    X2[0] = x2_0
    V2[0] = v2_0
    A2[0] = a2_0

    # i = 0
    t.append(t[0]+dt)
    # Dof 1
    if m1 > 0:
        G1[0] = (mag1*math.sin(forceFreq1*t[0])) - k2*(X1[0]-X2[0]) - c2*(V1[0]-V2[0])
        X1[1] = X1[0] + dt*V1[0] + ((dt**2)/2)*A1[0]
        V1[1] = (-2*m1*X1[0] + 2*m1*X1[1] + G1[0]*dt**2 - k1*dt**2*X1[1])/(2*m1*dt + c1*dt**2)
        A1[1] = (G1[0] - c1*V1[1] - k1*X1[1])/m1
        #print(f"g1 :{G1[0]}\nX1: {X1[1]}\nV1: {V1[1]}\nA1: {A1[1]}")
    # Dof 2
    if m2 > 0:
        G2[0] = mag2*math.sin(forceFreq2*t[0]) + k2*(X1[0]) + c2*(V1[0])
        X2[1] = X2[0] + dt*V2[0] + ((dt**2)/2)*A2[0]
        V2[1] = (-2*m2*X2[0] + 2*m2*X2[1] + G2[0]*dt**2 - k2*dt**2*X2[1])/(2*m2*dt + c2*dt**2)
        A2[1] = (G2[0] - c2*V2[1] - k2*X2[1])/m2
        #print(f"g1 :{G2[0]}\nX1: {X2[1]}\nV1: {V2[1]}\nA1: {A2[1]}")

    t.append(t[1]+dt)
    # i = 1
    # Dof 1
    if m1 > 0:
        G1[1] = mag1*math.sin(forceFreq1*t[1]) - k2*(X1[1]-X2[1]) - c2*(V1[1]-V2[1])
        X1[2] = 2*X1[1] - X1[0] + dt**2*A1[1]
        V1[2] = (-2*m1*X1[1] + 2*m1*X1[2] + G1[1]*dt**2 - k1*dt**2*X1[2])/(2*m1*dt + c1*dt**2)
        A1[2] = (G1[1] - c1*V1[2] - k1*X1[2])/m1
        #print(f"g1 :{G1[1]}\nX1: {X1[2]}\nV1: {V1[2]}\nA1: {A1[2]}")
    # Dof 2
    if m2 >0:
        G2[1] = mag2*math.sin(forceFreq2*t[1]) + k2*(X1[1]) + c2*(V1[1])
        #print(f"G2[1]: {G2[1]}\nforce: {mag2*math.sin(forceFreq2*t[1])}")
        X2[2] = 2*X2[1] - X2[0] + dt**2*A2[1]
        V2[2] = (-2*m2*X2[1] + 2*m2*X2[2] + G2[1]*dt**2 - k2*dt**2*X2[2])/(2*m2*dt + c2*dt**2)
        A2[2] = (G2[1] - c2*V2[2] - k2*X2[2])/m2
        #print(f"A2 : {A2[2]}")
        #print(f"g1 :{G2[1]}\nX1: {X2[1]}\nV1: {V2[1]}\nA1: {A2[1]}")
  
    
    # loop for numerical analysis
    for i in range(2,num_of_steps):
        #print (f"i={i}")
        t.append(t[i]+dt)
        if m1 > 0:
            G1[i] = mag1*math.sin(forceFreq1*t[i]) - (k2*(X1[i]-X2[i])) - (c2*(V1[i]-V2[i]))
            #print(f"force1:{mag1*math.sin(forceFreq1*t[i])}")
            X1[i+1] = (2*X1[i]) - (X1[i-1]) + ((dt**2)*A1[i])
            V1[i+1] = ((-2*m1*X1[i]) + (2*m1*X1[i+1]) + (G1[i]*dt**2) - (k1*(dt**2)*X1[i+1]))/((2*m1*dt) + (c1*dt**2))
            A1[i+1] = (G1[i] - (c1*V1[i+1]) - (k1*X1[i+1]))/m1
        if m2 > 0:
            G2[i] = mag2*math.sin(forceFreq2*t[i]) + (k2*(X1[i])) + (c2*(V1[i]))
            #print(f"force2:{mag2*math.sin(forceFreq2*t[i])}")
            X2[i+1] = 2*X2[i] - X2[i-1] + dt**2*A2[i]
            V2[i+1] = (-2*m2*X2[i] + 2*m2*X2[i+1] + G2[i]*dt**2 - k2*dt**2*X2[i+1])/(2*m2*dt + c2*dt**2)
            A2[i+1] = (G2[i] - c2*V2[i+1] - k2*X2[i+1])/m2
        
    
    return X1,X2,V1,V2,A1,A2,G1,G2,t



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    
     

    x1_0 = 0.0 #float(input("Enter initial displacement: "))
    v1_0 = 0.0 #float(input("Enter initial velocity: "))
    a1_0 = 0.0 #float(input("Enter initial acceleration: "))
    x2_0 = 0.0
    v2_0 = 0.0
    a2_0 = 0.0
    mag1 = 0.0
    mag2 = 20.0 #float(input("Enter force magnitude: "))
    forceFreq1 = 0.0 #float(input("Enter force frequency: "))
    forceFreq2 = 10.0
    m1 = 20.0 #float(input("Enter mass: "))
    m2 = 40.0
    c1 = 0.0 #float(input("Enter damper: "))
    c2 = 0.0
    k1 = 4000.0 #float(input("Enter spring: "))
    k2 = 2000.0
    tfin = 10.0 #float(input("Enter final time: "))

    X1,X2,V1,V2,A1,A2,G1,G2,t = Solver(x1_0,v1_0,a1_0,x2_0,v2_0,a2_0,mag1,mag2,forceFreq1,forceFreq2,m1,m2,c1,c2,k1,k2,tfin)
    #print(X2)

    
    plt.plot(t,X1)
    
    if m2>0:
        plt.plot(t,X2)
    else:
        print("no second Dof")
    plt.show()

