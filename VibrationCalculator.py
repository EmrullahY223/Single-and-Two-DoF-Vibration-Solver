"""
Project 3 kod dosyasıdır
numerik yöntemler ile titreşim analizi için yazılmıştır.
"""
import math
def step_size(k1,m1,tfin):

    dt = float(1/(((k1/m1)**0.5)/(2*math.pi))/120)
    print(f"dt : {dt}")
    num_of_steps = round((tfin/dt)+1)
    return num_of_steps, dt
def forceCalc(type1,type2,mag1,mag2,freq1,freq2,times1,times2,num_of_steps,dt):
    F1 = [0] * (num_of_steps+1)
    F2 = [0] * (num_of_steps+1)
    t  = [0] 
    num_of_stepsBetweenSeconds = int(1/dt)
    for i in range(0,num_of_steps):
        t.append(t[i]+dt)
    if type1 == 1:
        for i in range(0,len(F1)):
            F1[i] = mag1*math.sin(freq1*t[i])
            
    elif type1 ==2:
        for i in times1:
            print(i)
            F1[round(i*num_of_stepsBetweenSeconds)] = mag1/dt
    else:
        print("dafuk")
    if type2 == 1:
        for i in range(0,len(F1)):
            F2[i] = mag2*math.sin(freq2*t[i])
            
    elif type2 ==2:
        for i in times2:
            F2[i*num_of_stepsBetweenSeconds] = mag2/dt
    else:
        print("dafuk2")
    return F1,F2
    
def Solver(x1_0,v1_0,a1_0,x2_0,v2_0,a2_0,F1,F2,m1,m2,c1,c2,k1,k2,num_of_steps):
    
    t  = [0]
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
    
    # Dof 1
    if m1 > 0:
        G1[0] = (F1[0]) - k2*(X1[0]-X2[0]) - c2*(V1[0]-V2[0])
        X1[1] = X1[0] + dt*V1[0] + ((dt**2)/2)*A1[0]
        V1[1] = (-2*m1*X1[0] + 2*m1*X1[1] + G1[0]*dt**2 - k1*dt**2*X1[1])/(2*m1*dt + c1*dt**2)
        A1[1] = (G1[0] - c1*V1[1] - k1*X1[1])/m1
        
    # Dof 2
    if m2 > 0:
        G2[0] = F2[0] + k2*(X1[0]) + c2*(V1[0])
        X2[1] = X2[0] + dt*V2[0] + ((dt**2)/2)*A2[0]
        V2[1] = (-2*m2*X2[0] + 2*m2*X2[1] + G2[0]*dt**2 - k2*dt**2*X2[1])/(2*m2*dt + c2*dt**2)
        A2[1] = (G2[0] - c2*V2[1] - k2*X2[1])/m2
    t.append(t[0]+dt)    
    # loop for numerical analysis
    for i in range(1,num_of_steps):
        #print (f"i={i}")
        t.append(t[i]+dt)
        if m1 > 0:
            G1[i] = F1[i] - (k2*(X1[i]-X2[i])) - (c2*(V1[i]-V2[i]))
            #print(f"force1:{mag1*math.sin(forceFreq1*t[i])}")
            X1[i+1] = (2*X1[i]) - (X1[i-1]) + ((dt**2)*A1[i])
            V1[i+1] = ((-2*m1*X1[i]) + (2*m1*X1[i+1]) + (G1[i]*dt**2) - (k1*(dt**2)*X1[i+1]))/((2*m1*dt) + (c1*dt**2))
            A1[i+1] = (G1[i] - (c1*V1[i+1]) - (k1*X1[i+1]))/m1
        if m2 > 0:
            G2[i] = F2[i] + (k2*(X1[i])) + (c2*(V1[i]))
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
    type1 = 2
    type2 = 2
    times1 = [0, 1.5]
    times2 = [0]
    mag1 = 5.0
    mag2 = 0.0 #float(input("Enter force magnitude: "))
    freq1 = 0.0 #float(input("Enter force frequency: "))
    freq2 = 0.0
    m1 = 10.0 #float(input("Enter mass: "))
    m2 = 0.0
    c1 = 89.40 #float(input("Enter damper: "))
    c2 = 0.0
    k1 = 20000.0 #float(input("Enter spring: "))
    k2 = 0.0
    tfin = 3.0 #float(input("Enter final time: "))

    num_of_steps,dt = step_size(k1,m1,tfin)
    F1,F2 = forceCalc(type1,type2,mag1,mag2,freq1,freq2,times1,times2,num_of_steps,dt)
    X1,X2,V1,V2,A1,A2,G1,G2,t = Solver(x1_0,v1_0,a1_0,x2_0,v2_0,a2_0,F1,F2,m1,m2,c1,c2,k1,k2,num_of_steps)
    #print(F1)

    
    plt.plot(t,X1)
    
    if m2>0:
        plt.plot(t,X2)
    else:
        print("no second Dof")
    plt.show()

