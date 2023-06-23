"""
Project 3 kod dosyasıdır
numerik yöntemler ile titreşim analizi için yazılmıştır.
"""
import math
def step_size(k1,m1,k2,m2,tfin):
    k = 0
    m = 0

    if k2 > k1:
        k = k2
    else:
        k = k1
    if m2 > m1:
        m = m2
    else:
        m = m1

    dt = float(1/(((k/m)**0.5)/(2*math.pi))/50)
    print(f"dt : {dt}")
    num_of_steps = round((tfin/dt)+1)
    return num_of_steps, dt

def MovingBase(baseMag,baseFreq,dt,num_of_steps,m2,movetype,stepSlope,stepEnd,times):
    Base = [0]*(num_of_steps+1)
    BaseSpeed = [0]*(num_of_steps+1)
    t  = [0]
    num_of_stepsBetweenSeconds = int(1/dt)
    if m2 <= 0:
        if movetype == 1:
            for i in range(0,num_of_steps):
                Base[i] = baseMag*math.sin(baseFreq*t[i])
                BaseSpeed[i] = baseMag*baseFreq*math.cos(baseFreq*t[i])
                t.append(t[i]+dt)
        elif movetype == 2:
            
            y = 0 * num_of_stepsBetweenSeconds
            for i in range(y,num_of_steps+1):
                Base[i] = stepSlope * t[i]
                print(i)
                BaseSpeed[i-1] = (Base[i]-Base[i-1])/dt
                t.append(t[i]+dt)
                if Base[i] > stepEnd:
                    break
        else:
            print(f"dafuk mate")
    else:
        print(f"can't have base mass with moving base!!")

    return Base,BaseSpeed

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
    
def Solver(x1_0,v1_0,a1_0,x2_0,v2_0,a2_0,F1,F2,m1,m2,c1,c2,k1,k2,num_of_steps,Base,BaseSpeed):
    
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
    if m2 <= 0:
        X2 = Base
        V2 = BaseSpeed
    
    elif m2 > 0:
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
    movetype = 1
    baseMag = 0.0
    baseFreq = 0.0
    stepSlope = 0.0
    stepEnd = 0.0
    times = [0]
    v2_0 = 0.0
    a2_0 = 0.0
    type1 = 1
    type2 = 2
    times1 = [0]
    times2 = [0]
    mag1 = 10.0
    mag2 = 0.0 #float(input("Enter force magnitude: "))
    freq1 = 82.0 #float(input("Enter force frequency: "))
    freq2 = 0.0
    m1 = 10.0 #float(input("Enter mass: "))
    m2 = 0.0
    c1 = 0.0 #float(input("Enter damper: "))
    c2 = 0.0
    k1 = 64000.0 #float(input("Enter spring: "))
    k2 = 0.0
    tfin = 10.0 #float(input("Enter final time: "))

    num_of_steps,dt = step_size(k1,m1,k2,m2,tfin)
    F1,F2 = forceCalc(type1,type2,mag1,mag2,freq1,freq2,times1,times2,num_of_steps,dt)
    Base,BaseSpeed = MovingBase(baseMag,baseFreq,dt,num_of_steps,m2,movetype,stepSlope,stepEnd,times)
    X1,X2,V1,V2,A1,A2,G1,G2,t = Solver(x1_0,v1_0,a1_0,x2_0,v2_0,a2_0,F1,F2,m1,m2,c1,c2,k1,k2,num_of_steps,Base,BaseSpeed)
    #print(F1)

    
    plt.plot(t,X1)
    #plt.plot(t,Base)
    if m2 <=0 :
        print("no second Dof/moving Base")
    plt.show()

