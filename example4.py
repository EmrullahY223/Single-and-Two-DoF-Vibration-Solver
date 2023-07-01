if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import VibrationCalculator as vibcalc
    import math

    x1_0 = 1 #float(input("Enter initial displacement: "))
    v1_0 = 0.0 #float(input("Enter initial velocity: "))
    a1_0 = 0.0 #float(input("Enter initial acceleration: "))
    x2_0 = 1
    movetype = 1
    baseMag = 0.0
    baseFreq = 0.0
    stepSlope = 0.0
    stepEnd = 0.0
    times = [0]
    v2_0 = 0.0
    a2_0 = 0.0
    type1 = 1
    type2 = 1
    times1 = [0]
    times2 = [0]
    mag1 = 0.0
    mag2 = 0.0 #float(input("Enter force magnitude: "))
    freq1 = 0.0 #float(input("Enter force frequency: "))
    freq2 = 0.0
    m1 = 50.0 #float(input("Enter mass: "))
    m2 = 300.0
    c1 = 0.0 #float(input("Enter damper: "))
    c2 = 1200.0
    k1 = 212.0e3 #float(input("Enter spring: "))
    k2 = 12.0e3
    tfin = 1.5 #float(input("Enter final time: "))
    num_of_division = 100

    num_of_steps,dt = vibcalc.step_size(k1 ,m1 ,k2 ,m2 ,tfin,num_of_division)
    F1,F2 = vibcalc.forceCalc(type1,type2,mag1,mag2,freq1 ,freq2 ,times1 ,times2 ,num_of_steps,dt)
    Base,BaseSpeed = vibcalc.MovingBase(baseMag,baseFreq,dt,num_of_steps,m2,movetype,stepSlope,stepEnd,times)
    X1,X2,V1,V2,A1,A2,G1,G2,t = vibcalc.Solver(x1_0 ,v1_0,a1_0 ,x2_0 ,v2_0 ,a2_0 ,F1,F2,m1,m2 ,c1 ,c2 ,k1 ,k2,num_of_steps,dt,Base,BaseSpeed)

    plt.title("Suspension Example")
    plt.plot(t,X2,label = 'Numerical Solution x1')
    plt.xlabel("time(sec)")
    plt.ylabel("Displacement(m)")
    plt.plot(t,X1,label = 'Numerical Solution x2',ls = ':',color = 'r')
    plt.legend()
    plt.show()

    exactSolution = list()