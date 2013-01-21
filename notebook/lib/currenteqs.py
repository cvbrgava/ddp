import re


def current(information,state,vth):
    """current equations verified for level 1 PMOS and NMOS """
    d=state[information[1]]
    g=state[information[2]]
    s=state[information[3]]
    p=re.compile('(l|w|=|u)')
    w=float(p.sub('',information[5]))
    l=float(p.sub('',information[6]))
    vth=float(vth)


    
    if(information[4]=='CMOSP'):
        if(information[0]=='c'):
            return 0
        elif(information[0]=='s'):
            cc=(43e-6) * (w / l) * ((g**2) + (1 - 2 * vth * 0.09) * (s**2) + (vth**2) + ((0.09 * vth - 1) * 2 * g * s) + (2 * g * vth) + (-2 * vth + 0.09 * vth**2) * s - (0.09 * d * vth**2) - (0.09 * d * g**2) -  (0.09 * d * s**2) + ( 2 * g * s * 0.09 * d) - (2 * g * vth * 0.09 * d) + ( 2* s * vth * 0.09 * d ) + (0.09 * s * g**2) + (0.09 * s**3) - (0.09 * 2 * g * s**2))
            return cc
        elif(information[0]=='l'):
            cc=(86e-6) * (w / l) * ((0.5 * s**2) - (0.5 * d**2) - (g * s) - (s * vth) + (g * d) + (d * vth) + (0.09 * 0.5 * s**3) - (0.09 * s * 0.5 * d**2) - (0.09 * g * s**2) - (s**2 * 0.09 * vth) + (2 * 0.09 * s * g * d) + (2 * 0.09 * s * d * vth) - (0.5 * d * 0.09 * s**2) + (0.5 * 0.09 * d**3) - (0.09 * g * d**2) - (0.09 * d**2 *vth))
            return  cc
        
    if(information[4]=='CMOSN'):
        if(information[0]=='c'):
            return 0
        elif(information[0]=='s'):
            cc = (0.5 * 285e-6) * (w / l) * ((g**2) + (1-2 * vth * 0.14) * (s**2) + (vth**2) + ((0.14 * vth - 1) * 2 * g * s) - (2 * g * vth) + (2 * vth - 0.14 * vth**2) * s + (0.14 * d * vth**2) + (0.14 * d * g**2) + (0.14 * d * s**2) - (2 * g * s * 0.14 * d) - (2 * g * vth * 0.14 * d) + (2 * s * vth * 0.14 * d) - (0.14 * s * g**2) - (0.14 * s**3) + (0.14 * 2 * g * s**2))
            return cc
        elif(information[0]=='l'):
            cc=(285e-6) * (w / l) * ((g * d) - (s * g) + (0.5 * s**2) - (vth * d) + (vth * s) - (0.5 * d**2) + (0.14 * g * d**2) - (2 * 0.14 * d * s * g) + (0.5 * 0.14 * d * s**2) - (0.14 * vth * d**2) + (2 * 0.14 * d * vth * s) - (0.5 * 0.14 * d**3) + (s**2 *g * 0.14) - (0.5 * 0.14 * s**3) - (0.14 * s**2 * vth) + (0.5 * 0.14 * s * d**2))
            return cc
