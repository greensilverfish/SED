#    Description:
#    ___________
#
#       This program takes observed magnitudes and spectral types and
#       calculates the extinction from the Rc-Ic color index returning a
#       dereddened observed SED, that expected from a star with input spectral 
#       type, as well as a face-on reprocessing disk model from HSKV (1992). 
#       
#       Watts cm^-2 mu^-1 were taken from LAH as follows - 
# 
#       UBVRcIc from Bessell 1979, AJ, V91, 589. 
#       JHKLM from Campins, Rieke, and Lebofsky 1985, AJ V90 896. 
#       NQ from Rieke, Lebofsky, and Low 1985, AJ, V90 904. 
#       and IRAS12,25,60, and 100 from IRAS Exp. Supp. (color-corrected!). 
# 
#       Intrinsic stellar SEDs were taken from - 
#
#       Teff B8-K5 from Schmitd-Kaler.  Teff K5-M6 from Bessell (1991, AJ, V101, 662).
#       J-H, H-K, V-K from Bessell and Brett (1988, PASP, V100, 1134).  
#       >K5 R-I, V-I Bessell (1991) <K5 R-I, V-I from Johnson (1966, ARAA, V4, 193).
#       With color corrections Johnson to Cousins from Besell (1979, AJ, V91, 589). 
#       B-V from Johnson B8-K5 Bessell (1991) K7-M6.
#       U-V from Johnson B8-M6.
#       (V-12) from Waters et al. (1987) A&A,  V172, 225.
#
#       Extinction law taken from Rieke and Lebofsky (ApJ, 1985, V288, 618).
#       
#       Reprocessing disk models were taken from Hillenbrand et al. (1992, APJ, 397, 613)
#       and presented here for spectral types A0, F0, G0, K0, K5, M0, and M3.
#  
#    Call Sequence:
#    ______________
#
#    > cc bigsed.c -o bigsed  -lm 
#
#    Parameters:
#    ___________
#
#       Input stellar effective temperature, distance [only for Log(L#)], 
#       U,B,V,Rc,Ic,J,H,K,IRAC-3.6,IRAC-4.5,IRAC-5.8,IRAC-8,MIPS-24,MIPS-70,60, and 100  IN MAGNITUDES!  In order
#       to convert Janskies to magnitudes, need flux for zero-mag star: 
#       
#       Band        U     B      V          Rc        Ic        J  
#       Fnu(Jy)     1810  4260   3640       3080      2550      1603  
#       Band        H     K    IRAC-3.6    IRAC-4.5  IRAC-5.8   IRAC-8
#       Fnu(Jy)     1075  667   280.9        179.7     115      64.13 
#       Band        MIPS-24    MIPS-70    IRAS60    IRAS100 
#       Fnu(Jy)     7.17       0.778      1.19      0.43  
#
#    Returns:
#    ________
#
#       Dereddened observed SED as well as expected stellar SED based on 
#       spectral type, normalized to match at I-band.
#
#    Notes:
#    ______
#
#    By:
#    ___
#
#    MRM 2.10.96     MPIA-Heidelberg
#       MRM 3-16-99     Steward Observatory - zero-points modified. 
#
###############################################################################/
import math

# Set constants   #/

# Units of Watts cm^-2 mu^-1 for a zero-magnitude star #/

flux = [4.19e-12, 6.59e-12, 3.60e-12, 2.25e-12, 1.23e-12, 3.129e-13, 1.133e-13, 4.283e-14,6.50e-15, 2.66e-15, 1.02e-15, 3.04e-16, 3.73e-18, 4.76e-20, 9.92e-20, 1.29e-20]

# Wavelength scale adopted U, B, V, Rc, Ic, J, H, K, L, M, N, IRAS12, Q, IRAS25, IRAS60, IRAS100

lam = [0.36, 0.44, 0.55, 0.64, 0.79, 1.235, 1.662, 2.159, 3.55, 4.49, 5.73, 7.87, 24.0, 70.0, 60.0, 100.0]

#  Extinction law adopted as Alam/Av from Reike and Lebofsky UBV, Rc, Ic, JHKLMN12  #/ 

extinc = [1.53, 1.32, 1.0, 0.82, 0.60, 0.333, 0.207, 0.133, 0.084, 0.071, 0.064, 0.065, 0.064, 0.00, 0.00, 0.00]

#  Read in source, dm, temp, U, B, V, Rc, Ic, Jc, Hc, Kc, L, M, N, 12, Q, 25, 60, 100 #/

file=open("/Users/Kryss/Desktop/sedinput.txt")
for line in file:    
    print line.replace("\n","").split("\t")
    row = line.replace("\n","").split("\t")
    print row
    name=row[0]
    lteff= float(row[1])
    dm= float(row[2])
    av= float(row[3])
    obs=[float(row[4]), float(row[5]), float(row[6]), float(row[7]), float(row[8]), float(row[9]), float(row[10]), float(row[11]), float(row[12]), float(row[13]), float(row[14]), float(row[15]), float(row[16]), float(row[17]), float(row[18]), float(row[19])]
    teff = math.pow(10.0,lteff)
    x = teff
# Bstars giving divsion by zero error so we are skiping them when we read in.  
    if name == "XMM1829561_0100217" or name == "HD_170634": 
        continue
    
    if (teff > 1.0):
    
    #  Use the derived fits (references in header).
    
        ub = 2.674702e1 - x*3.675133e-2 + x*x*1.996307e-5 - x*x*x*5.315118e-9 + x*x*x*x*7.388618e-13 - x*x*x*x*x*5.155066e-17 + x*x*x*x*x*x*1.426480e-21
        bv = 9.026432e1 - 1.116759e-1*x + x*x*5.917468e-5 -  x*x*x*1.703182e-8 + x*x*x*x*2.867242e-12 - x*x*x*x*x*2.826478e-16 + x*x*x*x*x*x*1.513376e-20 - x*x*x*x*x*x*x*3.400976e-25
        bcv = -101.3651 + x*8.971698e-2 + x*x*-3.290524e-5 + x*x*x*6.334656e-9 + x*x*x*x*-6.714475e-13 + x*x*x*x*x*3.705377e-17 + x*x*x*x*x*x*-8.309852e-22
        rico = 2.195985e1 - 1.499330e-2*x + x*x*4.177120e-6 -  x*x*x*5.805917e-10 + x*x*x*x*3.992331e-14 - x*x*x*x*x*1.084173e-18
        vico = 1.103107e2 - 1.147069e-1*x + x*x*5.191987e-5 -  x*x*x*1.296276e-8 + x*x*x*x*1.911211e-12 - x*x*x*x*x*1.657877e-16 + x*x*x*x*x*x*7.820634e-21 - x*x*x*x*x*x*x*1.546372e-25
        vko = 231.6392 - 0.2627116*x + x*x*1.32161e-4 - x*x*x*3.773421e-8 + x*x*x*x*6.658074e-12 - x*x*x*x*x*7.432943e-16 + x*x*x*x*x*x*5.129903e-20 - x*x*x*x*x*x*x*2.001820e-24 + x*x*x*x*x*x*x*x*3.381241e-29
        jho = 2.857311e2 + x*-4.551107e-1 + x*x*3.134182e-4 + x*x*x*-1.224157e-7 + x*x*x*x*2.994256e-11 + x*x*x*x*x*-4.764616e-15 + x*x*x*x*x*x*4.938796e-19 + x*x*x*x*x*x*x*-3.218916e-23 + x*x*x*x*x*x*x*x*1.197969e-27 + x*x*x*x*x*x*x*x*x*-1.940957e-32
        hko = 2.788605e0 + x*-1.731798e-3 + x*x*4.424067e-7 + x*x*x*-5.647841e-11 + x*x*x*x*3.559478e-15 + x*x*x*x*x*-8.828952e-20
        klo = 4.964018e0 - x*3.818125e-3 + x*x*1.191141e-6 - x*x*x*1.846119e-10 + x*x*x*x*1.411467e-14 - x*x*x*x*x*4.248989e-19
        v12 =  1.257185e1 - x*3.235974e-3 + x*x*2.866275e-7 - x*x*x*8.819869e-12
    #  Perform color transformations
    
        jhco =  jho*0.911 
        hkco =  hko*0.971 - 0.02
            
    #  Calculate Colors V-R
    
    #    ric = obs[2] - obs[3]
    #    print "ric=" ,ric
    
    #  Calculate Extinction
    
    #    av = 5*(ric - (vico-rico))
    #    print "av=" ,av
    
    #  Calculate absolute and dereddened Magnitudes for A0 and cooler used .76 for hotter than A0 used .82... mio is really mro, and Mi is really Mr
    
        mio = obs[3] - 0.76*av
        Mi = obs[3] - 0.76*av - dm
    
    # Calculate luminosities bcic is really bcrc, mboli is really mbolr, and v-r = vico -rico
    
        bcic = bcv + (vico-rico)
        mboli = Mi + bcic 
        llumi = 1.89 - 0.4*mboli
    
    # Calculate stellar SED normalized at R-band  - NOTE DISTANCE NOT USED! - 

        stellar = [ub + bv + (vico - rico) + mio]
        stellar.append(bv + (vico -rico) + mio)
        stellar.append((vico- rico) + mio)
        stellar.append(mio)
        stellar.append(mio- rico)
        stellar.append(stellar[2] - vko + hkco + jhco)
        stellar.append(stellar[5] - jhco)
        stellar.append(stellar[6] - hkco) 
        stellar.append(stellar[7] - klo) 
        stellar.append (stellar[8]) 
        stellar.append(stellar[9])  
        stellar.append(stellar[10])
        stellar.append(stellar[11]) 
        stellar.append(stellar[12]) 
        stellar.append(stellar[13]) 
        stellar.append(stellar[14])
    
    # Calculate star + face-on reprocessing disk from Hillenbrand et al. (1992) APJ V397 613
    
    if (12000.0 > teff and teff > 8600.0):
     
        repro= [stellar[0]]
        repro.append(stellar[1])
        repro.append(stellar[2])
        repro.append(stellar[3] - 0.49)
        repro.append(stellar[4] - 0.72)
        repro.append(stellar[5] - 1.09) 
        repro.append(stellar[6] - 1.41) 
        repro.append(stellar[7] - 1.86) 
        repro.append(stellar[8] - 2.60) 
        repro.append(stellar[9] - 3.09) 
        repro.append(stellar[10] - 3.48) 
        repro.append(stellar[11] - 4.00) 
        repro.append(stellar[12] - 5.98) 
        repro.append(stellar[13] - 8.01) 
        repro.append(stellar[14] - 7.51) 
        repro.append(stellar[15] - 8.43)
    
    
    if 8600.0 > teff and teff > 6400.0 :
    
        repro = [stellar[0]] 
        repro.append(stellar[1]) 
        repro.append(stellar[2]) 
        repro.append(stellar[3] - 0.26) 
        repro.append(stellar[4] - 0.42) 
        repro.append(stellar[5] - 0.72) 
        repro.append(stellar[6] - 0.99) 
        repro.append(stellar[7] - 1.39) 
        repro.append(stellar[8] - 2.07) 
        repro.append(stellar[9] - 2.54) 
        repro.append(stellar[10] - 2.93) 
        repro.append(stellar[11] - 3.45) 
        repro.append(stellar[12] - 5.43) 
        repro.append(stellar[13] - 7.46) 
        repro.append(stellar[14] - 6.91)
        repro.append(stellar[15] - 7.83)
    
    if (6400.0 > teff and teff > 5800.0):
        repro = [stellar[0]] 
        repro.append(stellar[1]) 
        repro.append(stellar[2]) 
        repro.append(stellar[3] - 0.17) 
        repro.append(stellar[4] - 0.30) 
        repro.append(stellar[5] - 0.55) 
        repro.append(stellar[6] - 0.79) 
        repro.append(stellar[7] - 1.16) 
        repro.append(stellar[8] - 1.81) 
        repro.append(stellar[9] - 2.27) 
        repro.append(stellar[10] - 2.66) 
        repro.append(stellar[11] - 3.18) 
        repro.append(stellar[12] - 5.16) 
        repro.append(stellar[13] - 7.19) 
        repro.append(stellar[14] - 6.60) 
        repro.append(stellar[15] - 7.51)
    
    
    if (5800.0 > teff and teff > 5000.0): 
    
        repro = [stellar[0]]
        repro.append(stellar[1]) 
        repro.append(stellar[2]) 
        repro.append(stellar[3] - 0.13) 
        repro.append(stellar[4] - 0.23) 
        repro.append(stellar[5] - 0.45) 
        repro.append(stellar[6] - 0.67) 
        repro.append(stellar[7] - 1.01) 
        repro.append(stellar[8] - 1.64) 
        repro.append(stellar[9] - 2.08) 
        repro.append(stellar[10] - 2.47) 
        repro.append(stellar[11] - 2.99) 
        repro.append(stellar[12] - 4.97) 
        repro.append(stellar[13] - 7.00) 
        repro.append(stellar[14] - 6.38) 
        repro.append(stellar[15] - 7.30)
    
    
    if (5000.0 > teff and teff > 4000.0):
     
        repro = [stellar[0]]
        repro.append(stellar[1]) 
        repro.append(stellar[2]) 
        repro.append(stellar[3] - 0.08) 
        repro.append(stellar[4] - 0.15) 
        repro.append(stellar[5] - 0.32) 
        repro.append(stellar[6] - 0.51) 
        repro.append(stellar[7] - 0.82) 
        repro.append(stellar[8] - 1.40) 
        repro.append(stellar[9] - 1.82) 
        repro.append(stellar[10] - 2.21) 
        repro.append(stellar[11] - 2.73) 
        repro.append(stellar[12] - 4.71) 
        repro.append(stellar[13] - 6.74) 
        repro.append(stellar[14] - 6.07) 
        repro.append(stellar[15] - 6.98)
    
    
    if (4000.0 > teff and teff > 3600.0):
     
        repro = [stellar[0]] 
        repro.append(stellar[1]) 
        repro.append(stellar[2]) 
        repro.append(stellar[3] - 0.05) 
        repro.append(stellar[4] - 0.11) 
        repro.append(stellar[5] - 0.26) 
        repro.append(stellar[6] - 0.42) 
        repro.append(stellar[7] - 0.70) 
        repro.append(stellar[8] - 1.25) 
        repro.append(stellar[9] - 1.66) 
        repro.append(stellar[10] - 2.05) 
        repro.append(stellar[11] - 2.57) 
        repro.append(stellar[12] - 4.55) 
        repro.append(stellar[13] - 6.58) 
        repro.append(stellar[14] - 5.87) 
        repro.append(stellar[15] - 6.78)
    
    if (3600.0 > teff and teff > 3000.0):
    
        repro = [stellar[0]] 
        repro.append(stellar[1]) 
        repro.append(stellar[2]) 
        repro.append(stellar[3] - 0.03) 
        repro.append(stellar[4] - 0.07) 
        repro.append(stellar[5] - 0.18) 
        repro.append(stellar[6] - 0.31) 
        repro.append(stellar[7] - 0.55) 
        repro.append(stellar[8] - 1.06) 
        repro.append(stellar[9] - 1.44) 
        repro.append(stellar[10] - 1.83) 
        repro.append(stellar[11] - 2.35) 
        repro.append(stellar[12] - 4.33) 
        repro.append(stellar[13] - 6.36) 
        repro.append(stellar[14] - 5.59) 
        repro.append(stellar[15] - 6.50)
    
    
    #  Calculate the observed and expected SED  and print out the results.
    
    #  Print out the header
     
    print(name)
    print(lteff, llumi, av, rico)
    outfile= open("/Users/Kryss/Desktop/SED/"+name+"out.txt", 'w+')
    outfile.write("**************************************" + name + "\n")
    outfile.write(str(lteff) + "\t"+str(llumi)+"\t"+str(rico)+"\t" + "\n")
    outfile.write("***************************************************************    \n")
    outfile.write("Lambda           Fobs          deredF                    F*                   Frepro    \n")
    outfile.write("======           ======             ======                =====            =======    \n")
    
    
    for i in range(0,16):
        obsSED = flux[i]/ math.pow(10.0, (obs[i]/2.5))
        derobsSED = flux[i]/ math.pow(10.0,((obs[i] - av*extinc[i])/2.5))
        starSED = flux[i]/ math.pow(10.0,(stellar[i]/2.5))
        reproSED = flux[i]/ math.pow(10.0,(repro[i]/2.5))
        valueout0 = lam[i]
        valueout1 = lam[i]*obsSED
        valueout2 = lam[i]*derobsSED
        valueout3 = lam[i]*starSED
        valueout4 = lam[i]*reproSED

        outfile.write(str(valueout0)+"\t"+str(valueout1)+"\t"+str(valueout2)+"\t"+str(valueout3)+"\t"+str(valueout4)+"\n")
        
    outfile.close()    
file.close()
    
 
