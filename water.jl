#main helmholz source: http://www.iapws.org/relguide/IAPWS95-2018.pdf
#auxiliary functions source: #http://www.teos-10.org/pubs/Wagner_and_Pruss_2002.pdf
using ForwardDiff

function _p0exp(T) #Vapor–pressure equation, eq 2.5, #SI units

    d=1-T/647.096
    a = [-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502]
    return 22.064e6*exp((647.096/T)*(a[1]*d + a[2]*d^1.5+ a[3]*d^3+ a[4]*d^3.5+ a[5]*d^4+ a[6]*d^7.5))
end
function _dp0dTexp(T) #Vapor–pressure equation derivative, eq 2.5a, #SI units
    d=1-T/647.096
    a = [-7.85951783,1.84408259,-11.7866497,22.6807411,-15.9618719,1.80122502]
p0 = _p0(T)
return-(p0/T)*(log(p0/22.064e6)+ a[1] + 1.5*a[2]*d^0.5+ 3*a[3]*d^2+ 3.5*a[4]*d^2.5+ 4*a[5]*d^3+ 7.5*a[6]*d^6.5)
end

function _rholsatexp(T) #Saturated liquid density equation, eq 2.6, #SI units
    d=1-T/647.096
    b = [1.99274064,1.09965342,-0.510839303,-1.75493479,-45.5170352,-6.74694450e05]
    return 322*(1.0+b[1]*d^(1.0/3.0)+b[2]*d^(2.0/3.0)+b[3]*d^(5.0/3.0)+b[4]*d^(16.0/3.0)+b[5]*d^(43.0/3.0)+b[6]*d^(110.0/3.0))
    
end
@inline function waterf0(delta,tau) #ideal helmholtz
    
#delta = rho/rhoc,
#tau = Tc/T
    
    n= [-8.3204464837497, 6.6832105275932, 3.00632,0.012436, 0.97315, 1.2795, 0.96956, 0.24873]
    gamma = [0, 0, 0, 1.28728967, 3.53734222, 7.74073708, 9.24437796,27.5075105]
    res = log(delta)+n[1]+n[2]*tau+n[3]*log(tau)
    
   for i = 4:8
     res=res+n[i]*log(-expm1(-gamma[i]*tau))
    end
    return res
end

@inline function waterfr(delta,tau) #residual helmholtz

#delta = rho/rhoc,
#tau = Tc/T

    nr1 = [0.12533547935523e-1, 0.78957634722828e1, -0.87803203303561e1,
    0.31802509345418, -0.26145533859358, -0.78199751687981e-2,
    0.88089493102134e-2]
    d1 = [1, 1, 1, 2, 2, 3, 4]
    t1 = [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1]
    
    nr2 = [-0.66856572307965, 0.20433810950965, -0.66212605039687e-4,
    -0.19232721156002, -0.25709043003438, 0.16074868486251,
    -0.4009282892587e-1, 0.39343422603254e-6, -0.75941377088144e-5,
    0.56250979351888e-3, -0.15608652257135e-4, 0.11537996422951e-8,
    .36582165144204e-6, -.13251180074668e-11, -.62639586912454e-9,
    -0.10793600908932, 0.17611491008752e-1, 0.22132295167546,
    -0.40247669763528, 0.58083399985759, 0.49969146990806e-2,
    -0.31358700712549e-1, -0.74315929710341, 0.47807329915480,
    0.20527940895948e-1, -0.13636435110343, 0.14180634400617e-1,
    0.83326504880713e-2, -0.29052336009585e-1, 0.38615085574206e-1,
    -0.20393486513704e-1, -0.16554050063734e-2, .19955571979541e-2,
    0.15870308324157e-3, -0.16388568342530e-4, 0.43613615723811e-1,
    0.34994005463765e-1, -0.76788197844621e-1, 0.22446277332006e-1,
    -0.62689710414685e-4, -0.55711118565645e-9, -0.19905718354408,
    0.31777497330738, -0.11841182425981]
    c2 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6,6, 6]
    d2 = [1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3,
    4, 4, 4, 5, 6, 6, 7, 9, 9, 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6]
    t2= [4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 
    10, 10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23,23, 10, 50, 44, 46, 50]
    gamma2 = fill(1,44)
    
    nr3 = [-0.31306260323435e2, 0.31546140237781e2, -0.25213154341695e4]
    d3 = fill(3,3)
    t3 = [0, 1, 4]
    alpha3 =  fill(20,3)
    beta3 = [150, 150, 250]
    gamma3 = [1.21, 1.21, 1.25]
    epsilon3 =  fill(1,3)
    
    nr4 =[-0.14874640856724, 0.31806110878444]
    a4 =[3.5, 3.5]
    b4 =[0.85, 0.95]
    B =[0.2, 0.2]
    C =[28, 32]
    D =[700, 800]
    A =[0.32,0.32]
    beta4 =[0.3, 0.3]

    
    res=0
    for i = 1:7
        res=res+(nr1[i]*delta^(d1[i])) * (tau^t1[i])
    end

    for i = 1:44
        res=res+(nr2[i]*delta^(d2[i])) * (tau^t2[i]) * exp(-delta^c2[i])
    end

    for i = 1:3
        res=res+(nr3[i]*delta^(d3[i])) * (tau^t3[i]) * exp(-alpha3[i]*abs2(delta-epsilon3[i])-beta3[i]*abs(tau-gamma3[i]))
    end
    delta1m2 = (delta-1)^2# 
    tau1m2 = (tau-1)^2
    
    for i = 1:2
        theta = (1-tau) + A[i]*delta1m2^(1/(2*beta4[i]))
        del = theta^2 + B[i]*delta1m2^a4[i]
        psi = exp(- C[i]*delta1m2 - D[i]*tau1m2)
        res = res+nr4[i]*del^b4[i]*psi
    end
        return res
end

waterf(rho,T) = (waterfr(rho/322.0,647.096/T)+waterf0(rho/322.0,647.096/T)) #total helmholtz function



#probable exportation of future helmholtz interfase for lavoisier
#Hemholtz(pvtxstate) = waterf(1/pvtxstate[2],pvtxstate[3])

function dwaterf(rho,T) #derivative of total helmholtz function
    d = [rho/322.0,647.096/T]
    f(d) = waterf(d[1],d[2])
    return ForwardDiff.gradient(f,d)
end 

function dwaterfr(rho,T) #derivative of residual helmholtz function
    d = [rho/322.0,647.096/T]
    f(d) = waterfr(d[1],d[2])
    return ForwardDiff.gradient(f,d)
end 

#pressure of water, using density (kg/m3) and Temperature(K), in Pa
pressure(rho,T) = rho*0.46151805*1000*T*(1+(rho/322.0)*dwaterfr(rho,T)[1])
